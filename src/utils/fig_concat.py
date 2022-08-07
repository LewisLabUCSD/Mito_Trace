"""
based on https://neuroscience.telenczuk.pl/?p=331
"""


from glob import glob
import pickle
import mplh.cluster_help as ch
import mplh.fig_utils as fu
import numpy as np
import os
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

import pandas as pd
import seaborn as sns
np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})

from svgutils.compose import *
import svgutils.transform as sg
import sys
from os.path import basename

class Svg(object):
    "svg files with a data object (the svg), width, height and coordinates"

    def __init__(self, data, dim, coords, name, global_scale=0.5):
        self.data = data
        self.width = dim[0] * global_scale
        self.height = dim[1] * global_scale
        self.x = coords[0]
        self.y = coords[1]
        self.global_scale = global_scale
        self.name = name

    def scale_width_to_reference(self, reference_width):
        """Proportionally scale the image to a given width."""
        scalings_factor = reference_width / self.width
        self.data.moveto(0, 0 )  # , scale=scalings_factor)
        self.data.scale(scalings_factor)

        self.width = self.width * scalings_factor
        self.height = self.height * scalings_factor

    def scale_height_to_reference(self, reference_height):
        """Proportionally scale the image to a given height."""
        scalings_factor = reference_height / self.height
        self.data.moveto(0, 0)  # , scale=scalings_factor)
        self.data.scale(scalings_factor)
        self.width = self.width * scalings_factor
        self.height = self.height * scalings_factor

    def scale_by_factor(self, scalings_factor):
        """Proportionally scale image by a scaling factor."""
        self.data.moveto(0, 0, scale=scalings_factor)
        self.width = self.width * scalings_factor
        self.height = self.height * scalings_factor

    def move(self, x, y):
        """Move the coordinates of an image."""
        self.data.moveto(x, y)
        self.x = x
        self.y = y


def get_size(svg_file):
    """Naively parse the svg text file to get the width and height."""
    with open(svg_file) as svg:
        for line in svg:
            if line.startswith('<svg'):
                try:
                    paramdict = {e.split('="')[0]: e.split('="')[1] for e in
                                 line[5:-2].split('" ')}
                except IndexError:
                    paramdict = {e.split("='")[0]: e.split("='")[1] for e in
                                 line[5:-2].split("' ")}
                break
    print(paramdict["width"])
    print("height")
    print(paramdict["height"])
    return int(np.ceil(float(
        paramdict["width"].replace('pt', '').replace('"', '')))), int(
        np.ceil(float(
            paramdict["height"].replace('pt', '').replace('"', ''))))


def rescale_height(svgs, curr_keys, ref_key):
    """Change the dimensions of the images to the desired combinations."""
    for key in curr_keys:
        svgs[key].scale_height_to_reference(svgs[ref_key].height)  # assert svgs[ref_key].height == svgs[curr_keys[0]].height == svgs[curr_keys[-1]].height


def rescale_width(svgs, curr_keys, ref_key):
    """Change the dimensions of the images to the desired combinations."""
    for key in curr_keys:
        svgs[key].scale_width_to_reference(svgs[
                                               ref_key].width)  # assert svgs[ref_key].height == svgs[curr_keys[0]].height == svgs[curr_keys[-1]].height


def change_positions(df, svgs, widths, heights):
    """Move the images to the desired positions."""
    for i_row, row in enumerate(df.index):
        for i_col, col in enumerate(df.columns):
            curr = df.loc[row, col].split(".svg")[0]
            svgs[curr].move(sum(widths[:i_col]), sum(heights[:i_row]))


def letter_annotations(svgs, df=None, only_grid_titles=True,
                       figure_names=None, transposed=False):
    if only_grid_titles:
        all_texts = []
        if transposed:  # make donors the columns and figure_names the rows
            for ind in df.index:
                nm = basename(figure_names[ind])
                value = svgs[df.loc[ind].iloc[0].split(".svg")[0]]
                all_texts.append(
                    sg.TextElement(value.x, value.y - 10, nm, size=15,
                                   weight="bold"))
            for col in df.columns:
                value = svgs[df[col].iloc[0].split(".svg")[0]]
                all_texts.append(
                    sg.TextElement(value.x, value.y - 10, col, size=15,
                                   weight="bold"))

        else:
            for ind in df.index:
                value = svgs[df.loc[ind].iloc[0].split(".svg")[0]]
                all_texts.append(
                    sg.TextElement(value.x, value.y, ind, size=15,
                                   weight="bold"))
            for col in df.columns:
                nm = basename(figure_names[col])
                value = svgs[df[col].iloc[0].split(".svg")[0]]
                all_texts.append(
                    sg.TextElement(value.x, value.y - 10, nm, size=15,
                                   weight="bold"))
        return all_texts
    """Add letters based on the location of the images."""
    return [sg.TextElement(value.x, value.y, value.name, size=15,
                           weight="bold") for key, value in
            svgs.items()]


#     return [sg.TextElement(value.x + 10, value.y + 15, value.name, size=15, weight="bold")
#             for key, value in svgs.items()]


def files_to_svg_dict(df, scale=0.5, verbose=False):
    """Convert a list of images to a dictionary.
    Mapping the image basename to the Svg class instance,
    setting the dimensions based on sizes and coordinates (0,0) by default
    """
    out = {}
    for row in df.index:
        for col in df.columns:
            s = df.loc[row, col]
            key = df.loc[row, col].split(".svg")[0]
            name = f"{row} {col.replace('dendro', '').replace('.svg', '')}"
            if verbose:
                print(s)
            out[key] = Svg(data=sg.fromfile(s).getroot().scale(scale),
                           dim=get_size(s), coords=(0, 0),
                           global_scale=scale, name=name)
    return out


def main_df(df, out_f, scale=0.25, figure_names=None, to_pdf=False,
            verbose=False, add_titles=True):
    svgs = files_to_svg_dict(df, scale=scale, verbose=verbose)  # df.values().flatten())

    for don, val in df.iterrows():
        print(don)
        curr_ref = val.iloc[0].split('.svg')[
            0]  # svgs[val.iloc[0].split('.svg')[0]]
        curr_keys = [x.split('.svg')[0] for x in val.values]
        rescale_height(svgs, curr_keys, curr_ref)

    col_max_widths = []

    for col in df.columns:
        curr_col_max = 0
        for ind in df.index:
            tmp_width = svgs[df.loc[ind, col].split(".svg")[0]].width
            if tmp_width > curr_col_max:
                curr_col_max = tmp_width
        col_max_widths.append(curr_col_max)

    row_max_heights = []  # should be the same in theory
    for ind in df.index:
        curr_row_max = 0
        for col in df.columns:
            tmp_height = svgs[df.loc[ind, col].split(".svg")[0]].height
            if tmp_height > curr_row_max:
                curr_row_max = tmp_height
        row_max_heights.append(curr_row_max)

    change_positions(df, svgs, col_max_widths, row_max_heights)

    full_width = sum(
        col_max_widths)  # [svgs[i].width for i in ['A', 'E']])
    full_height = sum(
        row_max_heights)  # [svgs[i].height for i in ['A', 'B', 'C', 'D']])
    fig = sg.SVGFigure(full_width + 20, full_height + 20)
    if add_titles:
        text = letter_annotations(svgs, df=df, figure_names=figure_names)
        fig.append(text)
    fig.append([s.data for s in svgs.values()])

    fig.save(out_f)

    cmd = f'inkscape --file={out_f} -d 300 --export-area-drawing --without-gui --export-png={out_f.replace(".svg", ".png")}'  # --export-area-drawing --without-gui --export-pdf=output.pdf'
    print(cmd)
    os.system(cmd)

    if to_pdf:
        cmd = f'inkscape --file={out_f} --export-area-drawing --without-gui --export-pdf={out_f.replace("svg", "pdf")} '
        print(cmd)
        os.system(cmd)


def main_df_transpose(df, out_f, scale=0.25, figure_names=None,
                      to_pdf=False, to_transpose=True, add_titles=True):
    svgs = files_to_svg_dict(df, scale=scale)  # df.values().flatten())

    df = df.transpose()
    for don, val in df.iterrows():
        curr_init_files = []
        curr_ref = val.iloc[0].split('.svg')[
            0]  # svgs[val.iloc[0].split('.svg')[0]]
        curr_keys = [x.split('.svg')[0] for x in val.values]
        rescale_width(svgs, curr_keys, curr_ref)

    col_max_widths = []
    for col in df.columns:
        curr_col_max = 0
        for ind in df.index:
            tmp_width = svgs[df.loc[ind, col].split(".svg")[0]].width
            if tmp_width > curr_col_max:
                curr_col_max = tmp_width
        col_max_widths.append(curr_col_max)

    row_max_heights = []  # should be the same in theory
    for ind in df.index:
        curr_row_max = 0
        for col in df.columns:
            tmp_height = svgs[df.loc[ind, col].split(".svg")[0]].height
            if tmp_height > curr_row_max:
                curr_row_max = tmp_height
        row_max_heights.append(curr_row_max)

    change_positions(df, svgs, col_max_widths, row_max_heights)

    full_width = sum(
        col_max_widths)  # [svgs[i].width for i in ['A', 'E']])
    full_height = sum(
        row_max_heights)  # [svgs[i].height for i in ['A', 'B', 'C', 'D']])
    fig = sg.SVGFigure(full_width + 20, full_height + 20)
    if add_titles:
        text = letter_annotations(svgs, df=df, figure_names=figure_names,
                                  transposed=True)
        fig.append(text)
    fig.append([s.data for s in svgs.values()])

    fig.save(out_f)

    cmd = f'inkscape --file={out_f} -d 300 --export-area-drawing --without-gui --export-png={out_f.replace(".svg", ".png")}'  # --export-area-drawing --without-gui --export-pdf=output.pdf'
    print(cmd)
    os.system(cmd)

    if to_pdf:
        cmd = f'inkscape --file={out_f} --export-area-drawing --without-gui --export-pdf={out_f.replace("svg", "pdf")} '
        print(cmd)
        os.system(cmd)

