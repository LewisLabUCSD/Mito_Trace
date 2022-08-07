
"""
based on https://neuroscience.telenczuk.pl/?p=331
"""

from svgutils.compose import *
import svgutils.transform as sg
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})



class Svg(object):
    "svg files with a data object (the svg), width, height and coordinates"

    def __init__(self, data, dim, coords, name, global_scale=0.5):
        self.data = data
        self.width = dim[0]*global_scale
        self.height = dim[1]*global_scale
        self.x = coords[0]
        self.y = coords[1]
        self.global_scale = global_scale
        self.name = name

    def scale_width_to_reference(self, reference_width):
        """Proportionally scale the image to a given width."""
        scalings_factor = reference_width / self.width
        self.data.moveto(0, 0, scale=scalings_factor)
        self.width = self.width * scalings_factor
        self.height = self.height * scalings_factor

    def scale_height_to_reference(self, reference_height):
        """Proportionally scale the image to a given height."""
        scalings_factor = reference_height / self.height
        self.data.moveto(0, 0)#, scale=scalings_factor)
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
                paramdict = {e.split('="')[0]: e.split('="')[1] for e in line[5:-2].split('" ')}
                break
    return int(paramdict["width"].replace('pt', '')), int(paramdict["height"].replace('pt', ''))



def rescale_height(svgs, curr_keys, ref_key):
    """Change the dimensions of the images to the desired combinations."""
    for key in curr_keys:
        svgs[key].scale_height_to_reference(svgs[ref_key].height)
    assert svgs[ref_key].height == svgs[curr_keys[0]].height == svgs[curr_keys[-1]].height


def change_positions(df, svgs, widths, heights):
    """Move the images to the desired positions."""
    for i_row, row in enumerate(df.index):
        for i_col, col in enumerate(df.columns):
            curr = df.loc[row, col].split(".svg")[0]
            svgs[curr].move(sum(widths[:i_col]), sum(heights[:i_row]))


def letter_annotations(svgs, df=None, only_grid_titles=True, figure_names=None):
    if only_grid_titles:
        all_texts = []
        for ind in df.index:
            value = svgs[df.loc[ind].iloc[0].split(".svg")[0]]
            all_texts.append(sg.TextElement(value.x, value.y, ind, size=15, weight="bold"))
        for col in df.columns:
            nm = figure_names[col]
            value = svgs[df[col].iloc[0].split(".svg")[0]]
            all_texts.append(sg.TextElement(value.x, value.y-10, nm, size=15, weight="bold"))
        return all_texts
    """Add letters based on the location of the images."""
    return [sg.TextElement(value.x, value.y, value.name, size=15, weight="bold")
            for key, value in svgs.items()]
#     return [sg.TextElement(value.x + 10, value.y + 15, value.name, size=15, weight="bold")
#             for key, value in svgs.items()]


def files_to_svg_dict(df, scale=0.5):
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
            out[key] = Svg(
                data=sg.fromfile(s).getroot().scale(scale),
                dim=get_size(s),
                coords=(0, 0),
                global_scale=scale,
                name=name)
    return out
#     return {
#         s.split('.svg')[0]: Svg(
#             data=sg.fromfile(s).getroot(),
#             dim=get_size(s),
#             coords=(0, 0))
#         for s in files}




# def main(out_f):
#     #svgs = files_to_svg_dict(["A.svg", "B.svg", "C.svg", "D.svg", "E.svg", "F.svg"])
#     rescale_height(svgs)
#     change_positions(svgs)
#     full_width = sum([svgs[i].width for i in ['A', 'E']])
#     full_height = sum([svgs[i].height for i in ['A', 'B', 'C', 'D']])
#     fig = sg.SVGFigure(full_width, full_height)
#     text = letter_annotations(svgs)
#     fig.append([s.data for s in svgs.values()])
#     fig.append(text)
#     fig.save(out_f)



def main_df(df, out_f, scale=0.25, figure_names=None):
    svgs = files_to_svg_dict(df, scale=scale)#df.values().flatten())

    for don, val in df.iterrows():
        curr_init_files = []
        curr_ref = val.iloc[0].split('.svg')[0] #svgs[val.iloc[0].split('.svg')[0]]
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

    row_max_heights = [] #should be the same in theory
    for ind in df.index:
        curr_row_max = 0
        for col in df.columns:
            tmp_height = svgs[df.loc[ind, col].split(".svg")[0]].height
            if tmp_height > curr_row_max:
                curr_row_max = tmp_height
        row_max_heights.append(curr_row_max)

    change_positions(df, svgs, col_max_widths, row_max_heights)

    full_width = sum(col_max_widths) #[svgs[i].width for i in ['A', 'E']])
    full_height = sum(row_max_heights) #[svgs[i].height for i in ['A', 'B', 'C', 'D']])
    fig = sg.SVGFigure(full_width, full_height)
    text = letter_annotations(svgs, df=df, figure_names=figure_names)
    fig.append([s.data for s in svgs.values()])
    fig.append(text)
    fig.save(out_f)

    cmd = f'inkscape --file={out_f} -d 300 --export-area-drawing --without-gui --export-png={out_f.replace(".svg", ".png")}' # --export-area-drawing --without-gui --export-pdf=output.pdf'
    print(cmd)
    os.system(cmd)
    cmd = f'inkscape --file={out_f} --export-area-drawing --without-gui --export-pdf={out_f.replace("svg", "pdf")} '
    print(cmd)
    os.system(cmd)




