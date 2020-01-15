import glob
import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import pysam
import time
from collections import defaultdict
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def extract_barcode_info(bam_f, out_f,rm_slash=False):
    samfile = pysam.AlignmentFile(bam_f, "rb")
    barcodes = set()
    corrected_barcodes = set()
    barcode_pairs = set()
    CR_read_number = defaultdict(int)
    CB_read_number = defaultdict(int)
    BC_read_number = defaultdict(int)

    for read in tqdm(samfile.fetch('MT')):
        d = dict(read.tags)
        if "CR" in d:
            barcodes.add(d["CR"])
            CR_read_number[d["CR"]] += 1
            if "CB" in d:
                if rm_slash:
                    CB = d["CB"].split('-')[0]
                else:
                    CB = d["CB"]
                barcode_pairs.add((CB, d["CR"]))
                corrected_barcodes.add(CB)
                CB_read_number[CB] += 1
            else:
                barcode_pairs.add(('-', d["CR"]))
        else:
            if "CB" in d:
                if rm_slash:
                    CB = d["CB"].split['-'][0]
                else:
                    CB = d["CB"]
                #pair_barcodes.add((CB, '-'))
                corrected_barcodes.add(CB)
                CB_read_number[CB] += 1
        if "BC" in d:
            BC_read_number[d["BC"]] += 1 
    samfile.close()
    pickle.dump([CR_read_number,CB_read_number,BC_read_number, barcodes, corrected_barcodes, barcode_pairs], file=open(out_f, "wb"))
    return CR_read_number,CB_read_number,BC_read_number, barcodes, corrected_barcodes, barcode_pairs


def main():
    bam_f = sys.argv[1]
    out_f = sys.argv[2]
    extract_barcode_info(bam_f, out_f, rm_slash=False)
    return


if __name__== '__main__':
    main()
