from tqdm import tqdm
import os
from os.path import join
import glob
import pickle
import numpy as np
import sys
import pandas as pd
#import click


def get_coverage(bam_dir, pileup_dir, barcode_f, reads_filter=-1,
                 maxBP=16571, base_qual=0, alignment_qual=0,
                 func_p="/home/isshamie/software/mito-genotyping/exampleProcessing/",
                 cellr_bc_f=None):
    if reads_filter == -1:
        reads_filter = maxBP

    print("Running get_coverage")
    if not os.path.exists(pileup_dir):
        os.mkdir(pileup_dir)
    CB_read_number = pickle.load(
        open(barcode_f, "rb"))


    CB_read_200 = set()
    for i in CB_read_number:
        if CB_read_number[i] >= reads_filter:
            CB_read_200.add(i)
    print(f"Number of cells: {len(CB_read_200)}")
    CB_read_200 = np.array(list(CB_read_200))

    for CB in tqdm(CB_read_200):
        bamfile = f"{join(bam_dir, 'CB_' + CB + '*.bam')}"
        if len(glob.glob(bamfile)) == 0:
            print(f"file not done yet {bamfile}")
            continue
        bamfile = glob.glob(bamfile)[0]
        outpre = join(pileup_dir,
                      os.path.basename(bamfile.replace(".bam", "")))
        sample = os.path.basename(bamfile).split("_")[-1]
        cmd = f"python {func_p}/01_pileup_counts.py {bamfile} {outpre} {maxBP} {base_qual} {sample} {alignment_qual}"
        #print(cmd)
        os.system(cmd)
    return


if __name__=="__main__":
    get_coverage(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))

#get_coverage(bam_dir="../data/raw/scBam",pileup_dir="../data/raw/scPileup", barcode_f="../data/scBam/barcode_data.p")