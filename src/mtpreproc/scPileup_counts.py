from tqdm import tqdm
import os
from os.path import join
import glob
import pickle
import numpy as np
import sys
from src.external.pileup_counts import run_pileup
from numpanpar import parallel_ar as parar


def parallel_run_pileup(CB_read_filt, bam_dir, pileup_dir, maxBP,
                        base_qual, alignment_qual):

    for CB in tqdm(CB_read_filt):
        bamfile = f"{join(bam_dir, 'CB_' + CB + '*.bam')}"
        if len(glob.glob(bamfile)) == 0:
            print(f"file not done yet {bamfile}")
            continue
        bamfile = glob.glob(bamfile)
        if len(bamfile) > 1:
            print(f"More than one cell matches {bamfile}")
            continue
        bamfile = bamfile[0]
        outpre = join(pileup_dir,
                      os.path.basename(bamfile.replace(".bam", "")))
        sample = os.path.basename(bamfile).split("_")[-1]
        print('sample', sample)
        run_pileup(bamfile, outpre, maxBP, base_qual, sample,
                   alignment_qual, use_strands=True)
    return


def get_coverage(bam_dir, pileup_dir, barcode_f, reads_filter=-1, base_qual=0,
                 maxBP=16571, alignment_qual=0, num_proc=4):
    """

    :param bam_dir: the split single cell files here
    :param pileup_dir:
    :param barcode_f:
    :param reads_filter:
    :param base_qual:
    :param maxBP:
    :param alignment_qual:
    :return:
    """

    if reads_filter < 0:
        reads_filter = 0

    print("Running get_coverage")
    if not os.path.exists(pileup_dir):
        os.mkdir(pileup_dir)

    [_, CB_read_number, _, _,_, _] = pickle.load(open(barcode_f, "rb"))

    CB_read_filt = set()
    for i in CB_read_number:
        if CB_read_number[i] >= reads_filter:
            CB_read_filt.add(i)
    print(f"Number of cells: {len(CB_read_filt)}")
    CB_read_filt = np.array(list(CB_read_filt))

    parar(CB_read_filt, parallel_run_pileup,
          func_args=(bam_dir, pileup_dir, maxBP, base_qual, alignment_qual),
          num_processes=num_proc)
    #
    # for CB in tqdm(CB_read_filt):
    #     bamfile = f"{join(bam_dir, 'CB_' + CB + '*.bam')}"
    #     if len(glob.glob(bamfile)) == 0:
    #         print(f"file not done yet {bamfile}")
    #         continue
    #     bamfile = glob.glob(bamfile)[0]
    #     outpre = join(pileup_dir,
    #                   os.path.basename(bamfile.replace(".bam", "")))
    #     sample = os.path.basename(bamfile).split("_")[-1]
    #     run_pileup(bamfile, outpre, maxBP, base_qual, sample, alignment_qual, use_strands=True)
    return


if __name__=="__main__":
    get_coverage(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))

