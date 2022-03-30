import click
from tqdm import tqdm
import os
from os.path import join, basename
import glob
import pickle
import numpy as np
import sys
from src.somvars.pileup_variants import run_pileup
from numpanpar import parallel_ar as parar


def parallel_run_pileup(bam_files, pileup_dir,
                        base_qual, alignment_qual):
    for bamfile in tqdm(bam_files):
        CB = basename(bamfile).split(".bam")[0]
        #print("CB", CB)
        #bamfile = f"{join(bamdir, 'CB_' + CB + '*.bam')}"
        if len(glob.glob(bamfile)) == 0:
            print(f"file not done yet {bamfile}")
            continue
        #bamfile = glob.glob(bamfile)
        #if len(bamfile) > 1:
        #    print(f"More than one cell matches {bamfile}")
        #    continue
        #bamfile = bamfile[0]
        outpre = join(pileup_dir, CB)
        run_pileup(bamfile, outpre, base_qual, alignment_qual)
    return


def get_coverage(bamdir, outdir, barcodes_f=None, base_qual=0,
                 alignment_qual=0, num_proc=4):
    """
    :param bamdir: the split single cell files here
    :param outdir:
    :param barcodes:
    :param base_qual:
    :param alignment_qual:
    :return:
    """
    from glob import glob
    print("Running get_coverage")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    bam_files = np.array(glob(bamdir+"/*CB*.bam"))
    if barcodes_f is not None:
        with open(barcodes_f, 'w') as f:
            barcodes = f.readlines()
        barcodes = [x.strip() for x in barcodes]

        new_bam_files = []
        for b in barcodes:
            if join(bamdir,f"CB_{b}.bam") in bam_files:
                new_bam_files.append(b)
        if len(new_bam_files)==0:
            raise ValueError("Barcodes are not matching")
        bam_files = new_bam_files

    print('bam_files', bam_files[:10])
    parar(bam_files, parallel_run_pileup,
          func_args=(outdir, base_qual, alignment_qual),
          num_processes=num_proc)
    return

@click.command()
@click.argument("bamdir", type=click.Path(exists=True))
@click.argument("outdir", type=click.Path())
@click.option("--barcodes", default=None)
@click.option("--bq", default=20)
@click.option("--aq", default=0)
@click.option("--nproc", default=8)
def main(bamdir, outdir, barcodes, bq, aq, nproc):
    get_coverage(bamdir, outdir, barcodes, bq, aq, nproc)
    return


if __name__ == "__main__":
    main()
    #get_coverage(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))

