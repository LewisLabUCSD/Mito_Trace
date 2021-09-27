import os
import click
import pandas as pd
import subprocess as subp
import glob
from os.path import join
from subprocess import *



def concat_files(directory, samplename, nt):
    cmd = f"find {directory} -type f -name *.{nt}.txt -exec cat {{}} \; > {samplename}_all.{nt}.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    cmd = f"find {directory} -type f -name *.{nt}.minus.txt -exec cat {{}} \; > {samplename}_all.{nt}.minus.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    return


#Samplename is full directory path plus the name: e.g. 'all_coverage/A'
@click.command()
@click.argument('directory', type=click.Path(exists=True))
@click.argument('samplename', type=click.STRING)
def main(directory, samplename):
    """This script takes all the A,C,G,T pileups for single cells and merges them

    :param directory: The directory where the single-cell bam files ares
    :param samplename: Samplename is full directory path plus the name: e.g. 'all_coverage/A'
    :return:
    """
    print(directory)
    print(samplename)

    ## Merge the two
    print("Merging the Forward and Reverse strand into one file")
    for nt in ["A", "C", "G", "T"]:
        print(nt)
        print(join(directory, f"{samplename}_all.{nt}.txt"))
        curr_fw = pd.read_csv(join(directory, f"{samplename}_all.{nt}.txt"), header=None)
        curr_rev = pd.read_csv(join(directory, f"{samplename}_all.{nt}.minus.txt"), header=None)
        curr_fw.columns = ["Position", "CB", "Count Fw", "BQ Fw"]
        curr_rev.columns = ["Position", "CB", "Count Rev", "BQ Rev"]
        curr = pd.merge(curr_fw, curr_rev, on=["Position","CB"], how="outer")
        curr = curr.sort_values(["Position", "CB"])
        curr = curr[["Position", "CB", "Count Fw", "BQ Fw", "Count Rev", "BQ Rev"]]
        curr.to_csv(join(directory, f"{samplename}_all.{nt}.strands.txt"), index=None)

    print('Coverage')
    curr_fw = pd.read_csv(join(directory,f"{samplename}_all.coverage.txt"), header=None, error_bad_lines=False)
    curr_rev = pd.read_csv(join(directory,f"{samplename}_all.coverage.minus.txt"), header=None, error_bad_lines=False)
    curr_fw.columns = ["Position", "CB", "Count Fw"]
    curr_rev.columns = ["Position", "CB", "Count Rev"]

    curr = pd.merge(curr_fw, curr_rev, on=["Position", "CB"], how="outer").fillna(0)
    curr = curr.astype({"Count Fw": int, "Count Rev": int})
    curr["Coverage"] = (curr["Count Fw"] + curr["Count Rev"])
    curr = curr.sort_values(["Position", "CB"])
    curr = curr[["Position", "CB", "Coverage", "Count Fw", "Count Rev"]]
    curr.to_csv(join(directory, f"{samplename}_all.coverage.strands.txt"), index=None)
    return


if __name__ == '__main__':
    main()

