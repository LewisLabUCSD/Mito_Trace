import os
import click
import pandas as pd
import glob


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
    # Combine all files

    print("A")
    os.system(f'cat {directory}/*.A.txt | gzip > "{samplename}_all.A.txt.gz"')
    os.system(
        f'cat {directory}/*.A.minus.txt | gzip > "{samplename}_all.A.minus.txt.gz"')
    print("C")
    os.system(
        f'cat {directory}/*.C.txt | gzip > "{samplename}_all.C.txt.gz"')
    os.system(
        f'cat {directory}/*.C.minus.txt | gzip > "{samplename}_all.C.minus.txt.gz"')

    print("G")
    os.system(
        f'cat {directory}/*.G.txt | gzip > "{samplename}_all.G.txt.gz"')
    os.system(
        f'cat {directory}/*.G.minus.txt | gzip > "{samplename}_all.G.minus.txt.gz"')

    print("T")
    os.system(
        f'cat {directory}/*.T.txt | gzip > "{samplename}_all.T.txt.gz"')
    os.system(
        f'cat {directory}/*.T.minus.txt | gzip > "{samplename}_all.T.minus.txt.gz"')

    print("Coverage")
    os.system(f'cat {directory}/*.coverage.txt | gzip > "{samplename}_all.coverage.txt.gz"')
    os.system(
        f'cat {directory}/*.coverage.minus.txt | gzip > "{samplename}_all.coverage.minus.txt.gz"')


    ## Merge the two
    print("Merging the Forward and Reverse strand into one file")
    for nt in ["A", "C", "G", "T"]:
        print(nt)
        curr_fw = pd.read_csv(f"{samplename}_all.{nt}.txt.gz")
        curr_rev = pd.read_csv(f"{samplename}_all.{nt}.minus.txt.gz")
        curr_fw.columns = ["Position", "CB", "Count Fw", "BQ Fw"]
        curr_rev.columns = ["Position", "CB", "Count Rev", "BQ Rev"]
        curr = pd.merge(curr_fw, curr_rev, on=["Position","CB"])
        curr.to_csv(f"{samplename}_all.{nt}.strands.txt.gz",compression='gzip', index=None)
    print('Coverage')
    curr_fw = pd.read_csv(f"{samplename}_all.coverage.txt.gz")
    curr_rev = pd.read_csv(f"{samplename}_all.coverage.minus.txt.gz")
    curr_fw.columns = ["Position", "CB", "Count Fw"]
    curr_rev.columns = ["Position", "CB", "Count Rev"]

    curr = pd.merge(curr_fw, curr_rev, on=["Position", "CB"], how="outer").fillna(0)
    curr = curr.astype({"Count Fw": int, "Count Rev": int})
    curr["Coverage"] = (curr["Count Fw"] + curr["Count Rev"])/2
    curr = curr[["Position", "CB", "Coverage", "Count Fw", "Count Rev"]]
    curr.to_csv(f"{samplename}_all.coverage.strands.txt.gz", compression='gzip', index=None)
    return


if __name__ == '__main__':
    main()

