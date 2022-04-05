##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory
# Split_script code modified from https://github.com/pezmaster31/bamtools/issues/135
### Python 3.6.8
import click
import pysam
import os
from os.path import join
from tqdm import tqdm
import numpy as np
import glob
import pickle
import sys

# Remove CB
# cmd = 'samtools view data/aml035_post_transplant_possorted_genome_bam.MT.bam | grep "CB:Z:" >> data/aml035_post_transplant_possorted_genome_bam.MT_CB.bam'
# os.system(cmd)
# cmd = "samtools view -S -b aml035_post_transplant_possorted_genome_bam.MT_CB.bam > aml035_post_transplant_possorted_genome_bam.MT_CB.bam"
# os.system(cmd)


def split(in_f, out_d, to_overwrite=False, use_cb_suffix=True):
    """ Splits bam file into single-cell bam files with CB as barcode.

    Input file is a CB sorted bam file, but if not, then will sort.
    If already sorted, suffix needs to be .CB.bam.
    Otherwise, it will add that as a suffix to in_f.
    :param in_f:
    :param out_d:
    :param to_overwrite:
    :return:
    """
    # Sort on CB
    if use_cb_suffix:
        if ".CB.bam" in in_f:
            unsplit_file = in_f
        else:
            unsplit_file = in_f.replace(".bam", "") + ".CB.bam"
    else:
        unsplit_file = in_f
    print(in_f, unsplit_file)

    if to_overwrite or (not os.path.exists(unsplit_file)):
        cmd = f"samtools sort -t CB {in_f} > {unsplit_file}"
        print(cmd)
        os.system(cmd)
    else:
        print("Already sorted. Not sorting")
    ### Input varibles to set
    # file to split on
    #unsplit_file = "data/aml035_post_transplant_possorted_genome_bam.MT.CB.bam"
    # where to place output files

    out_dir = out_d
    #out_dir = os.path.join(os.path.dirname(unsplit_file), out_d)
    if not os.path.exists(out_dir):
        print(f"creating directory {out_dir}")
        os.mkdir(out_dir)

    # variable to hold barcode index
    CB_hold = 'unset'
    itr = 0
    # read in upsplit file and loop reads by line
    count = 0
    count_perm = 0
    v = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(unsplit_file, "rb")
    pysam.set_verbosity(v)
    for read in tqdm(samfile.fetch(until_eof=True)):
        # barcode itr for current read
        try:
            #if "CB" in read.tags():
            CB_itr = read.get_tag('CB')
            # if change in barcode or first line; open new file  
            if (CB_itr!=CB_hold or itr==0):
                # close previous split file, only if not first read in file
                if (itr!=0):
                    try:
                        split_file.close()
                    except NameError:
                        print("split file not here yet")
                        continue
                CB_hold = CB_itr
                itr+=1
                #print(join(out_dir,"CB_{}.bam".format(CB_itr)))
                pysam.set_verbosity(v)
                try:
                    split_file = pysam.AlignmentFile(join(out_dir,"CB_{}.bam".format(CB_itr)), "wb", template=samfile)
                except PermissionError:
                    print("CB permission error")
                    count_perm += 1
                    continue

            # write read with same barcode to file
            split_file.write(read)
        except KeyError:
            count += 1
    split_file.close()
    samfile.close()

    print("no CBs: ", count)
    print("no CBs permission", count_perm)
    return


#@click.option("--rg", default=False) #in_f, out_d, to_overwrite=False, use_cb_suffix=True
#@click.option("--rdict", default=None) #in_f, out_d, to_overwrite=False, use_cb_suffix=True
@click.command()
@click.argument("in_f", type=click.Path(exists=True))
@click.argument("out_d", type=click.Path())
def main(in_f, out_d): #, rg, rdict):
    # if rdict is not None:
    #     rdict = rdict.split(";")
    #     rdict = {x.split(",")[0]:x.split(",")[1] for x in rdict}
    split(in_f, out_d) #, rg, rg_dict=rdict) #sys.argv[1], sys.argv[2], to_overwrite=False)


if __name__ == "__main__":
    main()

#split("../data/raw/14515X1.MT.bam", "scBam")
#get_coverage(bam_dir="../data/raw/scBam",pileup_dir="../data/raw/scPileup", barcode_f="../data/scBam/barcode_data.p")
