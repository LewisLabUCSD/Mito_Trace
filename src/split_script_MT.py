##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory
# Split_script code modified from https://github.com/pezmaster31/bamtools/issues/135
### Python 3.6.8
import pysam
import os
from os.path import join
from tqdm import tqdm
import numpy as np
import glob
import pickle
# Remove CB
# cmd = 'samtools view data/aml035_post_transplant_possorted_genome_bam.MT.bam | grep "CB:Z:" >> data/aml035_post_transplant_possorted_genome_bam.MT_CB.bam'
# os.system(cmd)
# cmd = "samtools view -S -b aml035_post_transplant_possorted_genome_bam.MT_CB.bam > aml035_post_transplant_possorted_genome_bam.MT_CB.bam"
# os.system(cmd)


def split(in_f, out_d,to_overwrite=False):
    
    # Sort on CB
    unsplit_file = in_f.replace(".bam","") + ".CB.bam"
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
    
    out_dir = os.path.join(os.path.dirname(unsplit_file), out_d)
    if not os.path.exists(out_dir):
        print(f"creating directory {out_dir}")
        os.mkdir(out_dir)
    
        
    # variable to hold barcode index
    CB_hold = 'unset'
    itr = 0
    # read in upsplit file and loop reads by line
    count = 0
    samfile = pysam.AlignmentFile( unsplit_file, "rb")
    for read in tqdm(samfile.fetch(until_eof=True)):
        # barcode itr for current read
        try:
            #if "CB" in read.tags():
            CB_itr = read.get_tag('CB')
            # if change in barcode or first line; open new file  
            if( CB_itr!=CB_hold or itr==0):
                # close previous split file, only if not first read in file
                if( itr!=0):
                    split_file.close()
                CB_hold = CB_itr
                itr+=1
                #print(join(out_dir,"CB_{}.bam".format(CB_itr)))
                split_file = pysam.AlignmentFile(join(out_dir,"CB_{}.bam".format(CB_itr)), "wb", template=samfile)
            # write read with same barcode to file
            split_file.write(read)
        except KeyError:
            count += 1


    split_file.close()
    samfile.close()

    print("no CBs: ", count)
    return


def get_coverage(bam_dir, pileup_dir, barcode_f, reads_filter=200, maxBP=16571, base_qual=0, alignment_qual=0, func_p="/home/isshamie/software/mito-genotyping/exampleProcessing/"):
    print("Running get_coverage")
    if not os.path.exists(pileup_dir):
        os.mkdir(pileup_dir)
    [CR_read_number,CB_read_number,BC_read_number, barcodes, corrected_barcodes, barcode_pairs] = pickle.load(open(barcode_f,"rb"))
    CB_read_200 = set()
    for i in CB_read_number:
        if CB_read_number[i] >= reads_filter:
            CB_read_200.add(i)
    print(f"Number of cells: {len(CB_read_200)}")
    CB_read_200 = np.array(list(CB_read_200))
    
    for CB in tqdm(CB_read_200):
        bamfile = f"{join(bam_dir,'CB_' + CB + '-*.bam')}"
        if len(glob.glob(bamfile)) == 0:
            print(f"file not done yet {bamfile}")
            continue
        bamfile = glob.glob(bamfile)[0]
        outpre = join(pileup_dir, os.path.basename(bamfile.replace(".bam","")))
        sample = os.path.basename(bamfile).split("_")[-1]  
        cmd = f"python {func_p}/01_pileup_counts.py {bamfile} {outpre} {maxBP} {base_qual} {sample} {alignment_qual}"
        print(cmd)
        os.system(cmd)
    return
    

#split("../data/raw/14515X1.MT.bam", "scBam")

#get_coverage(bam_dir="../data/raw/scBam",pileup_dir="../data/raw/scPileup", barcode_f="../data/scBam/barcode_data.p")
