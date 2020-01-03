##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### Python 3.6.8
import pysam
import os

# Remove CB
# cmd = 'samtools view data/aml035_post_transplant_possorted_genome_bam.MT.bam | grep "CB:Z:" >> data/aml035_post_transplant_possorted_genome_bam.MT_CB.bam'
# os.system(cmd)
# cmd = "samtools view -S -b aml035_post_transplant_possorted_genome_bam.MT_CB.bam > aml035_post_transplant_possorted_genome_bam.MT_CB.bam"
# os.system(cmd)

# Sort on CB
cmd = "samtools sort -t CB data/aml035_post_transplant_possorted_genome_bam.MT.bam > data/aml035_post_transplant_possorted_genome_bam.MT.CB.bam"


### Input varibles to set
# file to split on
unsplit_file = "data/aml035_post_transplant_possorted_genome_bam.MT.CB.bam"
# where to place output files
out_dir = "data/scBam"

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
count = 0
samfile = pysam.AlignmentFile( unsplit_file, "rb")
for read in samfile.fetch(until_eof=True):
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
            split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format( itr), "wb", template=samfile)
                # write read with same barcode to file
        split_file.write( read)
    except KeyError:
        print("no CB")
        count += 1
    

split_file.close()
samfile.close()

print("no CBs: ", count)

