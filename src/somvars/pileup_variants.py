#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
from os.path import basename
import pandas as pd
import pysam
pysam.set_verbosity(0)
from os.path import join
# bamfile = sys.argv[1] # Filepath to raw bamfile
# outpre = sys.argv[2] # Prefix / basename for raw file output
# maxBP = sys.argv[3] # Maximum length of mtDNA genome
# base_qual = float(sys.argv[4]) # Minimum base quality to be considered for pileup
# sample = sys.argv[5] # Sample name to be considered for downstream analyses
# alignment_quality = float(sys.argv[6]) # minimum alignment quality required to be considered


# Export Functions
def writeSparseMatrix(mid, vec, outpre, maxBP, sample):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")

def writeSparseMatrix2(mid, vec1, vec2, outpre, maxBP, sample):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")


def run_pileup(bamfile, outpre, base_qual, alignment_quality): # vcf, vcf_chr_col="#CHROM", vcf_pos_col="POS", pad_bp=2):
	# BAQ
	# initialize with a pseudo count to avoid dividing by zero
	# countsA = [0.00000001] * n
	# countsC = [0.00000001] * n
	# countsG = [0.00000001] * n
	# countsT = [0.00000001] * n
	# qualA = [0.0] * n
	# qualC = [0.0] * n
	# qualG = [0.0] * n
	# qualT = [0.0] * n

	variants = {}
	#pileups_d = {"nt": [], "chr": [], "counts": [], "BQ": []}
	from collections import defaultdict

	cell = basename(bamfile).split(".bam")[0]
	if "CB_" in cell:
		cell = cell.split("CB_")[1]

	pileups_d = defaultdict(lambda: 0)
	bq_d = defaultdict(lambda: 0)

	try:
		bam2 = pysam.AlignmentFile(bamfile, "rb")
	except OSError:
		print(f"{bamfile} corrupt. Skipping..")
		return

	# for ind, val in vcf.iterrows():
	# 	chr = val[vcf_chr_col]
	# 	pos = val[vcf_pos_col]
	# 	curr_reads = bam2.fetch(chr, pos-pad_bp, pos+pad_bp)

	for read in bam2: #curr_reads:
		seq = read.seq
		ref_name = read.reference_name
		quality = read.query_qualities
		align_qual_read = read.mapping_quality

		for qpos, refpos in read.get_aligned_pairs(True):
			if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
				if(quality[qpos] > base_qual):
					pileups_d[(ref_name, refpos, seq[qpos])] += 1
					curr_c = pileups_d[(ref_name, refpos, seq[qpos])]
					prior_bq = bq_d[(ref_name, refpos, seq[qpos])]
					bq_d[(ref_name, refpos, seq[qpos])] = ((prior_bq*(curr_c-1)) + quality[qpos])/curr_c
					#bq_d["counts"][refpos] += 1
					#pileups_d["nt"].append(seq[qpos])
	if len(pileups_d) != 0:
		#print('pileups_d', pileups_d)
		pileups_df = pd.DataFrame(pileups_d, index=['count']).transpose().reset_index().rename({'level_0':"chr", 'level_1':"pos", 'level_2': "nt"}, axis=1)
		bq_df = pd.DataFrame(bq_d, index=['BQ']).transpose().reset_index().rename({'level_0':"chr", 'level_1':"pos", 'level_2': "nt"}, axis=1)#.rename({0:"chr", 1:"pos", 2: "BQ"},axis=1)
		#pileups_df.to_csv(outpre + "_pileups.tsv", sep="\t")
		pileups_df = pd.merge(pileups_df, bq_df, on=["chr", "pos", "nt"])
		pileups_df["cell"] = cell
		pileups_df.to_csv(outpre + "_pileups.bq.tsv", sep="\t")
	return


def main():
	bamfile = sys.argv[1]  # Filepath to raw bamfile
	outpre = sys.argv[2]  # Prefix / basename for raw file output
	base_qual = float(
		sys.argv[3])  # Minimum base quality to be considered for pileup
	alignment_quality = float(sys.argv[4])  # minimum alignment quality required to be considered

	run_pileup(bamfile, outpre, base_qual, alignment_quality)
	return


# import click
# @click.command()
# @click.argument("bam", type=click.Path(exists=True))
# @click.argument("outpre", type=click.STRING)
# @click.option("bq", default=20)
# @click.option("align_qual", default=27)
# def main(bamfile, outpre, base_qual, alignment_quality):
# 	run_pileup(bamfile, outpre, base_qual, alignment_quality)



if __name__ == "__main__":
	main()
