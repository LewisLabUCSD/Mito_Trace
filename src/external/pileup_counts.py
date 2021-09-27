#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import pysam
pysam.set_verbosity(0)
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


def run_pileup(bamfile, outpre, maxBP, base_qual, sample, alignment_quality, use_strands=False):
	n = int(maxBP)

	is_corrupt = False
	# BAQ
	# initialize with a pseudo count to avoid dividing by zero
	countsA = [0.00000001] * n
	countsC = [0.00000001] * n
	countsG = [0.00000001] * n
	countsT = [0.00000001] * n
	qualA = [0.0] * n
	qualC = [0.0] * n
	qualG = [0.0] * n
	qualT = [0.0] * n

	if use_strands:
		# initialize with a pseudo count to avoid dividing by zero
		neg_countsA = [0.00000001] * n
		neg_countsC = [0.00000001] * n
		neg_countsG = [0.00000001] * n
		neg_countsT = [0.00000001] * n
		neg_qualA = [0.0] * n
		neg_qualC = [0.0] * n
		neg_qualG = [0.0] * n
		neg_qualT = [0.0] * n
	try:
		bam2 = pysam.AlignmentFile(bamfile, "rb")
	except OSError:
		print(f"{bamfile} corrupt. Skipping..")
		return
	for read in bam2:
		seq = read.seq
		quality = read.query_qualities
		align_qual_read = read.mapping_quality
		reverse = read.is_reverse
		if use_strands:
			if reverse:
				for qpos, refpos in read.get_aligned_pairs(True):
					if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
						if (seq[qpos] == "A" and quality[
							qpos] > base_qual):
							neg_qualA[refpos] += quality[qpos]
							neg_countsA[refpos] += 1
						elif (seq[qpos] == "C" and quality[
							qpos] > base_qual):
							neg_qualC[refpos] += quality[qpos]
							neg_countsC[refpos] += 1
						elif (seq[qpos] == "G" and quality[
							qpos] > base_qual):
							neg_qualG[refpos] += quality[qpos]
							neg_countsG[refpos] += 1
						elif (seq[qpos] == "T" and quality[
							qpos] > base_qual):
							neg_qualT[refpos] += quality[qpos]
							neg_countsT[refpos] += 1
			else:
				for qpos, refpos in read.get_aligned_pairs(True):
					# print('qpos', qpos)
					# print('refpos', refpos)
					if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
						if (seq[qpos] == "A" and quality[
							qpos] > base_qual):
							qualA[refpos] += quality[qpos]
							countsA[refpos] += 1
						elif (seq[qpos] == "C" and quality[
							qpos] > base_qual):
							qualC[refpos] += quality[qpos]
							countsC[refpos] += 1
						elif (seq[qpos] == "G" and quality[
							qpos] > base_qual):
							qualG[refpos] += quality[qpos]
							countsG[refpos] += 1
						elif (seq[qpos] == "T" and quality[
							qpos] > base_qual):
							qualT[refpos] += quality[qpos]
							countsT[refpos] += 1

		else:
			for qpos, refpos in read.get_aligned_pairs(True):
				if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
					if(seq[qpos] == "A" and quality[qpos] > base_qual):
						qualA[refpos] += quality[qpos]
						countsA[refpos] += 1
					elif(seq[qpos] == "C" and quality[qpos] > base_qual):
						qualC[refpos] += quality[qpos]
						countsC[refpos] += 1
					elif(seq[qpos] == "G" and quality[qpos] > base_qual):
						qualG[refpos] += quality[qpos]
						countsG[refpos] += 1
					elif(seq[qpos] == "T" and quality[qpos] > base_qual):
						qualT[refpos] += quality[qpos]
						countsT[refpos] += 1


	meanQualA = [round(x/y,1) for x, y in zip(qualA, countsA)]
	meanQualC = [round(x/y,1) for x, y in zip(qualC, countsC)]
	meanQualG = [round(x/y,1) for x, y in zip(qualG, countsG)]
	meanQualT = [round(x/y,1) for x, y in zip(qualT, countsT)]

	countsA = [ int(round(elem)) for elem in countsA ]
	countsC = [ int(round(elem)) for elem in countsC ]
	countsG = [ int(round(elem)) for elem in countsG ]
	countsT = [ int(round(elem)) for elem in countsT ]
	# Allele Counts
	minBP = 0

	writeSparseMatrix2("A", countsA, meanQualA, outpre, maxBP, sample)
	writeSparseMatrix2("C", countsC, meanQualC, outpre, maxBP, sample)
	writeSparseMatrix2("G", countsG, meanQualG, outpre, maxBP, sample)
	writeSparseMatrix2("T", countsT, meanQualT, outpre, maxBP, sample)

	zipped_list = zip(list(countsA),list(countsC),list(countsG),list(countsT))
	sums = [sum(item) for item in zipped_list]
	writeSparseMatrix("coverage", sums, outpre, maxBP, sample)

	if use_strands: #Write the negative strand
		neg_meanQualA = [round(x / y, 1) for x, y in zip(neg_qualA, neg_countsA)]
		neg_meanQualC = [round(x / y, 1) for x, y in zip(neg_qualC, neg_countsC)]
		neg_meanQualG = [round(x / y, 1) for x, y in zip(neg_qualG, neg_countsG)]
		neg_meanQualT = [round(x / y, 1) for x, y in zip(neg_qualT, neg_countsT)]

		neg_countsA = [int(round(elem)) for elem in neg_countsA]
		neg_countsC = [int(round(elem)) for elem in neg_countsC]
		neg_countsG = [int(round(elem)) for elem in neg_countsG]
		neg_countsT = [int(round(elem)) for elem in neg_countsT]
		# Allele Counts
		writeSparseMatrix2("A.minus", neg_countsA, neg_meanQualA, outpre, maxBP,
						   sample)
		writeSparseMatrix2("C.minus", neg_countsC, neg_meanQualC, outpre, maxBP,
						   sample)
		writeSparseMatrix2("G.minus", neg_countsG, neg_meanQualG, outpre, maxBP,
						   sample)
		writeSparseMatrix2("T.minus", neg_countsT, neg_meanQualT, outpre, maxBP,
						   sample)

		zipped_list = zip(list(neg_countsA), list(neg_countsC),
						  list(neg_countsG), list(neg_countsT))
		sums = [sum(item) for item in zipped_list]
		writeSparseMatrix("coverage.minus", sums, outpre, maxBP, sample)


def main():
	bamfile = sys.argv[1]  # Filepath to raw bamfile
	outpre = sys.argv[2]  # Prefix / basename for raw file output
	maxBP = sys.argv[3]  # Maximum length of mtDNA genome
	base_qual = float(
		sys.argv[4])  # Minimum base quality to be considered for pileup
	sample = sys.argv[
		5]  # Sample name to be considered for downstream analyses
	alignment_quality = float(sys.argv[6])  # minimum alignment quality required to be considered

	run_pileup(bamfile, outpre, maxBP, base_qual, sample, alignment_quality)
	return


if __name__ == "__main__":
	main()
