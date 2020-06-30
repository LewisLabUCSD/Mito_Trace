# This script takes all the A,C,G,T pileups for single cells and merges them
directory=$1
samplename=$2
#Samplename is full directory path plus the name: e.g. 'all_coverage/A'

echo ${directory}
echo ${samplename}
# Combine all files
cat ${directory}/*.A.txt | gzip > "${samplename}_all.A.txt.gz"
cat ${directory}/*.C.txt | gzip > "${samplename}_all.C.txt.gz"
cat ${directory}/*.G.txt | gzip > "${samplename}_all.G.txt.gz"
cat ${directory}/*.T.txt | gzip > "${samplename}_all.T.txt.gz"
cat ${directory}/*.coverage.txt | gzip > "${samplename}_all.coverage.txt.gz"

