import sys
import os

## Taken from: https://www.biostars.org/p/13462/
def main(input_b, output_b):
    cmd = f"""samtools view -H {input_b} | \
       sed -e 's/SN:chr1/SN:1/' | sed -e 's/SN:chr2/SN:2/' | \
       sed -e 's/SN:chr3/SN:3/' | sed -e 's/SN:chr4/SN:4/' | \
       sed -e 's/SN:chr5/SN:5/' | sed -e 's/SN:chr6/SN:6/' | \
       sed -e 's/SN:chr7/SN:7/' | sed -e 's/SN:chr8/SN:8/' | \
       sed -e 's/SN:chr9/SN:9/' | sed -e 's/SN:chr10/SN:10/' | \
       sed -e 's/SN:chr11/SN:11/' | sed -e 's/SN:chr12/SN:12/' | \
       sed -e 's/SN:chr13/SN:13/' | sed -e 's/SN:chr14/SN:14/' | \
       sed -e 's/SN:chr15/SN:15/' | sed -e 's/SN:chr16/SN:16/' | \
       sed -e 's/SN:chr17/SN:17/' | sed -e 's/SN:chr18/SN:18/' | \
       sed -e 's/SN:chr19/SN:19/' | sed -e 's/SN:chr20/SN:20/' | \
       sed -e 's/SN:chr21/SN:21/' | sed -e 's/SN:chr22/SN:22/' | \
       sed -e 's/SN:chrX/SN:X/' | sed -e 's/SN:chrY/SN:Y/' | \
       sed -e 's/SN:chrM/SN:MT/' | samtools reheader - {input_b} > {output_b}""".replace('       ','')
    print(cmd)
    os.system(cmd)
    return

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

