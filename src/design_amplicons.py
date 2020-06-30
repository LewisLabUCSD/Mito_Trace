from Bio import SeqIO
import click
import pandas as pd
from src.config import ROOT_DIR
from os.path import join

def amplicons_over_entire_region(ref_fa, out_f, chromosome="MT", amplicon_length=150, pcr_length=(18,22)):
    record_dict = SeqIO.to_dict(SeqIO.parse(ref_fa, "fasta"))
    try:
        seq = record_dict[chromosome]
    except KeyError:
        print(f"{chromosome} not found in {ref_fa}. Please check if you have the proper name")
        return

    bed_df = pd.DataFrame(columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand"])
    #sequences = SeqIO.fas
    for ind, subseq in range(0,len(seq), pcr_length[1]+ amplicon_length):
        bed_df = pd.concat(bed_df, pd.DataFrame({"chrom":chromosome,
                                        "chromStart":subseq, "chromEnd":subseq+pcr_length[1], "name":ind, "score":0, "strand":"+"}))
        print(ind, subseq)

    bed_df.to_csv(out_f + ".bed", header=None, sep="\t", index=False)
    #SeqIO.write(sequences,open(out_f + ".fasta", "w"))
    return


@click.command(help="Create a bed file and fasta file of amplicon regions")
@click.option('--length', default=150, help='The length of the amplicon')
@click.option('--chromosome', default="MT", help='Which chromosome to target')
@click.option('--pcr_length', default=(18,22),type=click.Tuple(int,int), help='The PCR min and max')
@click.argument('ref_fa', type=click.Path(exists=True)) #help='file to append the command',
@click.argument('out_f', type=click.Path())#coverage_folder, f_save, maxBP = 16571, f_save_fig=None
def main_commandline(ref_fa, out_f, chromosome="MT", amplicon_length=150, pcr_length=(18,22)):
    amplicons_over_entire_region(ref_fa, out_f, chromosome, amplicon_length, pcr_length)
    return


def main(ref_fa, out_f, chromosome="MT", amplicon_length=150, pcr_length=(18,22)):
    amplicons_over_entire_region(ref_fa, out_f, chromosome, amplicon_length, pcr_length)
    return


if __name__ == "__main__":
    ref_fa = join(ROOT_DIR, "data/external/genomes/genome.fa")
    out_f = join(ROOT_DIR, "data/processed/amplicons")
    main(ref_fa, out_f)
