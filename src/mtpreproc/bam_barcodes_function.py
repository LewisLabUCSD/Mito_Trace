from tqdm import tqdm
import pysam
from collections import defaultdict
import pickle
import click


def extract_barcode_info(bam_f, out_f,rm_slash=False, mt_chr="MT"):
    """
    Args:
        bam_f:
        out_f:
        rm_slash:
    """
    print('mt_chr', mt_chr)
    samfile = pysam.AlignmentFile(bam_f, "rb")
    barcodes = set()
    corrected_barcodes = set()
    barcode_pairs = set()
    CR_read_number = defaultdict(int)
    CB_read_number = defaultdict(int)
    BC_read_number = defaultdict(int)

    for read in tqdm(samfile.fetch(mt_chr)):
        d = dict(read.tags)
        if "CR" in d:
            barcodes.add(d["CR"])
            CR_read_number[d["CR"]] += 1
            if "CB" in d:
                if rm_slash:
                    CB = d["CB"].split('-')[0]
                else:
                    CB = d["CB"]
                barcode_pairs.add((CB, d["CR"]))
                corrected_barcodes.add(CB)
                CB_read_number[CB] += 1
            else:
                barcode_pairs.add(('-', d["CR"]))
        else:
            if "CB" in d:
                if rm_slash:
                    CB = d["CB"].split['-'][0]
                else:
                    CB = d["CB"]
                #pair_barcodes.add((CB, '-'))
                corrected_barcodes.add(CB)
                CB_read_number[CB] += 1
        if "BC" in d:
            BC_read_number[d["BC"]] += 1 
    samfile.close()
    pickle.dump([CR_read_number,CB_read_number,BC_read_number, barcodes, corrected_barcodes, barcode_pairs], file=open(out_f, "wb"))

    return CR_read_number,CB_read_number,BC_read_number, barcodes, corrected_barcodes, barcode_pairs



@click.command(help="Extract Cell Barcodes from bam")
@click.argument('bam_f',  type=click.Path(exists=True))
@click.argument('out_f',  type=click.Path(exists=False))
@click.argument('mtchr',  type=click.STRING)
def main(bam_f, out_f, mtchr):
    # bam_f = sys.argv[1]
    # out_f = sys.argv[2]
    """
    Args:
        bam_f:
        out_f:
    """
    extract_barcode_info(bam_f, out_f, rm_slash=False, mt_chr=mtchr)

    return


if __name__== '__main__':
    main()