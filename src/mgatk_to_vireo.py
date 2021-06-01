import click
from src.utils.data_io import mgatk_to_vireo
from os.path import join
@click.command()
@click.argument("in_dir", type=click.Path(exists=True))
@click.argument("out_dir", type=click.Path(exists=True))
@click.argument("sample", type=click.STRING)
def main(in_dir, out_dir, sample):
    #samples = ','.split(samples)
    in_af = join(in_dir, f"{sample}.af.tsv")
    in_af_meta = join(in_dir, f"{sample}.af.mgatk.tsv")
    in_coverage = join(in_dir, f"{sample}.coverage.tsv")
    mgatk_to_vireo(in_af, in_af_meta, in_coverage, outdir=out_dir, out_name="cellSNP")
    return


if __name__ == '__main__':
    main()