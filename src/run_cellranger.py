import sys
import click
import os


def project_bam2fastq():
    in_f = "/data2/mito_lineage/data/raw/GSM3271867/Zhu_Single_Cell_10X_genomics_Human2_002_possorted_genome_bam.bam"
    out_dir = "/data2/mito_lineage/data/processed/GSM3271867/fastq"
    bam2fastq(in_f, out_dir)
    return



#@click.group(invoke_without_command=True)
@click.group()
@click.pass_context
def cellranger(ctx):
    # if ctx.invoked_subcommand is None:
    #     click.echo('I was invoked without subcommand')
    # else:
    #     click.echo('I am about to invoke %s' % ctx.invoked_subcommand)
    click.echo('I am about to invoke %s' % ctx.invoked_subcommand)

    return


@cellranger.command()
@click.argument('in_f', type=click.Path(exists=True))
@click.argument('out_dir', type=click.Path())
def bam2fastq(in_f, out_dir):
    """Run bamtofastq to convert bam to fastq file. Must be an output from 10x"""
    cmd = f"cellranger bamtofastq {in_f}  {out_dir}"
    print(cmd)
    os.system(cmd)
    return


@cellranger.command()
@click.argument('fastq_path', type=click.Path(exists=True))
@click.argument('id', type=click.STRING)
@click.argument('sample_name', type=click.STRING)
@click.option('--txn', type=click.Path(exists=True))
@click.option('--cells', type=click.INT, default=0)
@click.option('--cores', type=click.INT, default=16)
def count(fastq_path, id, sample_name, txn="/data2/genome/human_GRCh38/cellranger/refdata-cellranger-GRCh38-3.0.0", cells=0, cores=16):
    """
    Run cellranger count pipeline
    :param fastq_path: Path to fastq files
    :param id: A unique run id, used to name output folder
    :param sample_name: Prefix of the filenames of FASTQs to select.
    :param txn: Path of folder containing 10x-compatible reference
    :param cells: Expected number of cells.
    :return:
    """
    args = ""
    if cells != 0:
        args += f"--expect-cells={cells} "
    cmd = f"cellranger count --id={id} \
                   --transcriptome={txn} \
                   --fastqs={fastq_path} \
                   --sample={sample_name} \
                --localcores={cores} {args}"
    print(cmd)
    #os.system(cmd)
    return



@cellranger.command()
@click.option('--cores', type=click.INT, default=16)
@click.option('--path', type=click.STRING, default="isshamie")
def project_count(cores, path):
    """
    Run cellranger count with files defined in this function.
    :return:
    """
    print(path)
    if path == "isshamie":
        os.chdir("/data2/isshamie/mito_lineage/data/processed/GSM3271867")
    else:
        os.chdir("/data2/mito_lineage/data/processed/GSM3271867")
    id= "Human2_002"
    txn = "/data2/genome/human_GRCh38/cellranger/refdata-cellranger-GRCh38-3.0.0"
    sample_name = "bamtofastq"
    fastq_path = "/data2/mito_lineage/data/processed/GSM3271867/fastq/016_16_Human_P_h2_G8_MissingLibrary_1_CAL34ANXX"

    cmd = f"""umask g+x,g+w && cellranger count --id={id} \
                   --transcriptome={txn} \
                   --fastqs={fastq_path} \
                   --sample={sample_name} \
                --localcores={cores}"""
    print(cmd)
    os.system(cmd)
    return


if __name__ == "__main__":
    cellranger()
    # bam2fastq()
    #
