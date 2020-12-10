import os
import pandas as pd
from snakemake.utils import validate
from os.path import join, basename
from glob import glob
#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"
from natsort import natsorted




experiment_names = config

fq_names_r1 = list(map(lambda x: basename(x).replace('.concat.fastq.gz',''), glob(join(config["raw_folder"] ,config['dataset'], "*R1*concat.fastq.gz" ))))
print(fq_names_r1)

fq_names_r2 = list(map(lambda x: basename(x).replace('.concat.fastq.gz',''), glob(join(config["raw_folder"] ,config['dataset'], "*R2*concat.fastq.gz" ))))
print(fq_names_r2)

rule all:
    input:
        expand("data/processed/"+config['dataset']+"/trim/trimmed_{sample}_concat_fastqc.html",experiment_name=config["dataset"],
                                                                        sample=fq_names_r1),
        expand("data/processed/"+config['dataset']+"/trim/trimmed_{sample}_concat_fastqc.html",experiment_name=config["dataset"],
                                                                        sample=fq_names_r2),
        expand("data/processed/{experiment_name}/trim/bam/trimmed_{sample}.bam", sample=, experiment_name=config["dataset"])
#sample=config["samples"]),
#expand("{raw_f}.bam.bai", raw_f=RAW_SAMPLES),

#
# def get_pair_fq(wildcards):
#     """
#     This function list all fastq files into a list
#     """
#     fst_files = natsorted([f for f in os.listdir(join(config["raw_folder"],config["dataset"])) if '_R1' in f and (f.endswith("concat.fastq.gz") or f.endswith("concat.fq.gz"))])
#     snd_files = natsorted([f for f in os.listdir(join(config["raw_folder"],config["dataset"])) if '_R2' in f and (f.endswith("concat.fastq.gz") or f.endswith("concat.fq.gz"))])
#     if len(fst_files) == len(snd_files):
#         fastqFiles = [[f1,f2] for f1,f2 in zip(fst_files,snd_files)]
#     else:
#         raise ValueError('input has single end and paired end mixed')
#     print('fastqfiles', fastqFiles)
#     return fastqFiles
    #return list(map(lambda x: basename(x), glob(join(config["raw_dir"] ,"*fastq.gz" ))))

def get_pair_fq(wildcards):
    """
    This function list all fastq files into a list
    """
    return [join(config['raw_folder'], wildcards.experiment_name, wildcards.sample+"_R1.concat.fastq.gz"),
            join(config['raw_folder'], wildcards.experiment_name, wildcards.sample+"_R2.concat.fastq.gz")]
    # fst_files = natsorted([f for f in os.listdir(join(config["raw_folder"],wildcards.experiment_name)) if wildcards.sample + '_R1' in f and (f.endswith("concat.fastq.gz") or f.endswith("concat.fq.gz"))])
    # snd_files = natsorted([f for f in os.listdir(join(config["raw_folder"],wildcards.experiment_name)) if wildcards.sample + '_R2' in f and (f.endswith("concat.fastq.gz") or f.endswith("concat.fq.gz"))])
    # print('fastqfiles', fst_files[0], snd_files[0])
    # return [fst_files[0], snd_files[0]]
    # if len(fst_files) == len(snd_files):
    #     fastqFiles = [[f1,f2] for f1,f2 in zip(fst_files,snd_files)]
    # else:
    #     raise ValueError('input has single end and paired end mixed')
    #
    # return fastqFiles
    #return list(map(lambda x: basename(x), glob(join(config["raw_dir"] ,"*fastq.gz" ))))


rule cutadapt_pair:
    input: get_pair_fq
    output:
        tr_R1 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R1_concat.fastq.gz",
        tr_R2 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R2_concat.fastq.gz"
    params:
        nt_tag = 'AAGCAGTGGTATCAACGCAGAGTAC',
        tru_seq = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGG' #"AGATCGGAAGAG" #"AGATCGGAAGAG" #ACAGTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAG
    shell:
        'cutadapt -G {params.nt_tag} -A {params.tru_seq} --minimum-length 10 --cores=16 -o {output.tr_R1} -p {output.tr_R2} {input[0]} {input[1]}'
        # shell('cutadapt -g {params.nt_tag} -a {params.tru_seq} --minimum-length 20 --cores=16 -o tmp.1.fastq -p tmp.2.fastq {input[0]} {input[1]}'),
        # shell("cutadapt -g {params.nt_tag} -a {params.tru_seq} --minimum-length 20 --cores=16 -o {output.tr_R2} -p {output.tr_R1} tmp.2.fastq tmp.1.fastq"),
        # shell("rm tmp.1.fastq tmp.2.fastq")
#TSO Sequence: AAGCAGTGGTATCAACGCAGAGTACAT

rule qc_trim:
    input:
        tr_R1 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R1_concat.fastq.gz",
        tr_R2 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R2_concat.fastq.gz"
    output:
        qc1 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R1_concat_fastqc.html",
        qc2 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R2_concat_fastqc.html"
    run:
        shell("fastqc {input[0]}"),
        shell("fastqc {input[1]}")



# rule run_star:
#     input:
#         tr_R1 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R1_concat.fastq.gz",
#         tr_R2 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R2_concat.fastq.gz"
#     output:
#         bam_f = "data/processed/{experiment_name}/trim/bam/trimmed_{sample}.bam"
#     params:
#         star_dir= config["star_dir"]
#     shell:
#         """STAR --runThreadN 4  --genomeDir {params.star_dir}
#         --readFilesIn {tr_R2} --outFileNamePrefix {output.bam_f)}
#         --outSAMtype SAM   --readNameSeparator space --outStd Log --outSAMunmapped Within KeepPairs
#         --outSAMorder PairedKeepInputOrder --outSAMattrRGline ID:{ID}:0:1:HC23FDSXY:3
#         SM:{ID} LB:0.1 PU:{ID}:0:1:HC23FDSXY:3 PL:ILLUMINA  --outSAMmultNmax 18446744073709551615
#         --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3"""


rule run_star_wiggle:
    input:
        tr_R1 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R1_concat.fastq.gz",
        tr_R2 = "data/processed/{experiment_name}/trim/trimmed_{sample}_R2_concat.fastq.gz"
    output:
        bam_f = "data/processed/{experiment_name}/trim/bam/trimmed_{sample}.bam"
    params:
        star_dir= config["star_dir"]
    shell:
        """STAR --runThreadN 4  --genomeDir {params.star_dir}  
        --readFilesIn {tr_R2} --outFileNamePrefix {output.bam_f} 
       -–outSAMtype BAM SortedByCoordinate --outWigType bedGraph --readNameSeparator space --outStd Log  --outSAMattrRGline ID:{ID}:0:1:HC23FDSXY:3  
        SM:{ID} LB:0.1 PU:{ID}:0:1:HC23FDSXY:3 PL:ILLUMINA  --outSAMmultNmax 18446744073709551615 
        --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3"""

# --outSAMunmapped Within KeepPairs
        # --outSAMorder PairedKeepInputOrder

# rule cutadapt:
#     """Trim the 5' tag"""
#     input: config['raw_folder']+"{experiment_name}/{sample}.fastq.gz"
#     output: "data/processed/{experiment_name}/trim/trimmed_{sample}.fastq.gz"
#     log:
#         'logs/trim_{experiment_name}_{sample}.log'
#     shell: "cutadapt -g AAGCAGTGGTATCAACGCAGAGTAC --cores=16 -m 1 -a AGATCGGAAGAG -o {output} {input} "
#
# #AGATCGGAAGAG
#
# rule qc_trim:
#     input: "data/processed/{experiment_name}/trim/trimmed_{sample}.fastq.gz"
#     output: "data/processed/{experiment_name}/trim/trimmed_{sample}_fastqc.html"
#     shell:
#         "fastqc {input}"

# rule multiqc:
#     aggregate

# rule index_bam:
#     """Index the bam file"""
#     input: "data/processed/{sample}/{sample}.bam"
#     output: "data/processed/{sample}/{sample}.bam.bai"
#     shell: "samtools index {input}"

