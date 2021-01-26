""" Using sample PBMC_P to extract the full CBs

"""
from src.bam_barcodes_function import extract_barcode_info
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
from os.path import join

os.chdir(ROOT_DIR)
parameters = "parameters/2020_11_18_Croker_mito.yaml"

config = read_config_file(parameters)
sample="PBMC_P"

bam_f =  join(config["results"], f"{sample}/00_bam/{sample}.bam")
out_f = join(config["results"], f"{sample}/full_barcodes.p")
print(bam_f)
print(out_f)
extract_barcode_info(bam_f, out_f,rm_slash=False, mt_chr=None)