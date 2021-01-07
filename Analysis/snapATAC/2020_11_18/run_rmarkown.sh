#!/bin/bash
Rscript -e "rmarkdown::render('snapATAC_01.Rmd',params=list(configfile='snapATAC_2020_11_18_PBMC_J.yaml'), output_file='data/out/PBMC_J_01.pdf')"

Rscript -e "rmarkdown::render('snapATAC_02.Rmd',params=list(configfile='snapATAC_2020_11_18_PBMC_J.yaml'),output_file='data/out/PBMC_J_02.pdf')"
