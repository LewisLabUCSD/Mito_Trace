# Introduction:
This package includes:
1) scMT: Pipeline to process single-cell MT enriched sequencing run.   
2) sim-scMT: Simulation of cell growth followed by performance metrics and assessment.

## scMT Pipeline 
This pipeline calls MT variants across cells and clusters accordingly. 
Additional output files will be created depending on the mode and the parameters.
  

### Input Files:
```
- parameters.yaml
- samples.csv
- MT_ref.fa
- paired-end fastq files (OR cellranger out folder if cellranger already run) 
```
### Output files:
``` 
- cellranger out folder (see cellranger for details)
- sc_coverage.csv
- sc_{N}.coverage.csv
- cell_af.csv # Filtered cell-by-allele frequency
- cell_by_label.csv # Cell and predicted labels depending on which method
- OPTIONAL: cell_phylogeny.json # Phylogeny reconstruction if used.
- OPTIONAL: 
    demultiplex/
        - patient_variant.csv
        - coverage folder for each patient separately.
        - Overlap variants vs patient-specific overlap


# Figures:
- MT_cells_reads_thresh.png
- MT_position_coverage.png
- overlap_variants_with_without_ligand.png
- lineage_OnlyoverlapVariants.png
- lineage_allOveralps.png
- Cellranger QC across samples


```
----

There are a few different data types the pipeline can work with. 
They are both run using the 10x platform. These are:  
## A. mtscATAC-seq
This is a new sequencing data type from Lareau, Ludwig et al. 20. 
#### Additional Input:

#### Additional Output:
```
- peaks.tsv.gz # ATAC-seq peak files
- sc.coverage.tsv # Additional Columns:  Coverage +, Coverage -
- sc.{N}.coverage.tsv # Additional Columns:  Coverage +, Coverage -
- sc.cell_af_meta.tsv: Strand concordance values # Correlation variants on each strand

```

### Run:
 In order to run, use:
 `snakemake --snakefile snakefile --configfile {parameter_f}`
 
 `Parameter file:`



## B. 5'RNA-seq. 


Additionally, there are different methods to be run on them.  
i.  Single run per patient  
ii. Demultiplexed (NOT IMPLEMENTED YET)


## sim-scMT 




## Installation:
With conda:  
`conda install src`


## Modes: 


## Input: 


## Output:


## Pacakge Dependencies: 


## Future Directions: 