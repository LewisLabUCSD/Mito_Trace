# Mito_Trace  
Tracking cell-lineage using mutations in the mitochondrion genome.  
  

## Important Steps 
1. Run cellranger 10x on fastq files 
2. Preprocess: Convert bam file into single-cell pileup data
3. Filtering: Coverage and quality based filtering.
4. Variant calling: Filtering MT variants using Vireo or mgatk
5. Merging conditions
6. Demultiplexing: Assigning a donor to each cell, and None if it is unassigned.
7. Clonal detection: Assign cells to a clone and None if unassigned.
8. Clonal enrichment across conditions: Run enrichment analysis across conditions of the same sample/donor

### Ways the steps can be run:
There are different ways the pipeline process can be run, in different orders, 
such that the proper quality and parameters can be assessed.
Steps: 
[1->2] is done first, and [7->8] is done last. 
In the middle, it's either [3->5->4->6], [3->4->5->6] or [5->6->3->4]
