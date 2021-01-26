# Per rule explanation

This is an automatically generated list of all supported rules, their docstrings, and command. At the start of each workflow run a list is printed of which rules will be run. And while the workflow is running it prints which rules are being started and finished. This page is here to give an explanation to the user about what each rule does, and for developers to find what is, and isn't yet supported.

#### all
 Call variants using one of a few types. 

:type type = {'thresh', 'GMM', 'mgatk'}
Modes: 
thresh: Sets a simple threshold for number of cells and reads a
variant needs to be called a variant. 
GMM: Gaussian mixture model for quality of the variants. 
mgatk: This will use strand concordance, allele quality, and 
number of reads to determine the correct threshold. This only
works information on aligned read strands in scATAC. 
 
Include a model for known variants, such as transversion

#### model_call_lineages
 Looks at performance for cell-clone relationship
and MTvariant-clone relationship

If a tree and not cluster labels are given, will assume MT lineages
are based on variants, 

#### model_lineage_growth_estimate
 Will estimate the growth rate based on the number of cells in
the before and after experiments.

#### performance_growth_estimate
 Compare the growth of the largest clones between the disease and
non-disease samples. 

