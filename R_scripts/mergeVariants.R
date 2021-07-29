merge_variants <- function(SE_list, out_f) {
  SE_list = lapply(SE_list, readRDS)
  saveRDS(cbind(SE_list), out_f)
}

args <- commandArgs(trailingOnly = TRUE)
SE_list <- args[1:length(args)-1]
out_f <- args[length(args)]
print(args)
print("SE_list")
print(SE_list)
print("out_f")
print(out_f)
merge_variants(SE_list, out_f)

sessionInfo()
