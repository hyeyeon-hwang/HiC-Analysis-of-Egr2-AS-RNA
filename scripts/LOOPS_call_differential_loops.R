# File: LOOPS_find_differential_loops.R
# Last edited: 06 Apr 2022

library(optparse)
library(glue)
library(magrittr)
library(HiCDCPlus)

# Rscript --vanilla "${script_diffloops}" --sigindices "${sigindices}" --outdir "${outdir_diff_analysis}" --res "${res}" --distrange "50Kb-2Mb" --dmin 50000
# pval, log2fc 

args_list = list(
    make_option(c("-s","--sigindices"),
        type="character",
        default=NULL,
        help="Path to significant loops files",
        metavar="character"),
    make_option(c("--sigdir"),
        type="character",
        default=NULL,
        help="Path to output directory for significant loops analysis",
        metavar="character"),
    make_option(c("-o","--outdir"),
        type="character",
        default=NULL,
        help="Path to output directory for differential loops analysis",
        metavar="character"),
    make_option(c("-r","--res"),
        type="numeric", # cannot be "integer"
        default=NULL,
        help="Integer for resolution",
        metavar="integer"),
    make_option(c("--distrange"),
        type="character",
        default=NULL,
        help="String of distance range for loops, e.g. 50Kb-2Mb",
        metavar="character"),
    make_option(c("--dmin"),
        type="numeric",
        default=NULL,
        help="Integer for minimum distance for loops, e.g. 50000",
        metavar="integer")
)

args_parser = OptionParser(option_list=args_list)
args = parse_args(args_parser)
# print(glue("args = {args}"))

# Find differential interactions ------------------------------------------
start.diffloops <- Sys.time()
print(glue("********** START | Find differential interactions | {start.diffloops}"))

# hicdcdiff: Differential analysis using modified DESeq2 
hicdcdiff(
  input_paths = list(
    GFP = c(glue(args$sigdir,'/LentiGFP2_{args$res}.txt.gz'), glue(args$sigdir,'/LentiGFP3_{args$res}.txt.gz')),
    AS = c(glue(args$sigdir,'/LentiAS2_{args$res}.txt.gz'), glue(args$sigdir,'/LentiAS3_{args$res}.txt.gz'))),
  filter_file = args$sigindices,
  output_path = glue(args$outdir, '/'),
  #fitType = fitArg, # default local (local reg), parametric (par reg), mean (constant)
  binsize = args$res %>% as.numeric(),
  #diagnostics = TRUE,
  #DESeq.save = TRUE,
  Dmin = args$dmin %>% as.numeric(),
  #Dmax = distMax
)

end.diffloops <- Sys.time()
print(glue("********** END | Find differential interactions | {end.diffloops} | {end.diffloops - start.diffloops}"))


# Output files per chrm
# diff_resGFPoverAS_chr1.txt.gz

# diff_chr1 <- readr::read_table(glue(outdir, "/diff_analysis_chr1_chr2/diff_resGFPoverAS_chr1.txt.gz"))
# diff_chr1$log2FoldChange

