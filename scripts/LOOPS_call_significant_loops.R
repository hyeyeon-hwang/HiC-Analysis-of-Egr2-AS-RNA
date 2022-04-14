# File: LOOPS_find_significant_loops.R
# Last edited: 06 Apr 2022

library(optparse)
library(glue)
library(magrittr)
library(HiCDCPlus)

args_list = list(
    make_option(c("-d", "--datadir"),
        type="character",
        default=NULL,
        help="Path to data directory of *.matrix and *.bed HiC-Pro files for all samples",
        metavar="character"),
    make_option(c("-f","--features"),
        type="character",
        default=NULL,
        help="Path to features file previously generated",
        metavar="character"),
    make_option(c("-o","--outdir"),
        type="character",
        default=NULL,
        help="Path to output directory for significant loops analysis",
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
        metavar="integer"),
    make_option(c("-s", "--sampsize"),
        type="numeric",
        default=NULL,
        help="Sample size for downsampling, e.g. 0.01",
        metavar="numeric"),
    make_option(c("-q","--qval"),
        type="numeric",
        default=NULL,
        help="q-value for significant loops threshold",
        metavar="numeric")
)


args_parser = OptionParser(option_list=args_list)
args = parse_args(args_parser)
# print(glue("args = {args}"))

# Find significant interactions -------------------------------------------
start.sigloops <- Sys.time()
print(glue("********** START | Find significant interactions | {start.sigloops}"))

# Data files path ---------------------------------------------------------
matfile_paths <- c(
  glue(args$datadir, "/LentiAS2_", args$res, ".matrix"),
  glue(args$datadir, "/LentiAS3_", args$res, ".matrix"),
  glue(args$datadir, "/LentiGFP2_", args$res, ".matrix"),
  glue(args$datadir, "/LentiGFP3_", args$res, ".matrix")
)
bedfile_paths <- c(
  glue(args$datadir, "/LentiAS2_", args$res, "_abs.bed"),
  glue(args$datadir, "/LentiAS3_", args$res, "_abs.bed"),
  glue(args$datadir, "/LentiGFP2_", args$res, "_abs.bed"),
  glue(args$datadir, "/LentiGFP3_", args$res, "_abs.bed")
)

indexfile_all <- data.frame()
indexfile <- data.frame()
for (k in 1:length(matfile_paths)){
  matfile_path <- matfile_paths[k]
  bedfile_path <- bedfile_paths[k]
  gilist_outpath <- paste0(
    args$outdir,'/',
    gsub("^(.*[\\/])", "", gsub('.matrix','.txt.gz', matfile_path))
  )
  # Generate gi_list instance
  # bintolen_path = output of generate features

  # Exclude chrY b/c of error
  # Error in h(simpleError(msg, call)) : error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': unable to find an inherited method for function 'regions' for signature '"NULL"'
  gi_list <- generate_bintolen_gi_list(
    bintolen_path = args$features,
    gen = "Rnorvegicus",
    gen_ver = "rn6",
    chrs = sapply(c(1:20, 'X'), function(x) {glue("chr", x)}) %>% as.character()
  )

  # Add HiC-Pro counts to gi_list
  gi_list <- add_hicpro_matrix_counts(
    gi_list = gi_list,
    absfile_path = bedfile_path,
    matrixfile_path = matfile_path
  )

  # Expand features for modeling
  gi_list <- expand_1D_features(gi_list)
  # HiC-DC downsamples rows for modeling
  # ssize = stratified sampling size, default = 0.01
  # can decrease for large chroms, increase recommended if model fails to converge

  set.seed(1010) # HiC-DC downsamples rows for modeling
  gi_list <- HiCDCPlus(gi_list, Dmin=as.numeric(args$dmin))  
  # gi_list <- HiCDCPlus(gi_list, Dmin=50000)  

  # Write results to a text file: Lenti{exp}_10000.txt.gz
  # saveRDS(gi_list, file = glue(gilist_outpath, ".rds"))

  print(glue("seq(length(gi_list)) = {seq(length(gi_list))}"))
  # i = 1:num of chroms
  indexfile_sample_all <- data.frame()
  #indexfile_sample_unique <- data.frame()
  for (i in seq(length(gi_list))){
    gi_list_noNA <- gi_list[[i]][!is.na(gi_list[[i]]$qvalue)]

    indexfile <- unique(rbind(
      indexfile,
      #indexfile_sample_all[c('seqnames1','start1','start2')]
      as.data.frame(gi_list_noNA[gi_list_noNA$qvalue<=as.numeric(args$qval)])[c('seqnames1','start1','start2')]
      #as.data.frame(gi_list_noNA[gi_list_noNA$qvalue<=0.05])[c('seqnames1','start1','start2')]
      # 'seqnames1', 'start1', 'end1', 'width1', 'strand1', 'gc1', 'len1',
      # 'seqnames2', 'start2', 'end2', 'width2', 'strand2', 'gc2', 'len2',
      # 'D', 'counts', 'gc', 'len')]
    ))
    indexfile_sample_all <- rbind(indexfile_sample_all,
      as.data.frame(gi_list_noNA[gi_list_noNA$qvalue<=as.numeric(args$qval)])[c('seqnames1','start1','start2')]
    )
    indexfile_all <- rbind(indexfile_all,
      as.data.frame(gi_list_noNA[gi_list_noNA$qvalue<=as.numeric(args$qval)])[c('seqnames1','start1','start2')]
    )
    #data.table::fwrite(
    #  indexfile_sample_all,
    #  glue(args$outdir, '/{args$res}_{args$distrange}_analysis_sigindices_sample{i}_all.txt.gz'),
    #  sep='\t', row.names=FALSE, quote=FALSE
    #)
    #data.table::fwrite(
    #  indexfile_sample_unique,
    #  glue(args$outdir, '/{args$res}_{args$distrange}_analysis_sigindices_sample{i}_unique.txt.gz'),
    #  sep='\t', row.names=FALSE, quote=FALSE
    #)
    
  } # rbind, unique, for
  print(glue('Writing gi_list to {gilist_outpath}'))
  gi_list_write(gi_list, fname = gilist_outpath)
  print(glue('Wrote gi_list to {gilist_outpath}'))

  
  sigindices_sample_outpath <- glue(args$outdir,'/supplementary/{args$res}_{args$distrange}_analysis_sigindices_{k}.txt.gz')
  colnames(indexfile_sample_all) <- c('chr','startI','startJ')
  data.table::fwrite(
    indexfile_sample_all,
    sigindices_sample_outpath,
    sep='\t', row.names=FALSE, quote=FALSE
  )
  
} # end: for

# Save index file---union of significant interactions at 10kb
# dim(indexfile) 
colnames(indexfile) <- c('chr','startI','startJ')
colnames(indexfile_all) <- c('chr','startI','startJ')
# Output index file of significant interaction indices
sigindices_outpath <- glue(args$outdir,'/{args$res}_{args$distrange}_analysis_sigindices.txt.gz')
data.table::fwrite(
  indexfile,
  sigindices_outpath,
  sep='\t', row.names=FALSE, quote=FALSE
)
sigindices_all_outpath <- glue(args$outdir,'/supplementary/{args$res}_{args$distrange}_analysis_sigindices_notunique.txt.gz')
data.table::fwrite(
  indexfile_all,
  sigindices_all_outpath,
  sep='\t', row.names=FALSE, quote=FALSE
)

end.sigloops <- Sys.time()
print(glue("********** END | Find significant interactions | {end.sigloops} | {end.sigloops - start.sigloops}"))















