# File: LOOPS_generate_features.R
# Last edited: 06 Apr 2022 

library(optparse)
library(glue)
library(magrittr)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(HiCDCPlus)

args_list = list(
    make_option(c("-o", "--outdir"),
                type="character",
                default=NULL,
                help="Path to output directory for features file",
                metavar="character"),
    make_option(c("-r", "--res"),
                type="numeric", 
                default=NULL,
                help="Integer for resolution",
                metavar="integer"),
    make_option(c("-c", "--chrs"),
                type="character",
                default="all",
                help="One chromosome (e.g. 'chr5'). Default is 'all' for all chromosomes",
                metavar="character")
)

args_parser = OptionParser(option_list=args_list)
args = parse_args(args_parser)

# Generate features -------------------------------------------------------
features_prefix_path <- glue(args$outdir, "/rn6_{args$res}_GATC_GANTC_features")

if (args$chrs == "all") {
     construct_features(
          output_path = features_prefix_path,
          gen = "Rnorvegicus",
          gen_ver = "rn6",
          sig = c("GATC", "GANTC"),
          bin_type = "Bins-uniform",
          binsize = args$res,
          chrs = sapply(c(1:20, "X", "Y"), function(x) {glue("chr", x)}) %>% as.character())  
} else {
     construct_features(
          output_path = features_prefix_path,
          gen = "Rnorvegicus",
          gen_ver = "rn6",
          sig = c("GATC", "GANTC"),
          bin_type = "Bins-uniform",
          binsize = args$res,
          chrs = args$chrs)  
}














