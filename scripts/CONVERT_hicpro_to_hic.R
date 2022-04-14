options(scipen = 999)
#source("./scripts/generate_features.R")
library(HiCDCPlus)
library(magrittr)
library(glue)

indir_features <- "./output/CONVERT/hicdcplus_features"
indir_hicpro <-  "./data/hicProMerged" # OSCAR: "~/scratch/data/hicProMerged"
outdir_hic <- "./output/CONVERT/converted_hic/hic_all_intra"

valid_chrs <- c(1:20, "X") %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
conditions <- c("as", "gfp") # conditions <- c("as")

res_list <- c(10000, 100000)

# Generate rn6 chromosome size file 
get_chr_sizes(gen = "Rnorvegicus", gen_ver = "rn6", chrs = valid_chrs) %>%
    as.data.frame() %>%
    tibble::add_column(chrom = rownames(.), .before = 1) %>%
    readr::write_tsv(., file = "./output/CONVERT/rn6_chrom.sizes", col_names = FALSE)
rn6_chrom_sizes <- "./output/CONVERT/rn6_chrom.sizes"

# max_chrom_size <- get_chr_sizes(gen = "Rnorvegicus", gen_ver = "rn6", chrs = valid_chrs) %>% max()

convert_hicpro_to_hic <- function(all_intra){
    for(cond in conditions){
        for(res in res_list){
            path_bed <- glue("{indir_hicpro}/Lenti{toupper(cond)}/raw/{res}/Lenti{toupper(cond)}_{res}_abs.bed")
            path_mat <- glue("{indir_hicpro}/Lenti{toupper(cond)}/iced/{res}/Lenti{toupper(cond)}_{res}_iced.matrix")
            gi_list_bin <- NULL
            gi_list_hicpro <- NULL

            if (all_intra == FALSE){
                print("full data FALSE")
                gi_list_bin <- generate_bintolen_gi_list(
                    bintolen_path = glue("{indir_features}/rn6_{res}_GATC_GANTC_features_bintolen.txt.gz"),
                    chr = valid_chrs,
                    binsize = glue("{res}") %>% as.numeric(),
                    gen = "Rnorvegicus",
                    gen_ver = "rn6") 
                gi_list_hicpro <- add_hicpro_matrix_counts(
                    gi_list = gi_list_bin,
                    absfile_path = path_bed,
                    matrixfile_path = path_mat)
                hicdc2hic(gi_list = gi_list_hicpro,
                          hicfile = glue("{outdir_hic}/{cond}_{res}_merged_iced.hic"),
                          mode = "raw",
                          chrs = valid_chrs,
                          gen_ver = rn6_chrom_sizes,
                          memory = 25)
            } else {
                print("all intra data TRUE")
                print(glue("{cond} {res}"))
                gi_list_bin <- generate_bintolen_gi_list(
                    bintolen_path = glue("{indir_features}/rn6_{res}_GATC_GANTC_features_bintolen.txt.gz"),
                    chr = valid_chrs,
                    binsize = glue("{res}") %>% as.numeric(),
                    gen = "Rnorvegicus",
                    gen_ver = "rn6",
                    Dthreshold = get_chr_sizes(gen = "Rnorvegicus", gen_ver = "rn6", chrs = valid_chrs) %>% sum()
                )
                print("made gi_list_bin")
                gi_list_hicpro <- add_hicpro_matrix_counts(
                    gi_list = gi_list_bin,
                    absfile_path = path_bed,
                    matrixfile_path = path_mat,
                    add_inter = TRUE)
                print("made gi_list_hicpro")
                hicdc2hic(gi_list = gi_list_hicpro,
                          hicfile = glue("{outdir_hic}/{cond}_{res}_merged_iced_all_intra.hic"),
                          mode = "raw",
                          chrs = valid_chrs,
                          gen_ver = rn6_chrom_sizes,
                          memory = 120) # Memory in GB to allocate, will get segfault if not enough memory
                print("made hic file")
            }
        }
    }
}       
convert_hicpro_to_hic(all_intra = TRUE)
