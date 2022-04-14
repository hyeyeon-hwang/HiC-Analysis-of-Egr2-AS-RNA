library(glue)
library(magrittr)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(HiCDCPlus)

# Generate features -------------------------------------------------------

outdir <- "./output/CONVERT/hicdcplus_features"
if (!dir.exists(outdir)){dir.create(outdir)} else {print(glue("{outdir} exists!"))}

res_list <- c(10000, 100000)

path_features <- vector(mode = "character", length = length(res_list))

for(i in 1:length(res_list)){
    res <- res_list[i]
    path_features[i] <- construct_features(
        output_path = glue("{outdir}/rn6_{res}_GATC_GANTC_features"),
        gen = "Rnorvegicus",
        gen_ver = "rn6",
        sig = c("GATC", "GANTC"),
        bin_type = "Bins-uniform",
        binsize = res,
        chrs = sapply(c(1:20, "X"), function(x) {glue("chr", x)}) %>% as.character()
    )
}
