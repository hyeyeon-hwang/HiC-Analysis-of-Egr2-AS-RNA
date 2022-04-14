# Load libraries ------------------------------------------------------------------------------

library(glue)
library(magrittr)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(GenomicRanges)
library(janitor)
library(tibble)
library(xlsx)


rds_dir <- "./output/R_objects"

# DATA: ATAC-seq ------------------------------------------------------------------------------
data_dir <- "./data/ATACseq"

np_cols <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
             "signalValue", "pValue", "qValue", "peak")
as_peaks <- glue::glue("{data_dir}/as_NA_peaks.narrowPeak") %>%
    readr::read_tsv(., col_names = np_cols) # 60126 x 10

gfp_peaks <- glue::glue("{data_dir}/gfp_NA_peaks.narrowPeak") %>%
    readr::read_tsv(., col_names = np_cols) # 64896 x 10

saveRDS(as_peaks, glue::glue("{rds_dir}/as_peaks.rds"))
saveRDS(gfp_peaks, glue::glue("{rds_dir}/gfp_peaks.rds"))
# END of DATA ATAC-seq

# DATA: COREs ---------------------------------------------------------------------------------
cores_dir <- "./data/COREs"
cores_as1 <- glue(cores_dir, "/lentiAS1_COREs.bed") %>%
    read.table(., sep = "\t", col.names = c("chr", "start", "end")) %>%
    tibble::add_column(myid_core = sapply(c(1:nrow(.)), function(i){glue('c.as1.{i}')})) %>%
    tibble::add_column(coretype = sapply(c(1:nrow(.)), function(i){glue('as1')})) %>%
    tibble::add_column(coreset = sapply(c(1:nrow(.)), function(i){glue('as')}))
cores_as2 <- glue(cores_dir, "/lentiAS2_COREs.bed") %>%
    read.table(., sep = "\t", col.names = c("chr", "start", "end")) %>%
    tibble::add_column(myid_core = sapply(c(1:nrow(.)), function(i){glue('c.as2.{i}')})) %>%
    tibble::add_column(coretype = sapply(c(1:nrow(.)), function(i){glue('as2')})) %>%
    tibble::add_column(coreset = sapply(c(1:nrow(.)), function(i){glue('as')}))
cores_gfp1 <- glue(cores_dir, "/lentiGFP1_COREs.bed") %>%
    read.table(., sep = "\t", col.names = c("chr", "start", "end")) %>%
    tibble::add_column(myid_core = sapply(c(1:nrow(.)), function(i){glue('c.gfp1.{i}')})) %>%
    tibble::add_column(coretype = sapply(c(1:nrow(.)), function(i){glue('gfp1')})) %>%
    tibble::add_column(coreset = sapply(c(1:nrow(.)), function(i){glue('gfp')}))
cores_gfp2 <- glue(cores_dir, "/lentiGFP2_COREs.bed") %>%
    read.table(., sep = "\t", col.names = c("chr", "start", "end")) %>%
    tibble::add_column(myid_core = sapply(c(1:nrow(.)), function(i){glue('c.gfp2.{i}')})) %>%
    tibble::add_column(coretype = sapply(c(1:nrow(.)), function(i){glue('gfp2')})) %>%
    tibble::add_column(coreset = sapply(c(1:nrow(.)), function(i){glue('gfp')}))
# dim(cores_as1); dim(cores_as2); dim(cores_gfp1); dim(cores_gfp2) # 442, 503, 220, 398 x 6

cores_all <- cores_as1 %>%
    rbind(., cores_as2) %>%
    rbind(., cores_gfp1) %>%
    rbind(., cores_gfp2) %>%
    tibble::add_column(width = (.)$end - (.)$start)
# dim(cores_all) # 1563 x 7

# COREs > 10Kb have id's starting with capital 'C'
for(c in 1:nrow(cores_all)){
    if(cores_all$width[c] > 10000){
        cores_all$myid_core[c] <- gsub('c', 'C', cores_all$myid_core[c])
    }
}

saveRDS(cores_all, file="./output/R_objects/cores_all.rds")
# END of DATA COREs


# DATA: promoters -----------------------------------------------------------------------------
# RNAseq genes --------------------------------------------------------------------------------
genes_dir <- "../data/RNAseq"
mygenes <- read.csv(glue(genes_dir, "/", "genes_list.csv"))[,3]
# length(mygenes) # 8678
# Gene symbols to Entrez ID  ------------------------------------------------------------------
genome_name <- 'rn6'
if (genome_name == 'rn6') {txdb_obj <- TxDb.Rnorvegicus.UCSC.rn6.refGene}

valid_chrs <- c(1:20, 'X') %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()

# Set active seqlevels/chrs
# To reset to orig seqlevels: seqlevels0(txdb_obj)
seqlevels(txdb_obj) <- valid_chrs

gene_symbol_to_entrez_map <- as.list(org.Rn.egALIAS2EG) 
# length(gene_symbol_to_entrez_map) # 59568 (previously 62195)
gene_symbol_to_entrez_map <- gene_symbol_to_entrez_map[!is.na(gene_symbol_to_entrez_map)] 
# length(gene_symbol_to_entrez_map) # 59568, same length --> all genes map to Entrez gene id's

mygenes_to_entrez <- gene_symbol_to_entrez_map[mygenes] 
# length(mygenes_to_entrez) # 8678 (out of 8678) --> all gene symbols have an Entrez ID

# unique symbols = 8668 = length 8678 - 10 NA's = 8668
entrez_full_df <- data.frame(symbol = character(), entrez = character())
symbols <- names(mygenes_to_entrez)
for (i in 1:length(mygenes_to_entrez)){
    symbol_add <- symbols[i]
    entrez1_add <- mygenes_to_entrez[[i]]
    if(length(entrez1_add) == 1){
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add)
    }
    else if (length(entrez1_add) == 2){
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add[1]
        )
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add[2]
        )
    }
    else if (length(entrez1_add) == 3){
        print(glue('entrez1_add 3 = {entrez1_add}'))
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add[1]
        )
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add[2]
        )
        entrez_full_df <- entrez_full_df %>% tibble::add_row(
            symbol = symbol_add, entrez = entrez1_add[3]
        )
    }
    else {
        # 10 NA keys with NULL element
        print(mygenes_to_entrez[i])
    }
} # end for
# dim(entrez_full_df) # 8816 x 2 (prev 8818)

# Add (my) unique gene id = g{1:8818}
entrez_full_df <- entrez_full_df %>% tibble::add_column(
    myid_gene = sapply(c(1:nrow(entrez_full_df)), function(i){glue('g{i}')}))
# dim(entrez_full_df) # 8816 x 3


# Promoters -----------------------------------------------------------------------------------

ratPromoters <- promoters(txdb_obj) # length(ratPromoters) # 18859
mykeys <- entrez_full_df$entrez # 8816
mytxnames <- AnnotationDbi::select(
    txdb_obj,
    keys = mykeys,
    columns = c("TXNAME", "TXSTRAND", "TXCHROM"),
    keytype = "GENEID")
# dim(mytxnames) # 9416 x 4 (prev 9418)
# 'select()' returned many:many mapping between keys and columns
# mytxnames$GENEID %>% unique() %>% length() # 8756 (prev 8757)
# mytxnames$TXNAME %>% unique() %>% length() # 9254 (prev 9255)

# myp <- ratPromoters[mytxnames$TXNAME %in% ratPromoters$tx_name, ] # 18763
# longer_list[longer_list %in% shorter_list] -> elements of longer list also in shorter list
mypromoters <- ratPromoters[ratPromoters$tx_name %in% mytxnames$TXNAME,] # GRanges object, length 9368
# Some entrez id's don't have a transcript txname so they're automatically omitted
# p <- mytxnames[!(mytxnames$TXNAME %in% ratPromoters$tx_name),]
# length(mypromoters) # 9368 # from 9416 mytxnames

mypranges <- mypromoters@ranges %>% as.data.frame() # IRanges object
mypseq <- data.frame(chr = mypromoters@seqnames)
promoters_all <- dplyr::bind_cols(mypseq, mypranges) # 9368 x 5

promoters_entrez_map <- mytxnames[match(promoters_all$names, mytxnames$TXNAME),] # 9368 x 4
entrez_symbol_map <- entrez_full_df[match(promoters_entrez_map$GENEID, entrez_full_df$entrez),] # 9368 x 3

promoters_all <- promoters_all %>% 
    tibble::add_column(entrez_id = promoters_entrez_map$GENEID) %>%
    tibble::add_column(gene_symbol = entrez_symbol_map$symbol) %>%
    tibble::add_column(myid_gene = entrez_symbol_map$myid_gene) %>%
    tibble::add_column(myid_promoter = sapply(c(1:nrow(promoters_all)), function(i){glue('p.{i}')}))
# dim(promoters_all) # 9368 x 9 
# promoters_all$entrez_id %>% unique %>% length() # 8709 

gr_promoters_all <- GRanges(
    seqnames = promoters_all$chr,
    ranges = IRanges(promoters_all$start, width = 2200),
    names = promoters_all$names,
    entrez.id = promoters_all$entrez_id,
    gene.symbol = promoters_all$gene_symbol,
    gene.id = promoters_all$myid_gene,
    promoter.id = promoters_all$myid_promoter
) # 9368 ranges and 5 metadata

saveRDS(promoters_all, glue("{rds_dir}/promoters_all.rds"))
saveRDS(gr_promoters_all, glue("{rds_dir}/gr_promoters_all.rds"))
# END of DATA promoters


# DATA: Loops ---------------------------------------------------------------------------------
# Called loops files --------------------------------------------------------------------------
loopsets_dir <- "./output/LOOPS/loops_sets"
path <- list(
    glue(loopsets_dir, '/', 'significant_loops_not_differential_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_gain_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_lost_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_up_notgain_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_down_notlost_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_gain_sigp_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_lost_sigp_50Kb-2Mb.bed')
)
names(path) <- c('sig_notdiff_all', 'diff_gain', 'diff_lost', 'diff_up_notgain', 'diff_down_notlost',
                 'true_diff_gain', 'true_diff_lost')

loops_true_diff_gain <- read.csv(path$true_diff_gain, sep = "\t") %>% # 104 x 5
    tibble::add_column(loopset = 'Ldg') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ldg.true.{i}')})) 

loops_true_diff_lost <- read.csv(path$true_diff_lost, sep = "\t") %>% # 441 x 5
    tibble::add_column(loopset = 'Ldl') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ldl.true.{i}')})) 

loops_diff_gain <- read.csv(path$diff_gain, sep = "\t") %>% # 8665 
    tibble::add_column(loopset = 'Ldg') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ldg.{i}')})) 

loops_diff_lost <- read.csv(path$diff_lost, sep = "\t") %>% # 21240
    tibble::add_column(loopset = 'Ldl') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ldl.{i}')})) 

loops_diff_up_notgain <- read.csv(path$diff_up_notgain, sep = "\t") %>% # 42705
    tibble::add_column(loopset = 'Ldung') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ldung.{i}')})) 

loops_diff_down_notlost <- read.csv(path$diff_down_notlost, sep = "\t") %>% # 61893
    tibble::add_column(loopset = 'Lddnl') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Lddnl.{i}')})) 

loops_sig_notdiff_all <- read.csv(path$sig_notdiff_all, sep = "\t") %>% # 53675
    tibble::add_column(loopset = 'Ls') %>%
    tibble::add_column(myid_loop = sapply(c(1:nrow(.)), function(i){glue('Ls.{i}')})) 

# w/ true diff gain and lost: 188723 x 5
loops_all <- loops_sig_notdiff_all %>% 
    dplyr::full_join(., loops_diff_gain) %>% 
    dplyr::full_join(., loops_diff_lost) %>%
    dplyr::full_join(., loops_diff_up_notgain) %>%
    dplyr::full_join(., loops_diff_down_notlost) %>%
    dplyr::full_join(., loops_true_diff_gain) %>%
    dplyr::full_join(., loops_true_diff_lost)


loops_diff_gain_remove <- character() # 104 at end of for loop, ~5 sec
for(i in 1:nrow(loops_diff_gain)){
    for(j in 1:nrow(loops_true_diff_gain)){
        if( (loops_diff_gain$chr[i] == loops_true_diff_gain$chr[j]) &&
            (loops_diff_gain$startI[i] == loops_true_diff_gain$startI[j]) &&
            (loops_diff_gain$startJ[i] == loops_true_diff_gain$startJ[j])) {
            loops_diff_gain_remove <- loops_diff_gain_remove %>% append(loops_diff_gain$myid_loop[i])
        }
    }
}

loops_diff_lost_remove <- character() # 441 at end of for loop, ~10 sec
for(i in 1:nrow(loops_diff_lost)){
    for(j in 1:nrow(loops_true_diff_lost)){
        if( (loops_diff_lost$chr[i] == loops_true_diff_lost$chr[j]) &&
            (loops_diff_lost$startI[i] == loops_true_diff_lost$startI[j]) &&
            (loops_diff_lost$startJ[i] == loops_true_diff_lost$startJ[j])) {
            loops_diff_lost_remove <- loops_diff_lost_remove %>% append(loops_diff_lost$myid_loop[i])
        }
    }
}

# 188178 x 5
loops_all <- loops_all %>% 
    dplyr::filter(!myid_loop %in% c(loops_diff_gain_remove, loops_diff_lost_remove)) %>%
    dplyr::mutate(type = dplyr::case_when(
        stringr::str_detect(.$myid_loop, "true") == FALSE ~ "static",
        stringr::str_detect(.$myid_loop, "Ldg.true") == TRUE ~ "gain",
        stringr::str_detect(.$myid_loop, "Ldl.true") == TRUE ~ "lost"
    ))

# Check:
# loops_all %>% dplyr::filter(type == "static") %>% nrow() # 187633
# loops_all %>% dplyr::filter(type == "gain") %>% nrow() # 104
# loops_all %>% dplyr::filter(type == "lost") %>% nrow() # 441
# 187633 + 104 + 441 = 188178
saveRDS(loops_all, file=glue("{rds_dir}/loops_all.rds"))

loops_type_count_all_chrs <- tibble::tibble(
    chr = character(),
    total.loops.count = numeric(),
    static.loops.count = numeric(),
    gain.loops.count = numeric(),
    lost.loops.count = numeric())

valid_chrs <- c(1:20, "X") %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
for(chrom in valid_chrs){
    loops_type_count_all_chrs <- loops_type_count_all_chrs %>% tibble::add_row(
        chr = chrom,
        total.loops.count = loops_all %>% dplyr::filter(chr == chrom) %>% nrow(),
        static.loops.count = loops_all %>% dplyr::filter(chr == chrom) %>% dplyr::filter(type == "static") %>% nrow(),
        gain.loops.count = loops_all %>% dplyr::filter(chr == chrom) %>% dplyr::filter(type == "gain") %>% nrow(),
        lost.loops.count = loops_all %>% dplyr::filter(chr == chrom) %>% dplyr::filter(type == "lost") %>% nrow()
    )      
}

# readr::write_tsv(loops_type_count_all_chrs, file="./output/loops_data_processed/loops_type_count_all_chrs.txt")
saveRDS(loops_type_count_all_chrs, file=glue("{rds_dir}/loops_type_count_all_chrs.rds"))
# END of DATA loops


# DATA: TADs ----------------------------------------------------------------------------------
tad_dir <- "./output/TADS"
tad_cols <- c("startpos", "endpos", "TADlevel", "TADmean", "TADscore")
valid_chroms <- c(1:20, 'X') %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
ontad_as_tads <- list() 
ontad_gfp_tads <- list() 

for(chrom in valid_chroms) {
    ontad_as_tads[[chrom]] <- readr::read_tsv(
        glue("{tad_dir}/as_10000_merged_iced_all_intra_{chrom}.ontad.tad"), col_names=tad_cols) %>%
        tibble::add_column(chrom = chrom, .before=1) %>%
        tibble::add_column(condition = "as", .before=1)
    ontad_gfp_tads[[chrom]] <- readr::read_tsv(
        glue("{tad_dir}/gfp_10000_merged_iced_all_intra_{chrom}.ontad.tad"), col_names=tad_cols) %>%
        tibble::add_column(chrom = chrom, .before=1) %>%
        tibble::add_column(condition = "gfp", .before=1)
}
saveRDS(ontad_as_tads, file="./output/R_objects/ontad_as_tads.rds")
saveRDS(ontad_gfp_tads, file="./output/R_objects/ontad_gfp_tads.rds")


as_tads_all <- dplyr::bind_rows(ontad_as_tads) # 16911 x 7
gfp_tads_all <- dplyr::bind_rows(ontad_gfp_tads) # 16848 x 7
tads_all <- dplyr::bind_rows(as_tads_all, gfp_tads_all) # 33759 x 7


# DATA: TF footprinting -----------------------------------------------------------------------
# Preprocess TF footprinting data and save as R object ~1 min runtime

indir_tf <- './data/TF/footprints/'
tf_files <- list.files(indir_tf)

cols.fp <- c('chr','start','end','motif_name','V5','strand','peak_chr','peak_start','peak_end','binding_score')
# cols.fp <- c('chr','start','end','motif_name','strand','chr','V7','V8','V9','V10')
# chr, start of motif, end of motif, name of motif, strand, 
# chr, peak start, peak end, footprint score (binding strength)


a.fp.data <- glue('{indir_tf}/{tf_files[1]}') %>% read.table(., sep='\t', col.names=cols.fp) # 299270 x 10
b.fp.data <- glue('{indir_tf}/{tf_files[2]}') %>% read.table(., sep='\t', col.names=cols.fp) # 363518 x 10
c.fp.data <- glue('{indir_tf}/{tf_files[3]}') %>% read.table(., sep='\t', col.names=cols.fp) # 318314 x 10
d.fp.data <- glue('{indir_tf}/{tf_files[4]}') %>% read.table(., sep='\t', col.names=cols.fp) # 241371 x 10

# prep tf data (same in all files) line 39 - 93 same in all files
prep_fp <- function(fp, chrm){
    # For each TF, all chr for AS then all chr for GFP
    fp <- fp %>%
        dplyr::select('chr','start','end','motif_name','binding_score',
                      'V5','strand','peak_chr','peak_start','peak_end') %>%
        tibble::add_column(row_id = c(1:nrow(.))) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames('row_id')
    
    motifs <- fp$motif_name %>% unique() # 668
    fp.chrm <- fp %>% dplyr::filter(chr == chrm) # 75597 x 10
    
    fp.chrm.with.cond <- data.frame()
    for(m in 1:length(motifs)){
        #print(glue::glue('m = {m}'))
        fp.chrm.motif <- fp.chrm %>% dplyr::filter(motif_name == motifs[m]) 
        cond <- vector(mode="character", length=nrow(fp.chrm.motif))
        gfp_start_idx <- NULL
        if(nrow(fp.chrm.motif) != 0){
            for (i in 1:nrow(fp.chrm.motif)){
                
                if(!is.null(gfp_start_idx)) {
                    cond[gfp_start_idx:nrow(fp.chrm.motif)] <- 'GFP'
                    break
                } else if(i == 1) { 
                    cond[i] <- 'AS'
                } else if(fp.chrm.motif$start[i] > fp.chrm.motif$start[i-1]){ 
                    cond[i] <- 'AS'
                } else {
                    gfp_start_idx <- i
                }
            }
        }
        
        fp.chrm.motif <- fp.chrm.motif %>% tibble::add_column(condition = cond, .after=5)
        fp.chrm.with.cond <- dplyr::bind_rows(fp.chrm.with.cond,fp.chrm.motif)
    }
    fp.chrm.with.cond <- fp.chrm.with.cond %>%
        tibble::add_column(length = (.)$end - (.)$start, .after=6)
    return(fp.chrm.with.cond)
}

# Create and save tfs_all as tfs_all.rds
make_empty_tfs_df <- function(){
    data.frame(chr = character(),
               start = numeric(),
               end = numeric(),
               motif_name = character(),
               V5 = numeric(),
               strand = character(),
               peak_chr = character(),
               peak_start = numeric(),
               peak_end = numeric(),
               binding_score = numeric())
}



valid_chrs <- c(1:20, "X") %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
a.fp.all.chroms <- make_empty_tfs_df()
b.fp.all.chroms <- make_empty_tfs_df()
c.fp.all.chroms <- make_empty_tfs_df()
d.fp.all.chroms <- make_empty_tfs_df()

# ~50 sec
for(chrom in valid_chrs){
    a.fp.all.chroms <- a.fp.all.chroms %>%
        dplyr::bind_rows(prep_fp(fp = a.fp.data, chrm=chrom))
    b.fp.all.chroms <- b.fp.all.chroms %>%
        dplyr::bind_rows(prep_fp(fp = b.fp.data, chrm=chrom))
    c.fp.all.chroms <- c.fp.all.chroms %>%
        dplyr::bind_rows(prep_fp(fp = c.fp.data, chrm=chrom))
    d.fp.all.chroms <- d.fp.all.chroms %>%
        dplyr::bind_rows(prep_fp(fp = d.fp.data, chrm=chrom))
    #a.fp.chrom <- prep_fp(fp = a.fp.data, chrm=chr) 
    #b.fp.chrom <- prep_fp(fp = b.fp.data, chrm=chr) 
    #c.fp.chrom <- prep_fp(fp = c.fp.data, chrm=chr) 
    #d.fp.chrom <- prep_fp(fp = d.fp.data, chrm=chr) 
}

# Check: dims for *.fp.all.chroms match their *.fp.data 
# a.fp.all.chroms # 299270 x 12
# b.fp.all.chroms # 363518 x 12
# c.fp.all.chroms # 318314 x 12
# d.fp.all.chroms # 241371 x 12

# tfs_all <- make_empty_tfs_df()
tfs_all <- a.fp.all.chroms %>%
    dplyr::bind_rows(b.fp.all.chroms) %>%
    dplyr::bind_rows(c.fp.all.chroms) %>%
    dplyr::bind_rows(d.fp.all.chroms) %>%
    tibble::add_column(tf = .$motif_name %>%
                           stringr::str_split(., '_') %>%
                           sapply(., FUN = function(x){x[1]}) %>%
                           stringr::str_split(., pattern = 'var') %>%
                           sapply(., FUN = function(x){x[1]}))
# 1,222,473 x 12 = 241371 + 318314 + 363518 + 299270 = 1054936 + 167537

# ~6 sec
tfs_all_uniq <- tfs_all %>% unique() # 1,054,936 x 13
# tfs_all_dup <- tfs_all[duplicated(tfs_all),] # 167537 x 13
# tfs_all_dup_uniq <- tfs_all_dup %>% unique() # 131995 x 13

# GRanges object with 1054936 ranges and 7 metadata columns:
gr_tfs_all_uniq <- GRanges(
    seqnames = tfs_all_uniq$chr,
    ranges = IRanges(start=tfs_all_uniq$start, end = tfs_all_uniq$end),
    motif_name = tfs_all_uniq$motif_name,
    condition = tfs_all_uniq$condition,
    length = tfs_all_uniq$length,
    peak_start =  tfs_all_uniq$peak_start,
    peak_end = tfs_all_uniq$peak_end,
    binding_score = tfs_all_uniq$binding_score,
    tf = tfs_all_uniq$tf
) 


# Save as R objects
saveRDS(tfs_all, glue("{rds_dir}/tfs_all.rds"))
saveRDS(tfs_all_uniq, glue("{rds_dir}/tfs_all_uniq.rds"))
saveRDS(gr_tfs_all_uniq, glue("{rds_dir}/gr_tfs_all_uniq.rds"))
# END of DATA TF

# DATA: A/B compartments ----------------------------------------------------------------------
indir_eigen <- "./output/COMPARTMENTS"
valid_chrs <- c(1:20, "X") %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
binsize = 100000

prep_eigens <- function(cond){
    eigens <- list() 
    for (chrom in valid_chrs){
        
        # Adapted from HiCDCPlus::extract_hic_eigenvectors
        eigenfile <- glue("{indir_eigen}/{cond}/{cond}_all_intra_eigen_{chrom}.txt")
        out.df <- data.table::fread(eigenfile)
        out.df <- out.df %>% 
            dplyr::mutate(chr=chrom, 
                          start=binsize*seq(0,nrow(out.df)-1,1), 
                          end=binsize*seq(1,nrow(out.df),1)) %>%
            dplyr::rename("score"=.data$V1) %>%
            dplyr::select(.data$chr, .data$start, .data$end, .data$score)
        out.df[1,"start"]=1
        
        eigens[[chrom]] <- out.df
    }
    return(eigens)
}
as_eigens <- prep_eigens(cond = "as")
gfp_eigens <- prep_eigens(cond = "gfp")

saveRDS(as_eigens, file=glue("{indir_eigen}/{cond}/as_eigens.rds"))
saveRDS(gfp_eigens, file=glue("{indir_eigen}/{cond}/gfp_eigens.rds"))
saveRDS(as_eigens, file=glue("{rds_dir}/as_eigens.rds"))
saveRDS(gfp_eigens, file=glue("{rds_dir}/gfp_eigens.rds"))

as_eigens_chr5 <- as_eigens$chr5
gfp_eigens_chr5 <- gfp_eigens$chr5
# readr::write_tsv(as_eigens_chr5, file="./output/processed_files/as_eigens_chr5.txt")
# readr::write_tsv(gfp_eigens_chr5, file="./output/processed_files/gfp_eigens_chr5.txt")
# END of DATA A/B compartments


