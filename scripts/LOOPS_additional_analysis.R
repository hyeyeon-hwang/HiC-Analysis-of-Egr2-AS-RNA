# Last edited: 6 Apr 2022 

# Additional analysis

# Load R objects
genes.a1core.a2prom.uniq.fin <- readRDS("./output/R_objects/genes.a1core.a2prom.uniq.fin.rds")
genes.a1prom.a2core.uniq.fin <- readRDS("./output/R_objects/genes.a1prom.a2core.uniq.fin.rds")

# What is number of loops associated with genes? --------------------------------------------
genes.a1prom.a2core.uniq.fin$num_uniq_loops %>% sum() # 101
genes.a1core.a2prom.uniq.fin$num_uniq_loops %>% sum() # 126
# total = 227
a1prom.a2core.loops$myid_loop %>% length() %>% sum() # 90, all unique
a1core.a2prom.loops$myid_loop %>% length() %>% sum() # 111, all unique

# What is number of genes that mapped to anchor
# in one anchor and core in another anchor?
a1core.a2prom.loops$gene_symbol %>% unlist() %>% unique() %>% length() # 86
a1prom.a2core.loops$gene_symbol %>% unlist() %>% unique() %>% length() # 78


# What number of genes mapped to any anchor?
# there doesn't have to be a core in the other anchor

# uci = unique core id
uci <- table(genes.a1core.a2prom.uniq.fin$uniq_core_ids %>% unlist())
uci['C.as2.3'] # 1
uci['C.gfp1.11'] # 1
uci['C.gfp2.9'] # 1

uci2 <- table(genes.a1prom.a2core.uniq.fin$uniq_core_ids %>% unlist())
uci2['C.as2.3'] # NA 
uci2['C.gfp1.11'] # NA 
uci2['C.gfp2.9'] # NA



# Which genes overlap? ------------------------------------------------------------------------
all_genes <- c(genes.a1prom.a2core.uniq.fin$gene_symbol, genes.a1core.a2prom.uniq.fin$gene_symbol)
overlapping_genes <- all_genes[duplicated(all_genes)] # "Trmt11" "Hint3"


# What percent of COREs match to loop anchors? -------------------------------------------------

genes.a1core.a2prom.uniq.fin$uniq_core_ids %>% unlist() %>% length() # 133
ucores.a1c.a2p <- genes.a1core.a2prom.uniq.fin$uniq_core_ids %>% unlist() %>% unique() # 80

genes.a1prom.a2core.uniq.fin$uniq_core_ids %>% unlist() %>% length() # 144
ucores.a1p.a2c <- genes.a1prom.a2core.uniq.fin$uniq_core_ids %>% unlist() %>% unique() # 90

ucores_all <- c(ucores.a1c.a2p, ucores.a1p.a2c) # 170
ucores_all_uniq <- ucores_all %>% unique() # 157

regmatches(ucores_all_uniq, gregexpr('as',ucores_all_uniq )) %>% lengths() %>% sum() # 92
regmatches(ucores_all_uniq, gregexpr('gfp',ucores_all_uniq )) %>% lengths() %>% sum() # 65      

# Does mTOR promoter overlap with any loop anchors? -------------------------------

# mTOR promoter id = p.3271
mtor.prom.id <- genes.a1core.a2prom.uniq.fin %>% 
    dplyr::filter(gene_symbol == 'Mtor') %>% 
    dplyr::select(uniq_prom_ids) %>% 
    unlist()
mtor.prom.row <- promoters_all %>% dplyr::filter(myid_promoter == mtor.prom.id)
gr.mtor.prom <- GRanges(seqnames = mtor.prom.row$chr, 
                        ranges = IRanges(mtor.prom.row$start, width = 2200))

# is mTOR promoter p.3271 in any of the anchors of the 192 unique loop? NO.
a1prom.a2core.loops %>% dplyr::filter(grepl("p.3271", proms_per_loop_ids %>% unlist()))
a1prom.a2core.loops$proms_per_loop_ids %>% unlist() %>% grepl("p.3271", .) # all FALSE
a1prom.a2core.loops$proms_per_loop_ids %>% unlist() %>% grepl("p.3271", .) %>% sum() # 0

a1core.a2prom.loops$proms_per_loop_ids %>% unlist() %>% grepl("p.3271", .) # 3 TRUE -> mTOR anchors
a1core.a2prom.loops$proms_per_loop_ids %>% unlist() %>% grepl("p.3271", .) %>% sum() # 3



# Is mTOR promoter p.3271 in any loop anchors from the sig/diff loops? ------------------------

loops_all_chr5 <- loops_all %>% dplyr::filter(chr == 'chr5') # 8280 x 5
gr.loops.a1.chr5 <- GRanges(seqnames = loops_all_chr5$chr, 
                            ranges = IRanges(loops_all_chr5$startI, width = 10000),
                            myid_loop = loops_all_chr5$myid_loop)
gr.loops.a2.chr5 <- GRanges(seqnames = loops_all_chr5$chr, 
                            ranges = IRanges(loops_all_chr5$startJ, width = 10000),
                            myid_loop = loops_all_chr5$myid_loop)

mtor.prom.hits.a1 <- GenomicRanges::findOverlaps(
    gr.mtor.prom, gr.loops.a1.chr5,
    type='within', # promoter 2200 within loop anchor1 10000
    select='all', ignore.strand=TRUE) # 0 hits
mtor.prom.hits.a2 <- GenomicRanges::findOverlaps(
    gr.mtor.prom, gr.loops.a2.chr5,
    type='within', # promoter 2200 within loop anchor2 10000
    select='all', ignore.strand=TRUE) # 4 hits

mtor.prom.in.any.loopa2 <- gr.loops.a2.chr5[mtor.prom.hits.a2@to,]
mtor.prom.in.any.loopa2.loopids <- mtor.prom.in.any.loopa2$myid_loop
# GRanges object with 4 ranges and 1 metadata column:
# seqnames              ranges strand |   myid_loop
# <Rle>           <IRanges>  <Rle> | <character>
# [1]     chr5 165260000-165269999      * |    Ldg.3368
# [2]     chr5 165260000-165269999      * |    Ldl.7951
# [3]     chr5 165260000-165269999      * | Ldung.16390
# [4]     chr5 165260000-165269999      * | Lddnl.23589


# END of analysis -----------------------------------------------------------------------------
