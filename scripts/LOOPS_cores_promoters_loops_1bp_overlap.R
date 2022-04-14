# Last edited: 6 Apr 2022 

# Find COREs-promoters and promoters-COREs loops (overlap of at least 1 bp)
# Only difference with original script: Changed type='within' --> 'any' AND add maxgap = -1

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

# Output directories --------------------------------------------------------------------------
# outdir_main <- glue('./output/LOOPS/{format(Sys.time(), "%b%d")}')
outdir_main <- glue('./output/LOOPS')
if (!dir.exists(outdir_main)){dir.create(outdir_main)} else {print(glue("{outdir_main} exists!"))}

# outdir_pcloops <- glue('{outdir_main}/CORES_promoters_loops')
outdir_pcloops <- glue('{outdir_main}/CORES_promoters_loops_1bp_overlap')
if (!dir.exists(outdir_pcloops)){ dir.create(outdir_pcloops)} else {print(glue("{outdir_pcloops} exists!"))}
dir.create(glue('{outdir_pcloops}/a1core_a2prom'))
dir.create(glue('{outdir_pcloops}/a1prom_a2core'))

# Called loops files --------------------------------------------------------------------------

loopsets_dir <- "./output/LOOPS/loops_sets"
path <- list(
    glue(loopsets_dir, '/', 'significant_loops_not_differential_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_gain_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_lost_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_up_notgain_50Kb-2Mb.bed'),
    glue(loopsets_dir, '/', 'differential_loops_down_notlost_50Kb-2Mb.bed')
)
names(path) <- c('sig_notdiff_all', 'diff_gain', 'diff_lost', 'diff_up_notgain', 'diff_down_notlost')

# RNAseq genes --------------------------------------------------------------------------------

genes_dir <- "./data/RNAseq"
mygenes <- read.csv(glue(genes_dir, "/", "genes_list.csv"))[,3]
# length(mygenes) # 8678

# COREs ---------------------------------------------------------------------------------------

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
# dim(entrez_full_df) # 8816 x 2 

# Add (my) unique gene id = g{1:8816}
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
# dim(mytxnames) # 9416 x 4 
# 'select()' returned many:many mapping between keys and columns
# mytxnames$GENEID %>% unique() %>% length() # 8756 
# mytxnames$TXNAME %>% unique() %>% length() # 9254 

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

# Find promoters-cores and cores-promoters interactions --------------------------------------------

gr_promoters <- GRanges(
    seqnames = promoters_all$chr,
    ranges = IRanges(promoters_all$start, width = 2200),
    names = promoters_all$names,
    entrez_id = promoters_all$entrez_id,
    gene_symbol = promoters_all$gene_symbol,
    myid_gene = promoters_all$myid_gene,
    myid_promoter = promoters_all$myid_promoter
)

gr_cores <- GRanges(
    seqnames = cores_all$chr,
    ranges = IRanges(start=cores_all$start, end=cores_all$end),
    myid_core = cores_all$myid_core,
    mywidth = cores_all$width,
    coreset = cores_all$coreset)


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

# 188178 x 5
loops_all <- loops_sig_notdiff_all %>% 
    dplyr::full_join(., loops_diff_gain) %>% 
    dplyr::full_join(., loops_diff_lost) %>%
    dplyr::full_join(., loops_diff_up_notgain) %>%
    dplyr::full_join(., loops_diff_down_notlost) 

loops <- loops_all

gr_loops_a1 <- GRanges(seqnames = loops$chr, 
                       ranges = IRanges(loops$startI, width = 10000),
                       myid_loop = loops$myid_loop)
gr_loops_a2 <- GRanges(seqnames = loops$chr, 
                       ranges = IRanges(loops$startJ, width = 10000),
                       myid_loop = loops$myid_loop)


# Change type='within' --> 'any' AND add maxgap = -1 ------------------------------------------
indices.p.a1 <- GenomicRanges::findOverlaps(
    gr_promoters, gr_loops_a1,
    type='any', # promoters 2200 within loop anchor1 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
indices.p.a2 <- GenomicRanges::findOverlaps(
    gr_promoters, gr_loops_a2,
    type='any', # promoters 2200 within loop anchor2 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
# Changed type='within' --> 'any' AND add maxgap = -1 -----------------------------------------

prom.in.Loops.a1.promoters <- promoters_all[indices.p.a1@from,]
prom.in.Loops.a1.promoters.uniq <- prom.in.Loops.a1.promoters %>%
    unique() %>%
    tibble::add_column(loops_a1_per_prom = table(prom.in.Loops.a1.promoters$myid_promoter)[(.)$myid_promoter])

prom.in.Loops.a1.loops <- loops[indices.p.a1@to,] %>% 
    tibble::add_column(myid_promoter = prom.in.Loops.a1.promoters$myid_promoter) %>%
    #tibble::add_column(myid_gene = prom.in.Loops.a1.promoters$myid_gene) %>%
    tibble::add_column(gene_symbol = prom.in.Loops.a1.promoters$gene_symbol) 

prom_ids <- list()
gene_symbols <- list()
for(i in 1:nrow(prom.in.Loops.a1.loops)){
    prom_ids[i] <- list(prom.in.Loops.a1.loops[prom.in.Loops.a1.loops$myid_loop == prom.in.Loops.a1.loops$myid_loop[i],'myid_promoter'])
    gene_symbols[i] <- list(prom.in.Loops.a1.loops[prom.in.Loops.a1.loops$myid_loop == prom.in.Loops.a1.loops$myid_loop[i],'gene_symbol'])}

prom.in.Loops.a1.loops.uniq <- prom.in.Loops.a1.loops %>%
    tibble::add_column(proms_per_loop_a1 = table(prom.in.Loops.a1.loops$myid_loop)[(.)$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(proms_per_loop_ids = prom_ids) %>% 
    dplyr::select(., -gene_symbol) %>%
    tibble::add_column(gene_symbol = gene_symbols) %>%
    dplyr::select(., -myid_promoter) %>%
    unique()


prom.in.Loops.a2.promoters <- promoters_all[indices.p.a2@from,]
prom.in.Loops.a2.promoters.uniq <- prom.in.Loops.a2.promoters %>%
    unique() %>%
    tibble::add_column(loops_a2_per_prom = table(prom.in.Loops.a2.promoters$myid_promoter)[(.)$myid_promoter]) 


prom.in.Loops.a2.loops <- loops[indices.p.a2@to,] %>%
    tibble::add_column(myid_promoter = prom.in.Loops.a2.promoters$myid_promoter) %>%
    tibble::add_column(gene_symbol = prom.in.Loops.a2.promoters$gene_symbol) 

prom_ids <- list()
gene_symbols <- list()

for(i in 1:nrow(prom.in.Loops.a2.loops)){
    prom_ids[i] <- list(prom.in.Loops.a2.loops[prom.in.Loops.a2.loops$myid_loop == prom.in.Loops.a2.loops$myid_loop[i],'myid_promoter'])
    gene_symbols[i] <- list(
        prom.in.Loops.a2.loops[prom.in.Loops.a2.loops$myid_loop == prom.in.Loops.a2.loops$myid_loop[i],'gene_symbol'])
}

prom.in.Loops.a2.loops.uniq <- prom.in.Loops.a2.loops %>%
    tibble::add_column(proms_per_loop_a2 = table(prom.in.Loops.a2.loops$myid_loop)[.$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(proms_per_loop_ids = prom_ids) %>% 
    dplyr::select(., -gene_symbol) %>%
    tibble::add_column(gene_symbol = gene_symbols) %>%
    dplyr::select(., -myid_promoter) %>%
    unique()


# any promoter and any core 1bp overlap 
# Change type='within' --> 'any' AND add maxgap = -1 ------------------------------------------
pc_any <- GenomicRanges::findOverlaps(
    gr_promoters, gr_cores,
    type='any',
    maxgap = -1,
    select='all')
# Changed type='within' --> 'any' AND add maxgap = -1 -----------------------------------------

pc_any_proms <- gr_promoters[pc_any@from,] # 603 #  pc_any_proms %>% unique() %>% length() # 319
pc_any_cores <- gr_cores[pc_any@to,] # 603 # pc_any_cores %>% unique() %>% length() # 317

# prom = 2200, core > 2200
# 603 x 3
any_prom_in_any_core_gr2200 <- data.frame(gene_symbol = pc_any_proms$gene_symbol,
                                          prom_id = pc_any_proms$myid_promoter,
                                          core_id = pc_any_cores$myid_core)

#promoters_all[promoters_all$myid_promoter == 'p.3270',]
#promoters_all[promoters_all$myid_promoter == 'p.3573',]
#promoters_all[promoters_all$myid_promoter == 'p.3574',]

promoters_all[promoters_all$gene_symbol == 'Mad2l2',]
promoters_all[promoters_all$gene_symbol == 'Agtrap',]
promoters_all[promoters_all$gene_symbol == 'Fbxo6',]

any_prom_in_any_core_gr2200[any_prom_in_any_core_gr2200$core_id == 'C.as2.3',]
any_prom_in_any_core_gr2200[any_prom_in_any_core_gr2200$core_id == 'C.gfp1.11',]
any_prom_in_any_core_gr2200[any_prom_in_any_core_gr2200$core_id == 'C.gfp2.9',]


# find cores in loops a1
# find cores in loops a2
# intersect common loops 

# intersect prom loop a1 + core loops a2
# intersect core loop a1 + prom loops a2
# intersect core loop a1 + core loop a2

# Changed type='within' --> 'any' AND add maxgap = -1 ------------------------------------------
indices.c.a1 <- GenomicRanges::findOverlaps(
    gr_cores, gr_loops_a1,
    type='any', # promoters 2200 within loop anchor1 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
indices.c.a2 <- GenomicRanges::findOverlaps(
    gr_cores, gr_loops_a2,
    type='any', # promoters 2200 within loop anchor1 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
# Changed type='within' --> 'any' AND add maxgap = -1 -----------------------------------------

# 59 -> 35 unique cores
core.in.Loops.a1.cores <- cores_all[indices.c.a1@from,] # 59 x 4
core.in.Loops.a1.cores.uniq <- core.in.Loops.a1.cores %>%
    unique() %>% # 35 x 4
    tibble::add_column(loops_a1_per_core = table(core.in.Loops.a1.cores$myid_core)[(.)$myid_core])
# a1 coord in prom.in.Loops.a1.loops is same for different (adjacent) loops, so the same core will map
#    to multiple loops with the same a1 coord 

# 59 -> 40 unique loops
#prom.in.Loops.a1.loopsa1 <- gr_loops_a1[indices.p.a1@to]
core.in.Loops.a1.loops <- loops[indices.c.a1@to,] %>%
    tibble::add_column(myid_core = core.in.Loops.a1.cores$myid_core) %>%
    tibble::add_column(coreset = core.in.Loops.a1.cores$coreset)

core_ids <- list()
core_sets <- list()
as_cores_per_loop <- numeric()
gfp_cores_per_loop <- numeric()
for(i in 1:nrow(core.in.Loops.a1.loops)){
    core_ids[i] <- list(core.in.Loops.a1.loops[core.in.Loops.a1.loops$myid_loop == core.in.Loops.a1.loops$myid_loop[i],'myid_core'])
    core_sets[i] <- list(core.in.Loops.a1.loops[core.in.Loops.a1.loops$myid_loop == core.in.Loops.a1.loops$myid_loop[i],'coreset'])
    print(core_ids[i])
    print(core_sets[i])
    as_cores_count <- 0
    gfp_cores_count <- 0
    for(j in 1:length(core_sets[[i]])){
        if(core_sets[[i]][j] == 'as'){as_cores_count = as_cores_count + 1}
        else if(core_sets[[i]][j] == 'gfp'){gfp_cores_count = gfp_cores_count + 1}
    }
    as_cores_per_loop[i] <- as_cores_count
    gfp_cores_per_loop[i] <- gfp_cores_count
}

core.in.Loops.a1.loops.uniq <- core.in.Loops.a1.loops %>%
    tibble::add_column(cores_per_loop_a1 = table(core.in.Loops.a1.loops$myid_loop)[(.)$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(cores_per_loop_ids = core_ids) %>% 
    tibble::add_column(cores_per_loop_sets = core_sets) %>%
    tibble::add_column(as_cores_per_loop = as_cores_per_loop) %>%
    tibble::add_column(gfp_cores_per_loop = gfp_cores_per_loop) %>%
    dplyr::select(., -myid_core) %>%
    dplyr::select(., -coreset) %>%
    unique()
dim(core.in.Loops.a1.loops.uniq) # 1386 x 10

# In core.in.Loops.a1.loops, anchor 1 coordinates can be the same (with different anchor 2 coordinates) 
# That's why multiple loops with the same anchor 1 coordinates (but different anchor 2 coordinates) can be overlap
# with the same cores --> multiple loops per core

# Multiple cores that are < 10000 in width can map to the same loop (anchor 1). 

# For example,
# > core.in.Loops.a1.loops[c(1,17),]
# chr   startI   startJ myid_loop
# 5875   chr11 74310000 74740000   Ldg5875
# 5875.1 chr11 74310000 74740000   Ldg5875

# > core.in.Loops.a1.cores[c(1,17),]
# chr    start      end myid_core width
# 75  chr11 74311692 74316651  c.as1.75  4959
# 891 chr11 74311637 74316640  c.as2.89  5003

# Anchor 2
core.in.Loops.a2.cores <- cores_all[indices.c.a2@from,] # 87 x 5
core.in.Loops.a2.cores.uniq <- core.in.Loops.a2.cores %>% # 62 x 6
    unique() %>% 
    tibble::add_column(loops_a2_per_core = table(core.in.Loops.a2.cores$myid_core)[(.)$myid_core])

core.in.Loops.a2.loops <- loops[indices.c.a2@to,] %>% # 87 x 4
    tibble::add_column(myid_core = core.in.Loops.a2.cores$myid_core) %>%
    tibble::add_column(coreset = core.in.Loops.a2.cores$coreset)

core_ids <- list()
core_sets <- list()
as_cores_per_loop <- numeric()
gfp_cores_per_loop <- numeric()
for(i in 1:nrow(core.in.Loops.a2.loops)){
    core_ids[i] <- list(core.in.Loops.a2.loops[core.in.Loops.a2.loops$myid_loop == core.in.Loops.a2.loops$myid_loop[i],'myid_core'])
    core_sets[i] <- list(core.in.Loops.a2.loops[core.in.Loops.a2.loops$myid_loop == core.in.Loops.a2.loops$myid_loop[i],'coreset'])
    as_cores_count <- 0
    gfp_cores_count <- 0
    for(j in 1:length(core_sets[[i]])){
        if(core_sets[[i]][j] == 'as1' | core_sets[[i]][j] == 'as2'){as_cores_count = as_cores_count + 1}
        else if(core_sets[[i]][j] == 'gfp1' | core_sets[[i]][j] == 'gfp2'){gfp_cores_count = gfp_cores_count + 1}
    }
    as_cores_per_loop[i] <- as_cores_count
    gfp_cores_per_loop[i] <- gfp_cores_count
}
core.in.Loops.a2.loops.uniq <- core.in.Loops.a2.loops %>% # 50 x 5
    tibble::add_column(cores_per_loop_a2 = table(core.in.Loops.a2.loops$myid_loop)[(.)$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(cores_per_loop_ids = core_ids) %>%
    tibble::add_column(cores_per_loop_sets = core_sets) %>%
    tibble::add_column(as_cores_per_loop = as_cores_per_loop) %>%
    tibble::add_column(gfp_cores_per_loop = gfp_cores_per_loop) %>%
    dplyr::select(., -myid_core) %>%
    dplyr::select(., -coreset) %>%
    unique()
dim(core.in.Loops.a2.loops.uniq) # 1185 x 10


# Change type='within' --> 'any' AND add maxgap = -1 -------------------------------------------
indices.c.a1.v2 <- GenomicRanges::findOverlaps(
    gr_loops_a1, gr_cores,
    type='any', # promoters 2200 within loop anchor1 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
indices.c.a2.v2 <- GenomicRanges::findOverlaps(
    gr_loops_a2, gr_cores,
    type='any', # promoters 2200 within loop anchor1 10000
    maxgap = -1,
    select='all', ignore.strand=TRUE)
# Changed type='within' --> 'any' AND add maxgap = -1 ------------------------------------------


Loops.a1.in.core.cores <- cores_all[indices.c.a1.v2@to,] # 14 x 5
Loops.a1.in.core.cores.uniq <- Loops.a1.in.core.cores %>% # 10 x 6
    unique() %>% # 35 x 4
    tibble::add_column(loops_a1_per_core = table(Loops.a1.in.core.cores$myid_core)[(.)$myid_core])

Loops.a1.in.core.loops <- loops[indices.c.a1.v2@from,] %>% # 14 x 4
    tibble::add_column(myid_core = Loops.a1.in.core.cores$myid_core) %>%
    tibble::add_column(coreset = Loops.a1.in.core.cores$coreset)

core_ids <- list()
core_sets <- list()
as_cores_per_loop <- numeric()
gfp_cores_per_loop <- numeric()
for(i in 1:nrow(Loops.a1.in.core.loops)){
    core_ids[i] <- list(Loops.a1.in.core.loops[Loops.a1.in.core.loops$myid_loop == Loops.a1.in.core.loops$myid_loop[i],'myid_core'])
    core_sets[i] <- list(Loops.a1.in.core.loops[Loops.a1.in.core.loops$myid_loop == Loops.a1.in.core.loops$myid_loop[i],'coreset'])
    as_cores_count <- 0
    gfp_cores_count <- 0
    for(j in 1:length(core_sets[[i]])){
        if(core_sets[[i]][j] == 'as'){as_cores_count = as_cores_count + 1}
        else if(core_sets[[i]][j] == 'gfp'){gfp_cores_count = gfp_cores_count + 1}
    }
    as_cores_per_loop[i] <- as_cores_count
    gfp_cores_per_loop[i] <- gfp_cores_count
}
Loops.a1.in.core.loops.uniq <-  Loops.a1.in.core.loops %>% # 9 x 5
    tibble::add_column(cores_per_loop_a1_v2 = table(Loops.a1.in.core.loops$myid_loop)[(.)$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(cores_per_loop_ids_v2 = core_ids) %>%
    tibble::add_column(cores_per_loop_sets_v2 = core_sets) %>%
    tibble::add_column(as_cores_per_loop_v2 = as_cores_per_loop) %>%
    tibble::add_column(gfp_cores_per_loop_v2 = gfp_cores_per_loop) %>%
    dplyr::select(., -myid_core) %>%
    dplyr::select(., -coreset) %>%
    unique()
dim(Loops.a1.in.core.loops.uniq) # 9x 10



# multiple loops a1 per core:
# The same core can contain multiple loop anchor 1 (10000), the startI of the loops
#    can be the same with different startJ, or the startI can be different and just
#    be within the core (start and end)
# For example,
# # the same core:
# 4      chr17  43765228  43835963   c.as1.4   70735
# 4.1    chr17  43765228  43835963   c.as1.4   70735
# 4.2    chr17  43765228  43835963   c.as1.4   70735
# # different multiple loops, startI can be the same or different
# 7511   chr17  43770000  44740000   Ldg7511
# 7512   chr17  43770000  44750000   Ldg7512
# 7513   chr17  43810000  44740000   Ldg7513

# multiple cores per loop a1:
# the same loop anchor 1 (10000) can be in multiple cores (that have different widths)
# # For example, 164940000-164949999 in [3368    chr5 164940000 165260000   Ldg3368] is in 
#          chr     start       end myid_core   width
# 3100    chr5 164878361 164978443   c.as2.3  100082
# 1112    chr5 164884536 164952133 c.gfp1.11   67597
# 914     chr5 164884689 164978429  c.gfp2.9   93740


Loops.a2.in.core.cores <- cores_all[indices.c.a2.v2@to,] # 24 x 5
Loops.a2.in.core.cores.uniq <- Loops.a2.in.core.cores %>% # 15 x 6
    unique() %>% 
    tibble::add_column(loops_a2_per_core = table(Loops.a2.in.core.cores$myid_core)[(.)$myid_core])

Loops.a2.in.core.loops <- loops[indices.c.a2.v2@from,] %>% # 24 x 4
    tibble::add_column(myid_core = Loops.a2.in.core.cores$myid_core) %>%
    tibble::add_column(coreset = Loops.a2.in.core.cores$coreset)

core_ids <- list()
core_sets <- list()
as_cores_per_loop <- numeric()
gfp_cores_per_loop <- numeric()
for(i in 1:nrow(Loops.a2.in.core.loops)){
    core_ids[i] <- list(Loops.a2.in.core.loops[Loops.a2.in.core.loops$myid_loop == Loops.a2.in.core.loops$myid_loop[i],'myid_core'])
    core_sets[i] <- list(Loops.a2.in.core.loops[Loops.a2.in.core.loops$myid_loop == Loops.a2.in.core.loops$myid_loop[i],'coreset'])
    as_cores_count <- 0
    gfp_cores_count <- 0
    for(j in 1:length(core_sets[[i]])){
        if(core_sets[[i]][j] == 'as'){as_cores_count = as_cores_count + 1}
        else if(core_sets[[i]][j] == 'gfp'){gfp_cores_count = gfp_cores_count + 1}
    }
    as_cores_per_loop[i] <- as_cores_count
    gfp_cores_per_loop[i] <- gfp_cores_count
}
Loops.a2.in.core.loops.uniq <-  Loops.a2.in.core.loops %>% # 12 x 5
    tibble::add_column(cores_per_loop_a2_v2 = table(Loops.a2.in.core.loops$myid_loop)[(.)$myid_loop] %>% as.numeric()) %>%
    tibble::add_column(cores_per_loop_ids_v2 = core_ids) %>%
    tibble::add_column(cores_per_loop_sets_v2 = core_sets) %>%
    tibble::add_column(as_cores_per_loop_v2 = as_cores_per_loop) %>%
    tibble::add_column(gfp_cores_per_loop_v2 = gfp_cores_per_loop) %>%
    dplyr::select(., -myid_core) %>%
    dplyr::select(., -coreset) %>%
    unique()
dim(Loops.a2.in.core.loops.uniq) # 12 x 10


# Intersect - Promoters and COREs

# prom.in.Loops.a1.loops.uniq # 279 x 5
# prom.in.Loops.a2.loops.uniq # 273 x 5

# core.in.Loops.a1.loops.uniq # 40 x 5 # Ldg.a1.in.core.loops.uniq # 9 x 5
# core.in.Loops.a2.loops.uniq # 50 x 5 # Loops.a2.in.core.loops.uniq # 12 x 5

#core.Loops.a1.loops.uniq.id <- union(core.in.Loops.a1.loops.uniq$myid_loop, Ldg.a1.in.core.loops.uniq$myid_loop) # 46
#core.Loops.a2.loops.uniq.id <- union(core.in.Loops.a2.loops.uniq$myid_loop, Loops.a2.in.core.loops.uniq$myid_loop) # 61
core.Loops.a1.loops.uniq.join <- dplyr::full_join(core.in.Loops.a1.loops.uniq, Loops.a1.in.core.loops.uniq)
core.Loops.a2.loops.uniq.join <- dplyr::full_join(core.in.Loops.a2.loops.uniq, Loops.a2.in.core.loops.uniq)

pa1.ca2.loops.id <- intersect(prom.in.Loops.a1.loops.uniq$myid_loop, core.Loops.a2.loops.uniq.join$myid_loop) # 5
pa1.ca2.loops.temp1 <- prom.in.Loops.a1.loops.uniq[prom.in.Loops.a1.loops.uniq$myid_loop %in% pa1.ca2.loops.id,]
pa1.ca2.loops.temp2 <- core.Loops.a2.loops.uniq.join[core.Loops.a2.loops.uniq.join$myid_loop %in% pa1.ca2.loops.id,]
pa1.ca2.loops <- dplyr::left_join(pa1.ca2.loops.temp1, pa1.ca2.loops.temp2) # 5 x 7

ca1.pa2.loops.id <- intersect(core.Loops.a1.loops.uniq.join$myid_loop, prom.in.Loops.a2.loops.uniq$myid_loop) # 3
ca1.pa2.loops.temp1 <- core.Loops.a1.loops.uniq.join[core.Loops.a1.loops.uniq.join$myid_loop %in% ca1.pa2.loops.id,] # 3 x 6
ca1.pa2.loops.temp2 <- prom.in.Loops.a2.loops.uniq[prom.in.Loops.a2.loops.uniq$myid_loop %in% ca1.pa2.loops.id,] # 3 x 5
ca1.pa2.loops <- dplyr::left_join(ca1.pa2.loops.temp1, ca1.pa2.loops.temp2) # 3 x 7


# Ldg: 
# pa1.ca2.loops # 5 x 7
# ca1.pa2.loops # 3 x 7

# Ldl:
# pa1.ca2.loops # 3 x 7
# ca1.pa2.loops # 13 x 7

# Ldung:
# pa1.ca2.loops # 36 x 7
# ca1.pa2.loops # 32 x 7

# Lddnl:
# pa1.ca2.loops # 38 x 7
# ca1.pa2.loops # 48 x 7

# Ls:
# pa1.ca2.loops # 8 x 7
# ca1.pa2.loops # 15 x 7

# loops_all:
# pa1.ca2.loops # 90 x 18
# ca1.pa2.loops # 111 x 18
dim(pa1.ca2.loops)
dim(ca1.pa2.loops)


# Final dataframes ----------------------------------------------------------------------------

a1prom.a2core.loops <- pa1.ca2.loops %>% replace(is.na(.), 0) %>% # 90 x 17
    tibble::add_column(as_cores_per_loop_TOTAL = (.)$as_cores_per_loop + (.)$as_cores_per_loop_v2) %>%
    tibble::add_column(gfp_cores_per_loop_TOTAL = (.)$gfp_cores_per_loop + (.)$gfp_cores_per_loop_v2)

a1prom.a2core.loops.Ls <- a1prom.a2core.loops[a1prom.a2core.loops$loopset == 'Ls',] # 8 x 19
a1prom.a2core.loops.Ldg <- a1prom.a2core.loops[a1prom.a2core.loops$loopset == 'Ldg',] # 5 x 19
a1prom.a2core.loops.Ldung <- a1prom.a2core.loops[a1prom.a2core.loops$loopset == 'Ldung',] # 36 x 19
a1prom.a2core.loops.Ldl <- a1prom.a2core.loops[a1prom.a2core.loops$loopset == 'Ldl',] # 3 x 19
a1prom.a2core.loops.Lddnl <- a1prom.a2core.loops[a1prom.a2core.loops$loopset == 'Lddnl',] # 38 x 19 

a1prom.a2core.loops.dims <- data.frame(
    loop.set = c('Total', 'Sig.not.diff', 'Diff.gain', 'Diff.up.not.gain', 'Diff.lost', 'Diff.down.not.lost'),
    number.of.loops = c(nrow(a1prom.a2core.loops),nrow(a1prom.a2core.loops.Ls),
                        nrow(a1prom.a2core.loops.Ldg), nrow(a1prom.a2core.loops.Ldung),
                        nrow(a1prom.a2core.loops.Ldl), nrow(a1prom.a2core.loops.Lddnl)),
    num.loops.assoc.with.AS.cores = c(
        sum(a1prom.a2core.loops$as_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ls$as_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldg$as_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldung$as_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldl$as_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Lddnl$as_cores_per_loop_TOTAL > 0)),
    num.loops.assoc.with.GFP.cores = c(
        sum(a1prom.a2core.loops$gfp_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ls$gfp_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldg$gfp_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldung$gfp_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Ldl$gfp_cores_per_loop_TOTAL > 0),
        sum(a1prom.a2core.loops.Lddnl$gfp_cores_per_loop_TOTAL > 0))
)

a1core.a2prom.loops <- ca1.pa2.loops %>% replace(is.na(.), 0) %>% # 111 x 17
    tibble::add_column(as_cores_per_loop_TOTAL = (.)$as_cores_per_loop + (.)$as_cores_per_loop_v2) %>%
    tibble::add_column(gfp_cores_per_loop_TOTAL = (.)$gfp_cores_per_loop + (.)$gfp_cores_per_loop_v2)

a1core.a2prom.loops.Ls <- a1core.a2prom.loops[a1core.a2prom.loops$loopset == 'Ls',] # 15 x 19
a1core.a2prom.loops.Ldg <- a1core.a2prom.loops[a1core.a2prom.loops$loopset == 'Ldg',] # 3 x 19
a1core.a2prom.loops.Ldung <- a1core.a2prom.loops[a1core.a2prom.loops$loopset == 'Ldung',] # 32 x 19
a1core.a2prom.loops.Ldl <- a1core.a2prom.loops[a1core.a2prom.loops$loopset == 'Ldl',] # 13 x 19
a1core.a2prom.loops.Lddnl <- a1core.a2prom.loops[a1core.a2prom.loops$loopset == 'Lddnl',] # 48 x 19

a1core.a2prom.loops.dims <- data.frame(
    loop.set = c('Total', 'Sig.not.diff', 'Diff.gain', 'Diff.up.not.gain', 'Diff.lost', 'Diff.down.not.lost'),
    number.of.loops = c(nrow(a1core.a2prom.loops), nrow(a1core.a2prom.loops.Ls),
                        nrow(a1core.a2prom.loops.Ldg), nrow(a1core.a2prom.loops.Ldung),
                        nrow(a1core.a2prom.loops.Ldl), nrow(a1core.a2prom.loops.Lddnl)),
    num.loops.assoc.with.AS.cores = c(
        sum(a1core.a2prom.loops$as_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ls$as_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldg$as_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldung$as_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldl$as_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Lddnl$as_cores_per_loop_TOTAL > 0)),
    num.loops.assoc.with.GFP.cores = c(
        sum(a1core.a2prom.loops$gfp_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ls$gfp_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldg$gfp_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldung$gfp_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Ldl$gfp_cores_per_loop_TOTAL > 0),
        sum(a1core.a2prom.loops.Lddnl$gfp_cores_per_loop_TOTAL > 0))
)

# Check counts
# table(a1prom.a2core.loops$loopset)
# Lddnl   Ldg   Ldl Ldung    Ls 
# 89     9    15    92    27 
# table(a1core.a2prom.loops$loopset)
# Lddnl   Ldg   Ldl Ldung    Ls 
# 105     8    29    85    31 

saveRDS(a1prom.a2core.loops, "./output/R_objects/a1prom.a2core.loops.any1bp.overlap.rds")
saveRDS(a1core.a2prom.loops, "./output/R_objects/a1core.a2prom.loops.any1bp.overlap.rds")


# Write to file: promoter-CORE loops ----------------------------------------------------------

write.xlsx2(a1prom.a2core.loops, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_all.xlsx"))
write.xlsx2(a1prom.a2core.loops.Ls, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_Ls.xlsx"))
write.xlsx2(a1prom.a2core.loops.Ldg, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_Ldg.xlsx"))
write.xlsx2(a1prom.a2core.loops.Ldung, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_Ldung.xlsx"))
write.xlsx2(a1prom.a2core.loops.Ldl, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_Ldl.xlsx"))
write.xlsx2(a1prom.a2core.loops.Lddnl, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_Lddnl.xlsx"))
write.xlsx2(a1prom.a2core.loops.dims, file = glue(outdir_pcloops, '/a1prom_a2core/',"a1prom_a2core_loops_count.xlsx"))

write.xlsx2(a1core.a2prom.loops, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_all.xlsx"))
write.xlsx2(a1core.a2prom.loops.Ls, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_Ls.xlsx"))
write.xlsx2(a1core.a2prom.loops.Ldg, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_Ldg.xlsx"))
write.xlsx2(a1core.a2prom.loops.Ldung, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_Ldung.xlsx"))
write.xlsx2(a1core.a2prom.loops.Ldl, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_Ldl.xlsx"))
write.xlsx2(a1core.a2prom.loops.Lddnl, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_Lddnl.xlsx"))
write.xlsx2(a1core.a2prom.loops.dims, file = glue(outdir_pcloops, '/a1core_a2prom/',"a1core_a2prom_loops_count.xlsx"))


# Genes dataframe: promoter-CORE loops data transposed ----------------------------------------
dim(a1prom.a2core.loops) # 232 x 20
dim(a1core.a2prom.loops) # 258 x 20

genes.a1prom.a2core.temp <- a1prom.a2core.loops %>% 
    dplyr::select('myid_loop','proms_per_loop_ids','gene_symbol','loopset','cores_per_loop_ids','cores_per_loop_ids_v2')
genes.a1core.a2prom.temp <- a1core.a2prom.loops %>% 
    dplyr::select('myid_loop','proms_per_loop_ids','gene_symbol','loopset','cores_per_loop_ids','cores_per_loop_ids_v2')

# all promoters in each row but may be mapped to same gene
genes.a1prom.a2core <- tibble(
    gene_symbol = character(),prom_ids = character(),core_ids = list(),myid_loop = character(), loopset = character())
for(i in 1:nrow(genes.a1prom.a2core.temp)){
    proms <- genes.a1prom.a2core.temp$proms_per_loop_ids[[i]]
    cores <- list(c(genes.a1prom.a2core.temp$cores_per_loop_ids[[i]], genes.a1prom.a2core.temp$cores_per_loop_ids_v2[[i]]))
    names(cores) <- genes.a1prom.a2core.temp$myid_loop[[i]]
    symb <- genes.a1prom.a2core.temp$gene_symbol[[i]]
    for(p in 1:length(proms)){
        genes.a1prom.a2core <- genes.a1prom.a2core %>% 
            tibble::add_row(gene_symbol = symb[p],
                            prom_ids = proms[p], 
                            core_ids = cores,
                            myid_loop = genes.a1prom.a2core.temp$myid_loop[[i]],
                            loopset = genes.a1prom.a2core.temp$loopset[[i]])
    }
}
genes.a1prom.a2core <- data.frame(genes.a1prom.a2core) 

genes.a1prom.a2core.uniq <- tibble(gene_symbol = character(), 
                                   num_uniq_proms = numeric(), num_uniq_loops = numeric(), 
                                   num_uniq_cores = numeric(), num_uniq_large_cores = numeric(),
                                   num_as1_cores = numeric(), num_as2_cores = numeric(), 
                                   num_gfp1_cores = numeric(), num_gfp2_cores = numeric(),
                                   both_as_cores = character(), both_gfp_cores = character(), as_gfp_cores = character(),
                                   uniq_loopsets = list(),
                                   uniq_prom_ids = list(), uniq_loop_ids = list(), uniq_core_ids = list(),
                                   num_proms = numeric(), num_loops = numeric(),
                                   prom_ids = list(), loop_ids = list(), core_ids_per_loop = list())
symbols.a1prom.a2core.uniq <- unique(genes.a1prom.a2core$gene_symbol) # 78
for(i in 1:length(symbols.a1prom.a2core.uniq)){
    genes.subset <- genes.a1prom.a2core %>% dplyr::filter(gene_symbol == symbols.a1prom.a2core.uniq[i])
    proms.add <- list(genes.subset$prom_ids)
    loops.add <- list(genes.subset$myid_loop)
    num.proms.add <- length(proms.add[[1]])
    num.loops.add <- length(loops.add[[1]])
    
    num.uniq.proms.add <- length(unique(proms.add[[1]]))
    num.uniq.loops.add <- length(unique(loops.add[[1]]))
    uniq.proms.add <- list(unique(proms.add[[1]]))
    uniq.loops.add <- list(unique(loops.add[[1]]))
    uniq.loopsets.add <- genes.subset$loopset %>% unique() %>% list()
    
    cores.add <- genes.subset$core_ids[uniq.loops.add[[1]]] %>% list()
    uniq.cores.add <- genes.subset$core_ids[uniq.loops.add[[1]]] %>% unlist() %>% unique() %>% list()
    num.uniq.cores.add <- uniq.cores.add[[1]] %>% length()      
    num.uniq.large.cores.add <- regmatches(uniq.cores.add[[1]], gregexpr('C', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.as1.cores <- regmatches(uniq.cores.add[[1]], gregexpr('as1', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.as2.cores <- regmatches(uniq.cores.add[[1]], gregexpr('as2', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.gfp1.cores <- regmatches(uniq.cores.add[[1]], gregexpr('gfp1', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.gfp2.cores <- regmatches(uniq.cores.add[[1]], gregexpr('gfp2', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    
    both.as.cores <- 'no'
    both.gfp.cores <- 'no'
    as.gfp.cores <- 'no'
    if(num.as1.cores > 0 & num.as2.cores > 0){both.as.cores = 'yes'}
    if(num.gfp1.cores > 0 & num.gfp2.cores > 0){both.gfp.cores = 'yes'}
    if((num.as1.cores > 0 | num.as2.cores > 0) & (num.gfp1.cores > 0 | num.gfp2.cores > 0)){as.gfp.cores = 'yes'}
    
    genes.a1prom.a2core.uniq <- genes.a1prom.a2core.uniq %>%
        tibble::add_row(gene_symbol = genes.subset$gene_symbol[1],
                        num_uniq_proms = num.uniq.proms.add,
                        num_uniq_loops = num.uniq.loops.add,
                        num_uniq_cores = num.uniq.cores.add,
                        num_uniq_large_cores = num.uniq.large.cores.add,
                        num_as1_cores = num.as1.cores,
                        num_as2_cores = num.as2.cores,
                        num_gfp1_cores = num.gfp1.cores,
                        num_gfp2_cores = num.gfp2.cores,
                        both_as_cores = both.as.cores,
                        both_gfp_cores = both.gfp.cores,
                        as_gfp_cores = as.gfp.cores,
                        uniq_prom_ids = uniq.proms.add,
                        uniq_loop_ids = uniq.loops.add,
                        uniq_loopsets = uniq.loopsets.add,
                        uniq_core_ids = uniq.cores.add,
                        num_proms = num.proms.add,
                        num_loops = num.loops.add,
                        prom_ids = proms.add,
                        loop_ids = loops.add,
                        core_ids_per_loop = cores.add)
}
# genes.a1prom.a2core.uniq <- data.frame(genes.a1prom.a2core.uniq)
genes.a1prom.a2core.uniq.fin <- genes.a1prom.a2core.uniq %>% dplyr::select(., -c(num_proms, num_loops, prom_ids, loop_ids))
# dim(genes.a1prom.a2core.uniq.fin) # 179 x 17


# a1core.a2prom
genes.a1core.a2prom <- tibble(
    gene_symbol = character(),prom_ids = character(),core_ids = list(),myid_loop = character(), loopset = character())
for(i in 1:nrow(genes.a1core.a2prom.temp)){
    proms <- genes.a1core.a2prom.temp$proms_per_loop_ids[[i]]
    cores <- list(c(genes.a1core.a2prom.temp$cores_per_loop_ids[[i]], genes.a1core.a2prom.temp$cores_per_loop_ids_v2[[i]]))
    names(cores) <- genes.a1core.a2prom.temp$myid_loop[[i]]
    symb <- genes.a1core.a2prom.temp$gene_symbol[[i]]
    for(p in 1:length(proms)){
        genes.a1core.a2prom <- genes.a1core.a2prom %>% 
            tibble::add_row(gene_symbol = symb[p],
                            prom_ids = proms[p], 
                            core_ids = cores,
                            myid_loop = genes.a1core.a2prom.temp$myid_loop[[i]],
                            loopset = genes.a1core.a2prom.temp$loopset[[i]])
    }
}
genes.a1core.a2prom <- data.frame(genes.a1core.a2prom) # 106 x 5

genes.a1core.a2prom.uniq <- tibble(gene_symbol = character(), 
                                   num_uniq_proms = numeric(), num_uniq_loops = numeric(), 
                                   num_uniq_cores = numeric(), num_uniq_large_cores = numeric(),
                                   num_as1_cores = numeric(), num_as2_cores = numeric(), 
                                   num_gfp1_cores = numeric(), num_gfp2_cores = numeric(),
                                   both_as_cores = character(), both_gfp_cores = character(), as_gfp_cores = character(),
                                   uniq_loopsets = list(),
                                   uniq_prom_ids = list(), uniq_loop_ids = list(), uniq_core_ids = list(),
                                   num_proms = numeric(), num_loops = numeric(),
                                   prom_ids = list(), loop_ids = list(), core_ids_per_loop = list())
symbols.a1core.a2prom.uniq <- unique(genes.a1core.a2prom$gene_symbol) # 78
for(i in 1:length(symbols.a1core.a2prom.uniq)){
    genes.subset <- genes.a1core.a2prom %>% dplyr::filter(gene_symbol == symbols.a1core.a2prom.uniq[i])
    proms.add <- list(genes.subset$prom_ids)
    loops.add <- list(genes.subset$myid_loop)
    num.proms.add <- length(proms.add[[1]])
    num.loops.add <- length(loops.add[[1]])
    
    num.uniq.proms.add <- length(unique(proms.add[[1]]))
    num.uniq.loops.add <- length(unique(loops.add[[1]]))
    uniq.proms.add <- list(unique(proms.add[[1]]))
    uniq.loops.add <- list(unique(loops.add[[1]]))
    uniq.loopsets.add <- genes.subset$loopset %>% unique() %>% list()
    
    cores.add <- genes.subset$core_ids[uniq.loops.add[[1]]] %>% list()
    uniq.cores.add <- genes.subset$core_ids[uniq.loops.add[[1]]] %>% unlist() %>% unique() %>% list()
    num.uniq.cores.add <- uniq.cores.add[[1]] %>% length()      
    num.uniq.large.cores.add <- regmatches(uniq.cores.add[[1]], gregexpr('C', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.as1.cores <- regmatches(uniq.cores.add[[1]], gregexpr('as1', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.as2.cores <- regmatches(uniq.cores.add[[1]], gregexpr('as2', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.gfp1.cores <- regmatches(uniq.cores.add[[1]], gregexpr('gfp1', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    num.gfp2.cores <- regmatches(uniq.cores.add[[1]], gregexpr('gfp2', uniq.cores.add[[1]])) %>% lengths() %>% sum()
    
    both.as.cores <- 'no'
    both.gfp.cores <- 'no'
    as.gfp.cores <- 'no'
    if(num.as1.cores > 0 & num.as2.cores > 0){both.as.cores = 'yes'}
    if(num.gfp1.cores > 0 & num.gfp2.cores > 0){both.gfp.cores = 'yes'}
    if((num.as1.cores > 0 | num.as2.cores > 0) & (num.gfp1.cores > 0 | num.gfp2.cores > 0)){as.gfp.cores = 'yes'}
    
    genes.a1core.a2prom.uniq <- genes.a1core.a2prom.uniq %>%
        tibble::add_row(gene_symbol = genes.subset$gene_symbol[1],
                        num_uniq_proms = num.uniq.proms.add,
                        num_uniq_loops = num.uniq.loops.add,
                        num_uniq_cores = num.uniq.cores.add,
                        num_uniq_large_cores = num.uniq.large.cores.add,
                        num_as1_cores = num.as1.cores,
                        num_as2_cores = num.as2.cores,
                        num_gfp1_cores = num.gfp1.cores,
                        num_gfp2_cores = num.gfp2.cores,
                        both_as_cores = both.as.cores,
                        both_gfp_cores = both.gfp.cores,
                        as_gfp_cores = as.gfp.cores,
                        uniq_prom_ids = uniq.proms.add,
                        uniq_loop_ids = uniq.loops.add,
                        uniq_loopsets = uniq.loopsets.add,
                        uniq_core_ids = uniq.cores.add,
                        num_proms = num.proms.add,
                        num_loops = num.loops.add,
                        prom_ids = proms.add,
                        loop_ids = loops.add,
                        core_ids_per_loop = cores.add)
}
# genes.a1core.a2prom.uniq <- data.frame(genes.a1core.a2prom.uniq)
# dim(genes.a1core.a2prom.uniq) # 86 x 21
genes.a1core.a2prom.uniq.fin <- genes.a1core.a2prom.uniq %>% dplyr::select(., -c(num_proms, num_loops, prom_ids, loop_ids))
# dim(genes.a1core.a2prom.uniq.fin)
# 185 x 17

# Save transposed to genes in rows
write.xlsx2(genes.a1prom.a2core.uniq.fin, file = glue(outdir_pcloops, '/a1prom_a2core/',"genes_a1prom_a2core_uniq.xlsx"))
write.xlsx2(genes.a1core.a2prom.uniq.fin, file = glue(outdir_pcloops, '/a1core_a2prom/',"genes_a1core_a2prom_uniq.xlsx"))

saveRDS(genes.a1prom.a2core.uniq.fin, file = glue("./output/R_objects/", "genes.a1prom.a2core.uniq.fin.any1bp.rds"))
saveRDS(genes.a1core.a2prom.uniq.fin, file = glue("./output/R_objects/", "genes.a1core.a2prom.uniq.fin.any1bp.rds"))
write.xlsx2(genes.a1prom.a2core.uniq.fin, file = glue("./output/R_objects_files/", "genes_a1prom_a2core_uniq_any1bp.xlsx"))
write.xlsx2(genes.a1core.a2prom.uniq.fin, file = glue("./output/R_objects_files/", "genes_a1core_a2prom_uniq_any1bp.xlsx"))

# write.xlsx2(cores_all, file = glue(outdir_pcloops, "/cores_all_ids.xlsx"))

