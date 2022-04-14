library(plotgardener)
library(glue)
library(magrittr)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(GenomicRanges)
library(IRanges)

outdir <- "./output/FIGURES"

# R color palettes: https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/
# cheatsheet: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf

# Region to plot ------------------------------------------------------------------------------

# myparams <- pgParams(chrom = "chr5",
#                         chromstart = 164700000, #164700000, #163800000, #163800000, #164700000, # 163800000 # for all chr5 tads
#                         chromend = 165950000, #165950000, #166700000, # 166700000, #166000000,
#                         assembly = "rn6")

# myparams <- pgParams(chrom = "chr5",
#                      chromstart = 164500000, #164700000, #163800000, #163800000, #164700000, # 163800000 # for all chr5 tads
#                      chromend = 165970000, #165950000, #166700000, # 166700000, #166000000,
#                      assembly = "rn6")
# HiC -----------------------------------------------------------------------------------------

plot_hic <- function(myparams=NULL, cond=NULL, samp=NULL, mypal=NULL, z=NULL, x_add=NULL, y_add=NULL){
    hic_x = 1 + x_add; hic_y = -1 + y_add; hic_width = 7; hic_height = sqrt(7^2/2)
    
    hicfile <- NULL
    if(is.null(samp)){
        if (cond == "as"){
            hicfile <- "./output/CONVERT/converted_hic/as_10000_merged_iced_all_intra.hic"
        } else if (cond == "gfp") {
            hicfile <- "./output/CONVERT/converted_hic/gfp_10000_merged_iced_all_intra.hic"
        } 
    } else {
        hicfile <- glue("./output/CONVERT/converted_hic/{samp}_10000_iced_all_intra.hic")
    }
    
    
    if (!is.null(z)) {
        hic_plot <- plotHicTriangle(data = hicfile,
                                    params = myparams,
                                    norm = "NONE",
                                    zrange = z,
                                    x = hic_x, 
                                    y = hic_y,
                                    width = hic_width,
                                    height = hic_height,
                                    colorTrans = "log",
                                    palette = mypal,
                                    assembly = "rn6")
    } else if (is.null(z)){
        hic_plot <- plotHicTriangle(data = hicfile,
                                    params = myparams,
                                    norm = "NONE",
                                    x = hic_x, 
                                    y = hic_y,
                                    width = hic_width,
                                    height = hic_height,
                                    colorTrans = "log",
                                    palette = mypal,
                                    assembly = "rn6")
    }
    return(list(plot = hic_plot, x = hic_x, y = hic_y, width = hic_width, height = hic_height))
}


# Genes ---------------------------------------------------------------------------------------

plot_genes <- function(myplot=NULL, myparams=NULL, x_add=NULL, y_add=NULL){
    genes_height = 1.2
    y_new = myplot$y + myplot$height + y_add
    
    genes_plot <- plotGenes(
        params = myparams, 
        strandLabels = TRUE,
        fontsize = 20, 
        x = myplot$x, 
        y = y_new,
        width = myplot$width,
        height = genes_height,
        geneHighlights = data.frame(
            gene = c("Mtor","Egr2","Jun","Junb","Jund","Fosl1","Fos","Fosl2","Fosb"),
            color = "#290284"),
        geneBackground = c("black"),
        fill = c("black", "black"), 
        fontcolor = c("black", "black"))
    return(list(plot = genes_plot, height = genes_height, y = y_new))
}


plot_glabel <- function(myplot = NULL, y_add = 0.1, fsize = 20){
    y_new = myplot$y + myplot$height + y_add
    annoGenomeLabel(plot = myplot$plot, 
                    x = myplot$x, 
                    y = y_new, 
                    scale = "Mb",
                    fontsize = fsize,
                    just = c("left", "top"))
    
    return(list(y = y_new))
}



plot_ontads <- function(cond=NULL, myplot=NULL, res=NULL, myparams=NULL, level1_flag=FALSE){
    tads <- NULL
    if(cond == "as"){ 
        tads <- tads_as_all[[myparams$chrom]] } 
    else if (cond == "gfp"){ 
        tads <- tads_gfp_all[[myparams$chrom]] }
    
    if(level1_flag == TRUE){
        tads <- tads %>% dplyr::filter(TADlevel == 1)
    }
    
    gr_tads <- GRanges(seqnames = myparams$chrom,
                       ranges = IRanges(start = tads$startpos * res,
                                        end = tads$endpos * res))
    tads_to_plot <- gr_tads
    
    annoDomains(plot = myplot$plot,
                data = tads_to_plot,
                linecolor = "black")
    
    # tads_level2 <- tads %>% dplyr::filter(TADlevel == 2)
    # gr_tads_level2 <- GRanges(seqnames =  myparams$chrom,
    #                           ranges = IRanges(start = tads_level2$startpos * res,
    #                                            end = tads_level2$endpos * res))
    # annoDomains(plot = myplot$plot,
    #             data = gr_tads_level2,
    #             linecolor = "pink")
    # 
    # tads_level3 <- tads %>% dplyr::filter(TADlevel == 3)
    # gr_tads_level3 <- GRanges(seqnames =  myparams$chrom,
    #                           ranges = IRanges(start = tads_level3$startpos * res,
    #                                            end = tads_level3$endpos * res))
    # annoDomains(plot = myplot$plot,
    #             data = gr_tads_level3,
    #             linecolor = "blue")
    # 
    # tads_level4 <- tads %>% dplyr::filter(TADlevel == 4)
    # gr_tads_level4 <- GRanges(seqnames =  myparams$chrom,
    #                           ranges = IRanges(start = tads_level4$startpos * res,
    #                                            end = tads_level4$endpos * res))
    # annoDomains(plot = myplot$plot,
    #             data = gr_tads_level4,
    #             linecolor = "green")
    # 
    # 
    # tads_level5 <- tads %>% dplyr::filter(TADlevel == 5)
    # gr_tads_level5 <- GRanges(seqnames =  myparams$chrom,
    #                           ranges = IRanges(start = tads_level5$startpos * res,
    #                                            end = tads_level5$endpos * res))
    # annoDomains(plot = myplot$plot,
    #             data = gr_tads_level5,
    #             linecolor = "red")
    # 
    # tads_level6 <- tads %>% dplyr::filter(TADlevel == 6)
    # gr_tads_level6 <- GRanges(seqnames =  myparams$chrom,
    #                           ranges = IRanges(start = tads_level6$startpos * res,
    #                                            end = tads_level6$endpos * res))
    # annoDomains(plot = myplot$plot,
    #             data = gr_tads_level6,
    #             linecolor = "yellow")
}

plot_cores <- function(myplot = NULL, y_add = NULL, myparams=NULL){
    x_new <- myplot$x
    y_new <- myplot$y + myplot$height + y_add
    width_new <- myplot$width
    height_new <- 1
    
    cores_to_plot <- data.frame(
        chrom1 = cores_all$chr,
        start1 = cores_all$start,
        end1 = cores_all$end,
        chrom2 = cores_all$chr,
        start2 = cores_all$start,
        end2 = cores_all$end,
        coreset = cores_all$coreset) %>%
        dplyr::mutate(coreset_num = dplyr::case_when(
            endsWith(coreset, "as") ~ 1,
            endsWith(coreset, "gfp") ~ 0
        )) 
    mypairs <- plotPairs(data = cores_to_plot,
                         params = myparams,
                         #chrom = chrm,
                         #chromstart = chrm_s,
                         #chromend = chrm_e,
                         assembly = "rn6",
                         x = x_new, 
                         y = y_new,
                         width = width_new,
                         height = height_new,
                         fill = colorby(column = "coreset_num", 
                                        palette = colorRampPalette(c("#0C47D8","#DC2222"), bias=1))
                         #palette = colorRampPalette(c("blue","red"), bias=1))
                         #palette = colorRampPalette(c("#0C40BF","#F15454"), bias=1))
    )
    return(list(x = x_new, y = y_new, width = width_new, height = height_new))
}


plot_loops <- function(myplot = NULL, y_add = NULL, myparams=NULL, mtor_flag=FALSE){
    
    prep_loops_data <- function(loops_data, y_add, mtor_flag=FALSE){
        loops_to_plot <- loops_data %>% 
            dplyr::filter(chr == myparams$chrom) %>% 
            dplyr::filter(startI >= myparams$chromstart) %>%
            dplyr::filter(startJ <= myparams$chromend)
        
        if(mtor_flag == TRUE){
            loops_to_plot <- loops_to_plot %>% dplyr::filter(grepl("Mtor", gene_symbol))
        }
        
        loops_to_plot <- loops_to_plot %>%    
            dplyr::select(chr, startI, startJ, loopset) %>%
            tibble::add_column(chrom1 = (.)$chr) %>%
            tibble::add_column(start1 = (.)$startI) %>%
            tibble::add_column(end1 = (.)$start1 + 10000) %>% 
            tibble::add_column(chrom2 = (.)$chr) %>%
            tibble::add_column(start2 = (.)$startJ) %>%
            tibble::add_column(end2 = (.)$start2 + 10000) %>%
            tibble::add_column(length = floor((.)$start2 - (.)$start1)) %>%
            dplyr::select(-chr, -startI, -startJ) %>%
            dplyr::select(chrom1, start1, end1, chrom2, start2, end2, length)
        return(loops_to_plot)
    }
    
    loops_to_plot <- dplyr::bind_rows(
        prep_loops_data(loops_data = loops_pc, y_add=y_add),
        prep_loops_data(loops_data = loops_cp, y_add=y_add))
    
    if(mtor_flag == TRUE){
        loops_to_plot <- prep_loops_data(loops_data = loops_cp, y_add=y_add, mtor_flag=TRUE)
    }
    
    x_new = myplot$x 
    
    y_new = myplot$y + myplot$height + y_add
    width_new = myplot$width
    height_new = 0.5
    
    myarcs <- plotPairsArches(
        data = loops_to_plot,
        params = myparams,
        #fill = colorby("loopset_num", palette = colorRampPalette(c("#9CABDE","#6D83B5","#3A58A0"))),
        archHeight = loops_to_plot$length, 
        alpha = 0.4, 
        curvature = 10,
        flip = TRUE,
        x = x_new, 
        y = y_new,
        width = width_new,
        height = height_new,
        default.units = "inches",
        fill = "black",
        linecolor = "black"
    ) 
    return(list(x = x_new, y = y_new, width = width_new, height = height_new))
}

plot_tad_signal <- function(myplot = output_hic, cond=NULL, y_add = NULL){
    tads <- NULL
    ifelse(cond == "as", tads <- tads_as, tads <- tads_gfp)
    
    tad_signal <- tads$chr5$binSignal %>%
        dplyr::mutate(chrom = chr) %>%
        dplyr::mutate(start = from.coord) %>%
        dplyr::mutate(end = to.coord) %>%
        dplyr::mutate(end = to.coord-1) %>%
        dplyr::mutate(score = mean.cf) %>%
        dplyr::select(chrom, start, end, score, pvalue) 
    
    # astads_binsignal <- astads_binsignal %>% 
    #     dplyr::filter(start >= myparams$chromstart) %>% 
    #     dplyr::filter(end <= myparams$chromend) %>%
    #     dplyr::mutate(color = ifelse(pvalue < 0.05, "red", "black"))
    
    y_new <- myplot$y + myplot$height + y_add
    plotSignal(
        data = tad_signal,
        binSize = 50000,
        params = c(myparams),
        fill = "black",
        linecolor = "black",
        x = myplot$x, 
        y = y_new, 
        width = myplot$width,
        height = 1
    )
    return(list(y = y_new))
}


plot_eigens <- function(myplot = NULL, cond=NULL, y_add = NULL, myparams=NULL){
    eigens <- NULL
    ifelse(cond == "as", eigens <- eigens_as, eigens <- eigens_gfp)
    
    # dplyr::filter(start >= myparams$chromstart, end < myparams$chromend)
    y_new <- myplot$y + myplot$height + y_add
    plotSignal(
        data = eigens[[myparams$chrom]] %>%
            dplyr::mutate(chrom = chr, end = end - 1) %>%
            dplyr::select(chrom, start, end, score),
        negData = TRUE,
        #binSize = 100000,
        params = myparams,
        fill = "black",
        linecolor = "black",
        x = myplot$x, 
        y = y_new, 
        width = myplot$width,
        height = 1
    )
    return(list(y = y_new))
}

plot_peaks <- function(myplot = NULL, cond=NULL, y_add = NULL, myparams=NULL){
    peaks <- NULL
    ifelse(cond == "as", peaks <- peaks_as, peaks <- peaks_gfp)
    
    y_new <- myplot$y + myplot$height + y_add   
    plotSignal(
        data = peaks %>%
            dplyr::select(chrom, chromStart, chromEnd, signalValue) %>%
            dplyr::mutate(start = chromStart,
                          end = chromEnd - 1, 
                          score = signalValue) %>%
            dplyr::select(chrom, start, end, score),
        negData = FALSE,
        #binSize = 100000,
        params = c(myparams),
        fill = "black",
        linecolor = "black",
        x = myplot$x, 
        y = y_new, 
        width = myplot$width,
        height = 1
    )
    return(y_new)
}

plot_bigwig <- function(myplot = NULL, cond=NULL, y_add = NULL, myparams=NULL){
    y_new = myplot$y + myplot$height + y_add
    
    bw1 <- NULL
    bw2 <- NULL
    if(cond == "as"){
        bw1 <- bigwig_as1
        bw2 <- bigwig_as2
    } else {
        bw1 <- bigwig_gfp1
        bw2 <- bigwig_gfp2
    }
    
    plotSignal(
        data = bw1,
        params = myparams,
        fill = "black",
        linecolor = "black",
        x = myplot$x, 
        y = y_new, 
        width = myplot$width,
        height = 1
    )
    
    plotSignal(
        data = bw2,
        params = myparams,
        fill = "black",
        linecolor = "black",
        x = myplot$x, 
        y = y_new + 1.2, 
        width = myplot$width,
        height = 1
    )
}


# Read R data objects -------------------------------------------------------------------------

rds_dir <- "./output/R_objects"

tads_as_all <- readRDS("./output/R_objects/ontad_as_tads.rds")
tads_gfp_all <- readRDS("./output/R_objects/ontad_gfp_tads.rds")

cores_all <- readRDS(glue("./output/R_objects/cores_all.rds"))

loops_all <- readRDS(glue("./output/R_objects/loops_all.rds"))

loops_pc <- readRDS("./output/R_objects/a1prom.a2core.loops.any1bp.overlap.rds")
loops_cp <- readRDS("./output/R_objects/a1core.a2prom.loops.any1bp.overlap.rds")

eigens_as <- readRDS("./output/R_objects/as_eigens.rds")
eigens_gfp <- readRDS("./output/R_objects/gfp_eigens.rds")

peaks_as <- readRDS(glue("./output/R_objects/as_peaks.rds")) 
peaks_gfp <- readRDS(glue("./output/R_objects/gfp_peaks.rds")) 

bigwig_as1 <- readBigwig(file="./data/RNAseq/LentiAS_1_rna.bw")
bigwig_as2 <- readBigwig(file="./data/RNAseq/LentiAS_2_rna.bw")
bigwig_gfp1 <- readBigwig(file="./data/RNAseq/LentiGFP_1_rna.bw")
bigwig_gfp2 <- readBigwig(file="./data/RNAseq/LentiGFP_2_rna.bw")

# Read R data objects for each sample chr5, chr20 ---------------------------------------------

# tads_as2_all <- readRDS(file="./output/r_objects/ontad_as2_tads.rds")
# tads_as3_all <- readRDS(file="./output/r_objects/ontad_as3_tads.rds")
# tads_gfp2_all <- readRDS(file="./output/r_objects/ontad_gfp2_tads.rds")
# tads_gfp3_all <- readRDS(file="./output/r_objects/ontad_gfp3_tads.rds")

# Call plot functions -------------------------------------------------------------------------

make_plots <- function(cond=NULL, samp=NULL, mypal=NULL, z=NULL, myparams=NULL){
    output_hic <- NULL
    if(cond=="as"){
        output_hic <- plot_hic(cond=cond, mypal=mypal, z=z, myparams=myparams, x_add=8, y_add=0)
    }else{
        output_hic <- plot_hic(cond=cond, mypal=mypal, z=z, myparams=myparams, x_add=0, y_add=0)
    }
    
    # print(output_hic$plot$zrange)
    output_genes <- plot_genes(myplot=output_hic, y_add = 0.1, myparams=myparams)
    #output_tads <- plot_tads(cond = cond, myplot = output_hic)
    output_ontads_10k <- plot_ontads(cond = cond, myplot = output_hic, res=10000, myparams=myparams)
    #output_ontads_50k <- plot_ontads(cond = cond, myplot = output_hic, res=50000)
    output_cores <- plot_cores(myplot=output_hic, y_add = 0.7, myparams=myparams)
    output_loops <- plot_loops(myplot=output_hic, y_add = 1.8, myparams=myparams)
    #output_tad_signal <- plot_tad_signal(cond=cond, y_add = 2.1)
    output_eigens <- plot_eigens(myplot=output_hic, cond=cond, y_add=2.5, myparams=myparams)
    output_peaks <- plot_peaks(myplot=output_hic,cond=cond, y_add=3.7, myparams=myparams)
    output_bigwig <- plot_bigwig(myplot=output_hic,cond=cond, y_add=4.9, myparams=myparams)
    output_glabel <- plot_glabel(myplot = output_hic, y_add = 7.1, fsize = 20)
    #pageGuideHide()
}

# Mtor: chr5:165,262,812-165,374,967
# Egr2: chr20:22,451,169-22,462,018
# Jun: chr5:114,010,183-114,015,277
# Junb: chr19:26,091,971-26,095,756
# Jund: chr16:20,484,027-20,487,707

gene_coords <- data.frame(gene = character(), chrom = character(), start = numeric(), end = numeric())
gene_coords <- gene_coords %>% 
    tibble::add_row(gene = "Mtor", chrom = "chr5", start = 165261818, end = 165373967) %>%
    tibble::add_row(gene = "Egr2", chrom = "chr20", start = 22454461, end = 22458753) %>%
    tibble::add_row(gene = "Jun", chrom = "chr5", start = 114011186, end = 114014277) %>%
    tibble::add_row(gene = "Junb", chrom = "chr19", start = 26092972, end = 26094756) %>%
    tibble::add_row(gene = "Jund", chrom = "chr16", start = 20485028, end = 20486707) %>%
    tibble::add_row(gene = "Fosl1", chrom = "chr1", start = 220826560, end = 220835066) %>%
    tibble::add_row(gene = "Fos", chrom = "chr6", start = 109300433, end = 109303299) %>%
    tibble::add_row(gene = "Fosl2", chrom = "chr6", start = 25598936, end = 25616995) %>%
    tibble::add_row(gene = "Fosb", chrom = "chr1", start = 80214691, end = 80221417)

genes_list <- list()
for(mygene in gene_coords$gene){
    gene_row <- gene_coords %>% dplyr::filter(gene==mygene)
    mychrom <- gene_row$chrom
    mystart <- gene_row$start
    myend <- gene_row$end
    genes_list[[mygene]] <- data.frame(chrom = mychrom, start = mystart, end = myend) %>%
        tibble::add_row(chrom = mychrom, start = mystart-1e5, end = myend+1e5) %>%
        tibble::add_row(chrom = mychrom, start = mystart-1e6, end = myend+1e6)
}

genes_list$Mtor <- genes_list$Mtor %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-6.5e5, end = .$end[1]+6e5)
genes_list$Egr2 <- genes_list$Egr2 %>%
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-6e5, end = .$end[1]+9e5)
genes_list$Jun <- genes_list$Jun %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-8e5, end = .$end[1]+1.3e6)
genes_list$Junb <- genes_list$Junb %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-1.2e6, end = .$end[1]+3e5)
genes_list$Jund <- genes_list$Jund %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-7.4e5, end = .$end[1]+5.1e5)
genes_list$Fosl1 <- genes_list$Fosl1 %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-5e5, end = .$end[1]+5e5)
genes_list$Fos <- genes_list$Fos %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-4e5, end = .$end[1]+3.1e5)
genes_list$Fosl2 <- genes_list$Fosl2 %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-3e5, end = .$end[1]+4e5)
genes_list$Fosb <- genes_list$Fosb %>% 
    tibble::add_row(chrom = .$chrom[1], start = .$start[1]-3e5, end = .$end[1]+3.6e5)


get_params <- function(gene=NULL, idx=NULL){
    pgParams(chrom = genes_list[[gene]][idx,]$chrom,
             chromstart = genes_list[[gene]][idx,]$start,
             chromend = genes_list[[gene]][idx,]$end,
             assembly = "rn6")
}


reds <- c("#FFF3F3","#FFE0DC","#FFD8D5","#EE1A03") 
blues <- c("#E1DAF5","#D8C9FF", "#BED8FF", "#0943F5")

# gene="Egr2"; idx=3
plot_fig <- function(gene=NULL, idx=NULL){
    fig_params <- get_params(gene=gene, idx=idx)
    #pdf(file = glue::glue("{outdir}/3Mar2022_FINAL_figure_{gene}_{fig_params$chrom}_{fig_params$chromstart}_{fig_params$chromend}.pdf"), width = 17, height = 11)
    
    pageCreate(width = 17, height = 12, default.units = "inches")
    make_plots(cond="gfp", mypal=colorRampPalette(blues, bias=1), myparams=fig_params)
    make_plots(cond="as", mypal=colorRampPalette(reds, bias=1), myparams=fig_params)
    pageGuideHide()
    
    #dev.off()
}

# Plot all genes AS and GFP
# for(gene in names(genes_list)){
#     for(idx in 1:4){
#         plot_fig(gene=gene, idx=idx)
#     }
# }

# Plot Egr2 TAD levels
# plot_fig(gene="Egr2", idx=3)

# Plot Mtor with bigwig
plot_fig(gene="Mtor", idx=4)
plot_fig(gene="Mtor", idx=3)
plot_fig(gene="Mtor", idx=2)
plot_fig(gene="Mtor", idx=1)

plot_fig_mtor_loops <- function(cond=NULL, mypal=NULL){
    myparams <- pgParams(chrom = "chr5", chromstart = 164840000, chromend = 165390000, assembly = "rn6")
    pdf(glue::glue("{outdir}/Mar10_Mtor_3loops_{cond}_{myparams$chrom}_{myparams$chromstart}_{myparams$chromend}.pdf"), width = 9, height = 7)
    pageCreate(width = 9, height = 7, default.units = "inches")
    output_hic <- plot_hic(cond=cond, mypal=mypal, z=NULL, myparams=myparams, x_add=0, y_add=0)
    output_genes <- plot_genes(myplot=output_hic, y_add = 0.1, myparams=myparams)
    output_ontads_10k <- plot_ontads(cond=cond, myplot = output_hic, res=10000, myparams=myparams, level1_flag=TRUE)
    output_cores <- plot_cores(myplot=output_hic, y_add = 0.7, myparams=myparams)
    output_loops <- plot_loops(myplot=output_hic, y_add = 1.8, myparams=myparams, mtor_flag=TRUE)
    output_glabel <- plot_glabel(myplot = output_hic, y_add = 2.5, fsize = 20)
    pageGuideHide()
    dev.off()
}
plot_fig_mtor_loops(cond="gfp", mypal=colorRampPalette(blues, bias=1))
plot_fig_mtor_loops(cond="as", mypal=colorRampPalette(reds, bias=1))
# saveRDS(genes_list, file = "./output/r_objects/figures_genes_list.rds")

# plot_level_tads <- function(myplot=NULL, gr_tads=NULL, level=NULL, color=NULL){
#     
# }


plot_ontads_samp <- function(cond=NULL, samp=NULL, myplot=NULL, res=NULL, myparams=NULL){
    tads <- NULL
    if(samp == "as2"){
        tads <- tads_as2_all[[myparams$chrom]] } 
    else if(samp == "as3"){ 
        tads <- tads_as3_all[[myparams$chrom]] } 
    else if (samp == "gfp2"){ 
        tads <- tads_gfp2_all[[myparams$chrom]] }
    else if (samp == "gfp3"){ 
        tads <- tads_gfp3_all[[myparams$chrom]] }
    
    # gr_tads <- GRanges(seqnames = myparams$chrom,
    #                    ranges = IRanges(start = tads$startpos * res,
    #                                     end = tads$endpos * res))
    # tads_to_plot <- gr_tads 
    # annoDomains(plot = myplot$plot,
    #             data = tads_to_plot,
    #             linecolor = "black")
    
    tad_level_colors <- c(
        "black", "green3", "red", "blue", "turquoise1",
        "black", "yellow2", "magenta", "tan1")
    
    for(level in tads$TADlevel %>% unique()){
        level_tads <- tads %>% dplyr::filter(TADlevel == level) 
        gr_tads <- GRanges(seqnames = myparams$chrom, 
                           ranges = IRanges(start = level_tads$startpos * res, end = level_tads$endpos * res))
        
        level_color <- tad_level_colors[level + 1]
        annoDomains(plot = myplot$plot,
                    data = gr_tads,
                    linecolor = level_color)
        
        #plot_level_tads(myplot=myplot, gr_tads=gr_tads, color=level_color)
    }
}

make_plots_samp <- function(cond=NULL, z=NULL, samp=NULL, mypal=NULL, myparams=NULL){
    output_hic <- NULL
    if(samp=="as2"){
        output_hic <- plot_hic(cond=cond, samp=samp, mypal=mypal, z=z, myparams=myparams, x_add=0, y_add=0)
    } else if (samp == "as3"){
        output_hic <- plot_hic(cond=cond, samp=samp, mypal=mypal, z=z, myparams=myparams, x_add=8, y_add=0)
    } else if (samp == "gfp2") {
        output_hic <- plot_hic(cond=cond, samp=samp, mypal=mypal, z=z, myparams=myparams, x_add=0, y_add=6)
    } else if (samp == "gfp3") {
        output_hic <- plot_hic(cond=cond, samp=samp, mypal=mypal, z=z, myparams=myparams, x_add=8, y_add=6)
    }
    
    output_genes <- plot_genes(myplot=output_hic, y_add = 0.1, myparams=myparams)
    output_ontads_samp <- plot_ontads_samp(cond = cond, samp=samp, myplot = output_hic, res=10000, myparams=myparams)
    output_glabel <- plot_glabel(myplot = output_hic, y_add = 1.5, fsize = 20)
}


plot_fig_samp <- function(gene=NULL, idx=NULL){
    fig_params <- get_params(gene=gene, idx=idx)
    pdf(file = glue::glue("{outdir}/Mar6_figure_allsamples_{gene}_tad_levels_{fig_params$chrom}_{fig_params$chromstart}_{fig_params$chromend}.pdf"), width = 17, height = 11)
    
    pageCreate(width = 17, height = 12, default.units = "inches")
    make_plots_samp(cond="as", samp="as2", mypal=colorRampPalette(reds, bias=1), myparams=fig_params)
    make_plots_samp(cond="as", samp="as3", mypal=colorRampPalette(reds, bias=1), myparams=fig_params)
    make_plots_samp(cond="gfp", samp="gfp2", mypal=colorRampPalette(blues, bias=1), myparams=fig_params)
    make_plots_samp(cond="gfp", samp="gfp3", mypal=colorRampPalette(blues, bias=1), myparams=fig_params)
    pageGuideHide()
    
    dev.off()
}
# gene="Egr2"; idx=4

# plot_fig_samp(gene="Egr2", idx=4)
# plot_fig_samp(gene="Mtor", idx=4)






# Check OK: 15 hierachical tads in plot and in tibble 15 x 7
# tads_as_all$chr5 %>% 
#     dplyr::mutate(startpos = .$startpos*10000, endpos = .$endpos*10000) %>%
#     dplyr::filter(startpos >= 164600000) %>%
#     dplyr::filter(endpos <= 165960000) 




# make_plots(cond="as", mypal=colorRampPalette(rev(brewer.pal(n = 9, "RdBu")))) 
# make_plots(cond="gfp", mypal=colorRampPalette(brewer.pal(n = 9, "RdBu"))) 
# make_plots(cond="as", mypal=colorRampPalette(brewer.pal(n = 9, "Reds")), z = c(0.0001,247.6424) )
# make_plots(cond = "as", mypal = colorRampPalette(brewer.pal(n = 9, "Reds")), z = c(1,200)) 
# make_plots(cond = "gfp", mypal = colorRampPalette(brewer.pal(n = 9, "Blues")), z = c(1,200)) 




# TODO: go to B compartments -> get # of tads in B compartments
# based on coordinates -> no hierachical tads? checks
# genome wide, and these specific gene regions
# quantify that COREs are at boundaries 
# need to check mean narrowPeak signal in these regions.
# do for each 
# if mean is greater, then that compartment is A, ow it's B

# Egr how close are to TAD boundaries
# midpoint of COREs to TAD boundaries
# for all COREs 
# get start and end coordinates of TAD boundaries 
# CORE midpoints tend to be close to TAD boundaries 

# wilcoxon for alL_chroms don't do mean of mean, just combine all info
# wilcoxon for some of these gene regions, essentially a subset of each chrom
