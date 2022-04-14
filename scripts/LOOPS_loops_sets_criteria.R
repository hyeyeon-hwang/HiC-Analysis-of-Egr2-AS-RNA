library(glue)
library(magrittr)
library(readr)
library(optparse)
library(kableExtra)

args_list = list(
    make_option(c("-d", "--diffdir"),
        type="character",
        default=NULL,
        help="Path to data directory of *txt.gz files of differential loops for all chromosomes",
        metavar="character"),
    make_option(c("-s","--sigloops"),
        type="character",
        default=NULL,
        help="Path to significant loops file previously generated",
        metavar="character"),
    make_option(c("-o","--outdir"),
        type="character",
        default=NULL,
        help="Path to output directory for loops sets",
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
    make_option(c("-p","--pval"),
        type="numeric",
        default=NULL,
        help="p-value threshold for differential and significant loops",
        metavar="numeric"),
    make_option(c("-l","--log2fc"),
        type="numeric",
        default=NULL,
        help="Log2 fold change threshold for user-defined differential loops",
        metavar="numeric")
)


args_parser = OptionParser(option_list=args_list)
args = parse_args(args_parser)
# print(glue("args = {args}"))


findDiffLoops <- function(
    indir = args$diffdir, 
    outdir = args$outdir, 
    distance = args$distrange, 
    sigindices = args$sigloops %>% data.table::fread(), 
    log2fc.cutoff = args$log2fc, 
    pval.cutoff = args$pval,
    diffChr = c(1:20, 'X')){
 
    # To extract pval and log2fc of sig and diff loops --------------------------------------------
    # indir = "./output/LOOPS/differential_loops"
    # outdir = "./output/LOOPS/pval_logfc"
    # distance = 10000
    # sigindices = "./output/LOOPS/significant_loops/10000_50Kb-2Mb_analysis_sigindices.txt.gz" %>% data.table::fread()
    # log2fc.cutoff = 1
    # pval.cutoff = 0.05
    # diffChr = c(1:20, 'X')
    
    outdir.pval.logfc <- "./output/LOOPS/pval_logfc"
    diffloops.files <- sapply(diffChr, function(i) {glue('{indir}/diff_resASoverGFP_chr{i}.txt.gz')}) %>% as.character()
  
    # log2FoldChange is not 0
    diffloops <- data.frame() # 134503 
    for (f in diffloops.files){ diffloops <- rbind(diffloops, f %>% data.table::fread())}
    
    dloops.gain <- diffloops[(diffloops$log2FoldChange > log2fc.cutoff) & (diffloops$pvalue < pval.cutoff),] # 104
    dloops.lost <- diffloops[(diffloops$log2FoldChange < -log2fc.cutoff) & (diffloops$pvalue < pval.cutoff),] # 441
    
    dloops.static.up <- diffloops[(diffloops$log2FoldChange > 0) & (diffloops$log2FoldChange <= log2fc.cutoff),] # 42705
    dloops.static.gain.notsig <- diffloops[(diffloops$log2FoldChange > log2fc.cutoff) & (diffloops$pvalue >= pval.cutoff),] # 8561
    dloops.static.down <- diffloops[(diffloops$log2FoldChange < 0) & (diffloops$log2FoldChange >= -log2fc.cutoff),] # 61893
    dloops.static.lost.notsig <- diffloops[(diffloops$log2FoldChange < -log2fc.cutoff) & (diffloops$pvalue >= pval.cutoff),] # 20799
    # dloops.static total: 133958 x 9
    dloops.static <- dloops.static.up %>%
        dplyr::bind_rows(dloops.static.gain.notsig) %>%
        dplyr::bind_rows(dloops.static.down) %>%
        dplyr::bind_rows(dloops.static.lost.notsig) # 133958 x 9
    
    int.dstatic <- dloops.static %>% tibble::add_column(type = 'static', .before = 1) # 133958
    int.gain <- dloops.gain %>% tibble::add_column(type = 'gain', .before = 1) # 104
    int.lost <- dloops.lost %>% tibble::add_column(type = 'lost', .before = 1) # 441
 
    #dloops.static.coords <- dloops.static %>% dplyr::select(chr, startI, startJ) # 133958 x 3
    sloops.nodloops <- dplyr::setdiff(sigindices, diffloops %>% dplyr::select(chr, startI, startJ)) # 53675 x 3 
    sdloops.static.total <-  dloops.static %>%
        dplyr::select(chr, startI, startJ) %>%
        dplyr::bind_rows(sloops.nodloops) # 187633 x 3 (= 188178 - (104 + 441))
  
    int.sstatic <- sloops.nodloops %>% 
        tibble::add_column(type = 'static', .before = 1) %>%
        tibble::add_column(log2FoldChange = 0) %>%
        tibble::add_column(pvalue = NA) # 53675 x 6
  
    int.all <- int.gain %>% dplyr::select(type, chr, startI, startJ, log2FoldChange, pvalue) %>%
        dplyr::bind_rows(int.lost %>% dplyr::select(type, chr, startI, startJ, log2FoldChange, pvalue)) %>%
        dplyr::bind_rows(int.dstatic %>% dplyr::select(type, chr, startI, startJ, log2FoldChange, pvalue)) %>%
        dplyr::bind_rows(int.sstatic) # 188178 x 6
    
    readr::write_tsv(int.all, glue('{outdir.pval.logfc}/int.all.tsv'))
    
    valid_chrs <- c(1:20, 'X') %>% sapply(., function(i) {glue('chr{i}')}) %>% as.character()
    for(chrm in valid_chrs){
        int.all %>% 
            filter(chr == chrm) %>% 
            readr::write_tsv(., glue('{outdir.pval.logfc}/int.{chrm}.tsv'))
        int.all %>% 
            filter(chr == chrm) %>%
            filter(!is.na(pvalue)) %>%
            readr::write_tsv(., glue('{outdir.pval.logfc}/int.noNA.{chrm}.tsv'))
  }
  
  # a1prom.a2core.loops 
  # a1core.a2prom.loops
  
  int.pc <- a1prom.a2core.loops %>% dplyr::select(chr, startI, startJ, loopset) # 90
  int.cp <- a1core.a2prom.loops %>% dplyr::select(chr, startI, startJ, loopset) # 111
  int.pc.cp <- int.pc %>% dplyr::bind_rows(int.cp) # 201
  int.pc.cp.uniq <- int.pc.cp %>% unique() # 192 
  # This means that 9 loops had promoter and core in both anchors
  # even fewer loop anchors

  dplyr::intersect(int.all %>% dplyr::select(chr, startI, startJ), 
                   int.pc.cp.uniq %>% dplyr::select(chr, startI, startJ)) # 192, which we want
  
  int192.df <- data.frame(type = character(), 
                          chr = character(),
                          startI = numeric(),
                          startJ = numeric(),
                          log2FoldChange = numeric(),
                          pvalue = numeric())
  for(i in 1:nrow(int.pc.cp.uniq)) {
       #i = 11 # 13
       loopi <- int.pc.cp.uniq[i,]
       #print(i)
       #print(loopi$loopset)
       loop.add <- int.all %>% 
            dplyr::filter(chr == loopi$chr) %>%
            dplyr::filter(startI == loopi$startI) %>%
            dplyr::filter(startJ == loopi$startJ)
       #print(loop.add$type)
       int192.df <- int192.df %>% tibble::add_row(loop.add)
  }
  
  int192.df %>% 
       readr::write_tsv(., glue('{outdir.pval.logfc}/int192.tsv'))
  int192.df %>% 
       filter(!is.na(pvalue)) %>%
       readr::write_tsv(., glue('{outdir.pval.logfc}/int192.noNA.tsv'))
  
  # > int.pc.cp.uniq %>% filter(loopset == 'Ldg') %>% nrow() # 8
  # > int.pc.cp.uniq %>% filter(loopset == 'Ldl') %>% nrow() # 16
  # > int.all %>% filter(type == 'gain') %>% nrow() # 104
  # > int.all %>% filter(type == 'lost') %>% nrow() # 441
  
  a1prom.a2core.loops %>% filter(loopset == 'Ldg') # 5
  a1core.a2prom.loops %>% filter(loopset == 'Ldg') # 3
  a1prom.a2core.loops %>% filter(loopset == 'Ldl') # 3
  a1core.a2prom.loops %>% filter(loopset == 'Ldl') # 13
  
  
  # end of recording pval and logfc to {outdir.pval.logfc}
  
  
  diffloops.up <- diffloops[diffloops$log2FoldChange > 0,] # new: 51370, 52737 x 9
  diffloops.down <- diffloops[diffloops$log2FoldChange < 0,] # new: 83133, 84958 x 9
  # diffloops.zero <- diffloops[diffloops$log2FoldChange == 0,] # 0 x 9
  
  # can also do: diffloops[diffloops$log2FoldChange > log2fc.cutoff, ] AND < -log2fc.cutoff
  diffloops.gain <- diffloops.up[diffloops.up$log2FoldChange > log2fc.cutoff,] # #new:8665, 8716
  diffloops.lost <- diffloops.down[diffloops.down$log2FoldChange < -log2fc.cutoff,] # new: 21240, 21361
  
  diffloops.gain.sigp <- diffloops.gain[diffloops.gain$pvalue < pval.cutoff, ] # 104
  diffloops.lost.sigp <- diffloops.lost[diffloops.lost$pvalue < pval.cutoff, ] # 441
  
  # in first that is not in second
  # in diffloops.up but not in diffloops.gain
  diffloops.up.notgain <- dplyr::setdiff(diffloops.up, diffloops.gain) # new: 42705 x 9, 44021 x 9, 44021 + 8716 = 52737
  diffloops.down.notlost <- dplyr::setdiff(diffloops.down, diffloops.lost) # new: 61893 x 9,  63597 x 9, 63597 + 21361 = 84958
  # significant loops that are not in differential loops
  sigloops.notdiff <- dplyr::setdiff(sigindices, diffloops[,1:3]) # new: 53675 x 3, 53798 x 3
  # all differential loops are in the significant loops
  # sigdiff.union <- dplyr::union(sigindices, diffloops[,1:3]) # 191493 x 4 = dim(sigindices)
  
  
  # outdir <- "./output/LOOPS/loops_sets"

  readr::write_tsv(sigindices, glue('{outdir}/significant_loops_all_{distance}.bed'))
  readr::write_tsv(diffloops[,1:3], glue('{outdir}/differential_loops_all_{distance}.bed'))
  readr::write_tsv(diffloops.gain[,1:3], glue('{outdir}/differential_loops_gain_{distance}.bed'))
  readr::write_tsv(diffloops.lost[,1:3], glue('{outdir}/differential_loops_lost_{distance}.bed'))
  # differential but doesn't meet our criteria of diff 
  # significant but not in list of differential (HiC-DC+'s list of differential)
  readr::write_tsv(diffloops.up.notgain[,1:3], glue('{outdir}/differential_loops_up_notgain_{distance}.bed'))
  readr::write_tsv(diffloops.down.notlost[,1:3], glue('{outdir}/differential_loops_down_notlost_{distance}.bed'))
  readr::write_tsv(sigloops.notdiff, glue('{outdir}/significant_loops_not_differential_{distance}.bed'))
  
  readr::write_tsv(diffloops.gain.sigp[,1:3], glue('{outdir}/differential_loops_gain_sigp_{distance}.bed'))
  readr::write_tsv(diffloops.lost.sigp[,1:3], glue('{outdir}/differential_loops_lost_sigp_{distance}.bed'))
  
  
  dimDf <- data.frame(
    Dataset = c(
      glue('significant_loops_all_{distance}.bed'),
      glue('differential_loops_all_{distance}.bed'), 
      glue('significant_loops_not_differential_{distance}.bed'),
      
      glue('differential_loops_gain_{distance}.bed'),
      glue('differential_loops_lost_{distance}.bed'),
      
      glue('differential_loops_up_notgain_{distance}.bed'),
      glue('differential_loops_down_notlost_{distance}.bed'),
      
      glue('differential_loops_gain_sigp_{distance}.bed'),
      glue('differential_loops_lost_sigp_{distance}.bed')),
    
    NumberOfLoops = c(
      sigindices %>% dim() %>% .[1],
      diffloops %>% dim() %>% .[1],
      sigloops.notdiff %>% dim() %>% .[1],
      
      diffloops.gain %>% dim() %>% .[1],
      diffloops.lost %>% dim() %>% .[1],
      
      diffloops.up.notgain %>% dim() %>% .[1],
      diffloops.down.notlost %>% dim() %>% .[1],
      
      diffloops.gain.sigp %>% dim() %>% .[1], 
      diffloops.lost.sigp %>% dim() %>% .[1])
  ) # end: dimDf
  
  return(dimDf)
  
} # end: findDiffLoops()


dimDiffLoops <- findDiffLoops()
#dimDiffLoops.50Kb.2Mb

# > dimDiffLoops.50Kb.2Mb
# Dataset NumberOfLoops
# 1              significant_loops_all_50Kb-2Mb.bed        188178
# 2             differential_loops_all_50Kb-2Mb.bed        134503
# 3 significant_loops_not_differential_50Kb-2Mb.bed         53675
# 4            differential_loops_gain_50Kb-2Mb.bed          8665
# 5            differential_loops_lost_50Kb-2Mb.bed         21240
# 6      differential_loops_up_notgain_50Kb-2Mb.bed         42705
# 7    differential_loops_down_notlost_50Kb-2Mb.bed         61893
# 8       differential_loops_gain_sigp_50Kb-2Mb.bed           104
# 9       differential_loops_lost_sigp_50Kb-2Mb.bed           441

print(dimDiffLoops)

readr::write_tsv(dimDiffLoops, glue('{args$outdir}/loops_sets_table.tsv'))

#dimLoopsTable <- dimDiffLoops %>% 
#  kbl() %>%
#  kable_styling() %>%
#  save_kable(glue('{args$outdir}/loops_sets_table.png'))


