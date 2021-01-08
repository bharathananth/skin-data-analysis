library(limma)
library(magrittr)
library(tidyverse)
library(hgug4112a.db)
library(ggrepel)
library(ggforce)
library(statmod)
library(GO.db)
library(tibble)

##--------------Parameters----------------------------------------------------------------##
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 

setwd("~/Documents/POSTDOC/skin-data-analysis")
if (!file.exists("data/raw_MA_files.RDS")){
  ##--------------Read in image files-------------------------------------------------------##
  files <- list.files("/extra/Skin Data/raw data/", full.names = TRUE)
  m <- regexpr("(D|E)\\d+_P\\d+", files, perl = TRUE)
  names <- regmatches(files, m)
  images <- read.maimages(files = files, source = "agilent", green.only = TRUE, names = names, other.columns = "gIsWellAboveBG")
  
  ##--------annotation of probes------------------------------------------------------------##
  images$genes$Symbol <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "SYMBOL")
  images$genes$ENSEMBL <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "ENSEMBL")
  images$genes$EntrezID <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "ENTREZID")
  
  saveRDS(images, file = "data/raw_MA_files.RDS",compress = "gzip")
} else {
  images <- readRDS("data/raw_MA_files.RDS")
}

##--------background correct -------------------------------------------------------------##
y <- limma::backgroundCorrect(images, method = "normexp") #nomexp when background intensities are not known

##--------normalize between the different arrays------------------------------------------##
y <- normalizeBetweenArrays(y, method = "quantile")

##-----------------------------filtering out controls and lowly expressed----------------##
Control <- y$genes$ControlType==1L
#Control_M <- y$genes$ControlType==1 #This would also work

NoID <- is.na(y$genes$ENSEMBL)

IsExpr <- rowSums(y$other$gIsWellAboveBG>0) >= 77 #Tradeoff strictness-results (I could choose another threshold)
  #half of samples from each patient could be another option
  #TO REVISIT LATER ON -> make it a vble at beginning (also easier for reader later on)

y0 <- y[!Control & !NoID & IsExpr, ]   ## This data (expressed, identified, not controls) already has a gene annotation in it.
y0$genes <- y0$genes[, c("ProbeName", "Symbol", "ENSEMBL", "EntrezID")] #keep only these columns (ENSEMBL, EntrezID useful, SYMBOL is more humanly-readable)

yave <- avereps(y0, y0$genes[, "ENSEMBL"])  ## Averaging probes mapping to the same gene, group of probes "converted to" probeset)
rownames(yave$E) <- yave$genes$ProbeName

saveRDS(yave, file = "visualize/data/rawdata.rds")

##----------------------extracting sample details from column names ----------------------##
experiment <- data.frame(tissue = character(), time = integer(), subject = character()) %>%
              {strcapture("(\\w)(\\d+)_(\\w+)", colnames(y0$E), ., perl = TRUE)}

tissue <- factor(experiment$tissue)
time <- experiment$time
subject <- factor(experiment$subject)

saveRDS(experiment, "visualize/data/experiment.rds")

##--------------------constructing the multi-level design matrix-------------------------##
inphase <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #~ to fourier, to extract freq. info

##--------------------tissue-level analysis--------------------------------------------##
##---------are there rhythmics in at least tissue -------------------------------------##
##-------ANOVA approach (good for balanced designs)-------------------------------------##
##-------this looks for average population rhythms in a tissue--------------------------##
# This is a key part of the analysis. My results are going to depend on my hypothesis testing.
# I define how I will perform hypothesis testing in the following lines
  #model.matrix(~ subject + tissue) %>% head() -> my referenec is p100+DERMIS (H0 is no difference between subjects and p100Dermis)
  #model.matrix(~ 0+subject + tissue + inphase + outphase) -> subject, tissue + time are indep
design <- model.matrix(~ 0 + subject + tissue + tissue:inphase + tissue:outphase) 
# colon means interaction between categorical variables -> here, I say that the coefficients of inphase/outphase should be different for each tissue
# sex could be added here (the more terms in model.matrix, the more coefficients, also pvals will worsen)

wts <- arrayWeights(yave, model.matrix(~tissue + subject)) # simpler model used as suggested in userguide. (avoid?)
fit <- lmFit(yave, design, weights = wts) #exactly like lm, just 'better' pvals given in the next eBayes
  #H0 of a lm is that slopes=0

fit2 <- eBayes(fit, trend = TRUE, robust = TRUE) #calculates the pvals (std step in userguide)

rhy_indices <- which(grepl("phase",colnames(design)))
results <- topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
  set_colnames(gsub("\\.","_", colnames(.))) #for which genes is tissue:inphase and tissue:outphase != 0 (not caring about subject or tissue differences at this point)

saveRDS(results, file = "visualize/data/results.rds")

results %<>% #results %<>% == results <- results %>%
      dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #amplitude_D == cos_component**2 + sin_component**2, in log2 values??
             A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
             phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments, y and x, arctan just takes the angle
             phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)

{ggplot(data = dplyr::filter(results, pmax(A_D, A_E) > amp_cutoff) %>% dplyr::mutate(len = length(adj_P_Val))) + 
      #pmax gives the maximum of each zipped pair of values A_D[i], A_E[i]
  stat_ecdf(aes(x = adj_P_Val, len=len, y = ..y..*len), geom = "step") +  #we don't want fraction of genes (default of ecdf) but number, that's why we *len
  coord_cartesian(xlim = c(0.001, 0.05), ylim = c(0,2000)) + theme_bw() +
  xlab("False discovery rate") + ylab("number of rhythmic genes")} %T>% ggsave("figures/Noofgenes.pdf", .) #T pipe works like %>% but it returns the left-hand side instead of the right-hand

some_rhythm <- dplyr::filter(results, pmax(A_D, A_E) > amp_cutoff & adj_P_Val < fdr_cutoff) %$% ProbeName #some_rhythm in both A_D & A_E, note that amp is in log2?

clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")

{ggplot(data = dplyr::filter(results, ProbeName %in% some_rhythm), aes(x=2^A_D-1, y=2^A_E-1)) + #amps are in log2
  geom_point(aes(color=Symbol %in% clock_genes), size=1) + 
  geom_text_repel(aes(label=ifelse(Symbol %in% clock_genes, as.character(Symbol), "")), color="red", max.overlaps=Inf, box.padding = 0.5) + 
  scale_color_manual(values=c("black","red")) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none") + coord_cartesian(xlim=c(0,1), ylim = c(0,1)) +
  xlab("Relative Amplitude -- Dermis") + ylab("Relative Amplitude -- Epidermis")} %T>% ggsave("figures/Amplitude_compare.pdf", .)

{ggplot(data = dplyr::filter(results, ProbeName %in% some_rhythm), aes(x=phaseD, y=phaseE)) + 
  geom_point(aes(color=Symbol %in% clock_genes), size=1) + 
  geom_text_repel(aes(label=ifelse(Symbol %in% clock_genes, as.character(Symbol), "")), color="red", max.overlaps=Inf, box.padding = 0.5) + 
  scale_color_manual(values=c("black","red")) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none")  + geom_abline(slope = 1, alpha=0.2) + 
  scale_x_continuous(breaks = seq(-12,12,by = 6)) +
  scale_y_continuous(breaks = seq(-12,12,by = 6)) + 
  xlab("Phase -- Dermis") + ylab("Phase -- Epidermis")} %T>% ggsave("figures/Phase_compare.pdf", .)

##--------------Selecting only genes rhythmic somewhere--------------------------------##
# GO analysis of rhythmic genes (amp>, fdr< cutoff)
g <- goana(yave[yave$genes$ProbeName %in% some_rhythm,]$genes$EntrezID, universe = yave$genes$EntrezID) #universe is the 'background' universe for the analysis, if none provided, then all EntrezID will be used
topGO(g, n=20, truncate.term = "50") #body fluid secretion top term

# Pathway analysis of rhythmic genes (amp>, fdr< cutoff)
k <- kegga(yave[yave$genes$ProbeName %in% some_rhythm,]$genes$EntrezID, universe = yave$genes$EntrezID)
topKEGG(k, n=20, truncate.path = "50")

##------------working only with these rhythmic genes-----------------------------------##
##------------performing differential rhythmicity analysis---------------------------- ##
yrhy <- yave[yave$genes$ProbeName %in% some_rhythm, ] #expression levels of rhythmic genes
fit_rhy <- fit[yave$genes$ProbeName %in% some_rhythm, ] #fits of those rhythmic genes

colnames_proper <- gsub(":", "_", colnames(design))
contrast1 <- makeContrasts(tissueD_inphase - tissueE_inphase, tissueD_outphase - tissueE_outphase, levels = colnames_proper) #####?
rownames(contrast1) <- colnames(design)

fit_rhy2 <- contrasts.fit(fit = fit_rhy, contrasts = contrast1)
fit_rhy2 <- eBayes(fit_rhy2, trend = TRUE, robust = TRUE)

results2 <- topTable(fit_rhy2, number = Inf, sort.by = "none") %>%
  set_colnames(gsub("\\.","_", colnames(.)))

diff_rhythm <- dplyr::filter(results2, adj_P_Val < fdr_cutoff) %$% ProbeName

results_diff_rhythm <- results %>% dplyr::filter(ProbeName %in% diff_rhythm) %>%
  dplyr::mutate(Aratio = A_E-A_D, ph_diff = (phaseE-phaseD + 12) %% 24 - 12)

{ggplot(dplyr::filter(results, ProbeName %in% some_rhythm)) + 
  geom_rect(xmin = amp_cutoff, xmax=1, ymin=0, ymax = amp_cutoff, alpha=0.1, fill="grey90") + 
  geom_rect(xmin = 0, xmax=amp_cutoff, ymin=amp_cutoff, ymax = 1, alpha=0.1, fill="grey90") + 
  geom_point(aes(x=A_D, y =A_E, color=factor(ProbeName %in% diff_rhythm, levels=c(TRUE, FALSE))), size=0.75) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none") + scale_color_manual(values=c("red", "black")) + 
  xlab("Relative Amplitude -- Dermis") + ylab("Relative Amplitude -- Epidermis") + 
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))} %T>% ggsave("figures/Amplitude_compare_diffrhy.pdf", .)

ggplot(results_diff_rhythm) + geom_histogram(aes(x=Aratio),bins = 50, fill=NA, color="grey20") + theme_bw() +
  xlab("Amplitude_Epidermis - Amplitude_Dermis")

ggplot(results_diff_rhythm) + geom_histogram(aes(x=ph_diff),bins = 50, fill=NA, color="grey20") + theme_bw() +
  xlab("Phase_Epidermis - Phase_Dermis") + coord_polar(start = pi/50+pi)

ggplot(results_diff_rhythm) + geom_point(aes(x=Aratio, y=ph_diff))

fitted_values <- fit$coefficients[, rhy_indices] %*% t(design[, rhy_indices]) %>%
                 set_rownames(fit$genes$Symbol) %>%
                 set_colnames(colnames(yave$E))

in_both_rhythm <- results_diff_rhythm %>% dplyr::filter(pmin(A_D, A_E) > amp_cutoff)

mean_data <- yave$E %>% transform(ProbeName = yave$genes$ProbeName,
                     Symbol = yave$genes$Symbol) %>% as_tibble() %>%
         tidyr::gather(junk, value, -ProbeName, -Symbol) %>%
         tidyr::separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
         tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
         dplyr::filter(ProbeName %in% diff_rhythm) %>%
         dplyr::group_by(ProbeName, Symbol, tissue, time) %>%
         dplyr::summarise(value = mean(value)) %>%
         dplyr::mutate(value = value - mean(value))

page_no <- 8 #BA had 8 here, but I see 8 pages with n_pages (https://www.programmingwithr.com/how-to-make-your-facet-wrap-facet-grid-ggplots-span-across-multiple-pages-in-pdf/)
             #somehow I go back to 8 because 9 gives me problems when plotting
#{ggplot(mean_data) + geom_rect(xmin=4,xmax=36,ymin=-amp_cutoff,ymax=amp_cutoff, fill="grey80", alpha=0.1) +
#      geom_line(aes(x=time, y=value, group=tissue, color=tissue)) + 
#      facet_wrap_paginate(~Symbol, ncol=8, nrow = 9, page=page_no) + theme_bw(base_size = 8) + 
#      theme(aspect.ratio = 1, strip.text = element_text(size=6), legend.position = "top")
#} %>% ggsave(sprintf("figures/Diff_rhy_genes_%d.pdf", page_no),.,width = 7,height = 10) 
# this was BA's approach, but didn't work -> see below

# my workaround
for (i in 1:page_no){
  ggplot(mean_data) + geom_rect(xmin=4,xmax=36,ymin=-amp_cutoff,ymax=amp_cutoff, fill="grey80", alpha=0.1) +
    geom_line(aes(x=time, y=value, group=tissue, color=tissue)) + 
    facet_wrap_paginate(~Symbol, ncol=8, nrow = 9, page=i) + theme_bw(base_size = 8) + 
    theme(aspect.ratio = 1, strip.text = element_text(size=6), legend.position = "top") -> p
  p %>% ggsave(sprintf("figures/Diff_rhy_genes_%d.pdf", i),.,width = 7,height = 10)
}

##-----random effect approach (more powerful for unbalanced designs)--------------------##
# design3 <- model.matrix(~ 0 + tissue + tissue:inphase + tissue:outphase)
# 
# dupcor <- duplicateCorrelation(yave, design, block = experiment$subject)
# 
# fit <- lmFit(yave, design3, block = experiment$subject, correlation = dupcor$consensus.correlation)
# 
# fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)
# 
# results2 <- topTable(fit2, coef = 3:6, number = Inf, sort.by = "none")

#--------------------------------------
#--------------------------------------

### TO BE REVISED AFTER BA'S MEETING, I HAVE TO INCLUDE GENDER IN MY MODEL MATRIX PROBABLY
# Are there differences between dermis and epidermis?
rhythm_E <- dplyr::filter(results, A_E > amp_cutoff & adj_P_Val < fdr_cutoff) %$% ProbeName 
rhythm_D <- dplyr::filter(results, A_D > amp_cutoff & adj_P_Val < fdr_cutoff) %$% ProbeName 

# GO terms
g_E <- goana(yave[yave$genes$ProbeName %in% rhythm_E,]$genes$EntrezID, universe = yave$genes$EntrezID) #universe is the 'background' universe for the analysis, if none provided, then all EntrezID will be used
g_D <- goana(yave[yave$genes$ProbeName %in% rhythm_D,]$genes$EntrezID, universe = yave$genes$EntrezID) #universe is the 'background' universe for the analysis, if none provided, then all EntrezID will be used

topGO(g_E, n=20, truncate.term = "50") #circ. rhythm, detox, metals
topGO(g_D, n=20, truncate.term = "50")

# KEGG pathways
k_E <- kegga(yave[yave$genes$ProbeName %in% rhythm_E,]$genes$EntrezID, universe = yave$genes$EntrezID)
k_D <- kegga(yave[yave$genes$ProbeName %in% rhythm_D,]$genes$EntrezID, universe = yave$genes$EntrezID)

topKEGG(k_E, n=20, truncate.path = "50") #more metabolic-related pathways
topKEGG(k_D, n=20, truncate.path = "50") #lots of disease-related pathways! diabetes, malaria, trypanosomiasis, amoebiasis, HSV, legionella

#--------------------------------------
#--------------------------------------

# Are there differences between males and females? Look for genes in Y chromosome
Ychr_genes <- read.csv("./genes_chrY_bioMart.csv")$Gene.stable.ID #could be done directly with biomaRt in R but biomaRt can't get installed, see https://www.biostars.org/p/3650/

"""
But this below does not account for males and females actually!
I just see the genes related to the Y-chromosome vs the genes that are NOT in the y-chromosome

Bur from yrhy_males$E I can actually see which patients are males -> check males
SOMETHING IS WRONG, i get all patients being male + female xD
"""
### BUT THESE ARE NOT MALES AND FEMALES, they are just y-related genes vs not-y-related genes
yrhy_males   <- yave[yave$genes$ENSEMBL %in% Ychr_genes, ] 
yrhy_females <- yave[!yave$genes$ENSEMBL %in% Ychr_genes, ] 

all     <- unlist(strsplit(colnames(yrhy$E), "_"))[seq(2, length(unlist(strsplit(colnames(yrhy$E), "_"))), 2)] %>% unique
males   <- unlist(strsplit(colnames(yrhy_males$E), "_"))[seq(2, length(unlist(strsplit(colnames(yrhy_males$E), "_"))), 2)] %>% unique
females <- unlist(strsplit(colnames(yrhy_females$E), "_"))[seq(2, length(unlist(strsplit(colnames(yrhy_females$E), "_"))), 2)] %>% unique

fit_rhy_males   <- fit[yave$genes$ENSEMBL %in% Ychr_genes, ]
fit_rhy_females <- fit[!yave$genes$ENSEMBL %in% Ychr_genes, ]


fit_rhy2_males   <- contrasts.fit(fit = fit_rhy_males, contrasts = contrast1)
fit_rhy2_females <- contrasts.fit(fit = fit_rhy_females, contrasts = contrast1)
fit_rhy2_males   <- eBayes(fit_rhy2_males, trend = TRUE, robust = TRUE)
fit_rhy2_females <- eBayes(fit_rhy2_females, trend = TRUE, robust = TRUE)

results_males <- topTable(fit_rhy2_males, number = Inf, sort.by = "none") %>%
  set_colnames(gsub("\\.","_", colnames(.)))
results_females <- topTable(fit_rhy2_females, number = Inf, sort.by = "none") %>%
  set_colnames(gsub("\\.","_", colnames(.)))