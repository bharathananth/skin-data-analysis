library(limma)
library(magrittr)
library(tidyverse)
library(hgug4112a.db)
library(ggrepel)
library(ggforce)

##--------------Parameters----------------------------------------------------------------##
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2)

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
  readRDS("data/raw_MA_files.RDS")
}

##--------background correct -------------------------------------------------------------##
y <- limma::backgroundCorrect(images, method = "normexp")

##--------normalize between the different arrays------------------------------------------##
y <- normalizeBetweenArrays(y, method = "quantile")

##-----------------------------filtering out controls and lowly expressed----------------##
Control <- y$genes$ControlType==1L

NoID <- is.na(y$genes$ENSEMBL)

IsExpr <- rowSums(y$other$gIsWellAboveBG>0) >= 77

y0 <- y[!Control & !NoID & IsExpr, ]   ## This data already has a gene annotation in it.
y0$genes <- y0$genes[, c("ProbeName", "Symbol", "ENSEMBL","EntrezID")]

yave <- avereps(y0, y0$genes[, "ENSEMBL"])  ## Averaging probes mapping to the same gene
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
outphase <- sin(2*pi*time/24)

##--------------------tissue-level analysis--------------------------------------------##
##---------are there rhythmics in at least tissue -------------------------------------##
##-------ANOVA approach (good for balanced designs)-------------------------------------##
##-------this looks for average population rhythms in a tissue--------------------------##

design <- model.matrix(~ 0 + subject + tissue + tissue:inphase + tissue:outphase)

wts <- arrayWeights(yave, model.matrix(~tissue + subject)) # simpler model used as suggested in userguide.
fit <- lmFit(yave, design, weights = wts)

fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)

rhy_indices <- which(grepl("phase",colnames(design)))
results <- topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>%
           set_colnames(gsub("\\.","_", colnames(.)))

saveRDS(results, file = "visualize/data/results.rds")

results %<>% 
      mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2),
             A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
             phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi,
             phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)

{ggplot(data = filter(results, pmax(A_D, A_E) > amp_cutoff) %>% mutate(len = length(adj_P_Val))) + 
  stat_ecdf(aes(x = adj_P_Val, len = len, y = ..y..*len), geom = "step") +
  coord_cartesian(xlim = c(0.001, 0.05), ylim = c(0,2000)) + theme_bw() +
  xlab("False discovery rate") + ylab("number of rhythmic genes")} %T>% ggsave("figures/Noofgenes.pdf", .)

some_rhythm <- filter(results, pmax(A_D, A_E) > amp_cutoff & adj_P_Val < fdr_cutoff) %$% ProbeName

clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")

{ggplot(data = filter(results, ProbeName %in% some_rhythm), aes(x=2^A_D-1, y=2^A_E-1)) + geom_point(aes(color=Symbol %in% clock_genes), size=1) + 
  geom_text_repel(aes(label=ifelse(Symbol %in% clock_genes, as.character(Symbol), "")), color="red") + scale_color_manual(values=c("black","red")) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none") + coord_cartesian(xlim=c(0,1), ylim = c(0,1)) +
    xlab("Relative Amplitude -- Dermis") + ylab("Relative Amplitude -- Epidermis")} %T>% ggsave("figures/Amplitude_compare.pdf", .)

{ggplot(data = filter(results, ProbeName %in% some_rhythm), aes(x=phaseD, y=phaseE)) + geom_point(aes(color=Symbol %in% clock_genes), size=1) + 
  geom_text_repel(aes(label=ifelse(Symbol %in% clock_genes, as.character(Symbol), "")), color="red") + scale_color_manual(values=c("black","red")) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none")  + geom_abline(slope = 1, alpha=0.2) + scale_x_continuous(breaks = seq(-12,12,by = 6)) +
  scale_y_continuous(breaks = seq(-12,12,by = 6)) + xlab("Phase -- Dermis") + ylab("Phase -- Epidermis")} %T>% ggsave("figures/Phase_compare.pdf", .)

##--------------Selecting only genes rhythmic somewhere--------------------------------##
g <- goana(yave[yave$genes$ProbeName %in% some_rhythm,]$genes$EntrezID, universe = yave$genes$EntrezID)
topGO(g, n=20,truncate.term = "50")

k <- kegga(yave[yave$genes$ProbeName %in% some_rhythm,]$genes$EntrezID, universe = yave$genes$EntrezID)
topKEGG(k, n =20, truncate.path = "50")

##------------working only with these rhythmic genes-----------------------------------##
##------------performing differential rhythmicity analysis---------------------------- ##
yrhy <- yave[yave$genes$ProbeName %in% some_rhythm, ]
fit_rhy <- fit[yave$genes$ProbeName %in% some_rhythm, ]

colnames_proper <- gsub(":", "_", colnames(design))
contrast1 <- makeContrasts(tissueD_inphase - tissueE_inphase, tissueD_outphase - tissueE_outphase, levels = colnames_proper)
rownames(contrast1) <- colnames(design)

fit_rhy2 <- contrasts.fit(fit = fit_rhy, contrasts = contrast1)
fit_rhy2 <- eBayes(fit_rhy2, trend = TRUE, robust = TRUE)

results2 <- topTable(fit_rhy2, number = Inf, sort.by = "none") %>%
  set_colnames(gsub("\\.","_", colnames(.)))

diff_rhythm <- filter(results2, adj_P_Val < fdr_cutoff) %$% ProbeName

results_diff_rhythm <- results %>% filter(ProbeName %in% diff_rhythm) %>%
             mutate(Aratio = A_E-A_D, ph_diff = (phaseE-phaseD + 12) %% 24 - 12)

{ggplot(filter(results, ProbeName %in% some_rhythm)) + 
  geom_rect(xmin = amp_cutoff, xmax=1, ymin=0, ymax = amp_cutoff, alpha=0.1, fill="grey90") + 
  geom_rect(xmin = 0, xmax=amp_cutoff, ymin=amp_cutoff, ymax = 1, alpha=0.1, fill="grey90") + 
  geom_point(aes(x=A_D, y =A_E, color=factor(ProbeName %in% diff_rhythm, levels=c(TRUE, FALSE))), size=0.75) +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none") + scale_color_manual(values=c("red", "black")) + 
  xlab("Relative Amplitude -- Dermis") + ylab("Relative Amplitude -- Epidermis") + coord_cartesian(xlim=c(0,1),ylim=c(0,1))} %T>% ggsave("figures/Amplitude_compare_diffrhy.pdf", .)

ggplot(results_diff_rhythm) + geom_histogram(aes(x=Aratio),bins = 50, fill=NA, color="grey20") + theme_bw() +
  xlab("Amplitude_Epidermis - Amplitude_Dermis")

ggplot(results_diff_rhythm) + geom_histogram(aes(x=ph_diff),bins = 50, fill=NA, color="grey20") + theme_bw() +
  xlab("Phase_Epidermis - Phase_Dermis") + coord_polar(start = pi/50+pi)

ggplot(results_diff_rhythm) + geom_point(aes(x=Aratio, y=ph_diff))

fitted_values <- fit$coefficients[, rhy_indices] %*% t(design[, rhy_indices]) %>%
                 set_rownames(fit$genes$Symbol) %>%
                 set_colnames(colnames(yave$E))

in_both_rhythm <- results_diff_rhythm %>% filter(pmin(A_D, A_E) > amp_cutoff)

mean_data <- yave$E %>% transform(ProbeName = yave$genes$ProbeName,
                     Symbol = yave$genes$Symbol) %>% as_tibble() %>%
         gather(junk, value, -ProbeName, -Symbol) %>%
         separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
         separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
         filter(ProbeName %in% diff_rhythm) %>%
         group_by(ProbeName, Symbol, tissue, time) %>%
         summarise(value = mean(value)) %>%
         mutate(value = value - mean(value))

page_no <- 9
{ggplot(mean_data) + geom_rect(xmin=4,xmax=36,ymin=-amp_cutoff,ymax=amp_cutoff, fill="grey80", alpha=0.1) +
      geom_line(aes(x=time, y=value, group=tissue, color=tissue)) + 
      facet_wrap_paginate(~Symbol, ncol=8, nrow = 9, page = page_no) + theme_bw(base_size = 8) + 
      theme(aspect.ratio = 1, strip.text = element_text(size=6), legend.position = "top")
} %>% ggsave(sprintf("figures/Diff_rhy_genes_%d.pdf", page_no),.,width = 7,height = 10)



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






