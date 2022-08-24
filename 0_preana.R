suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hgug4112a.db))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(PCAtools))


dir.create("figures",showWarnings = FALSE)
dir.create("results",showWarnings = FALSE)

# R graphics stuff
scale_colour_discrete <- function(...) {
  scale_colour_brewer(..., palette="Dark2")
}

theme_custom <- function(base_size = 11, base_family = "Helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) + theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      panel.background = element_blank(),
      panel.spacing.y  = unit(1.5, "lines"),
      
      plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "lines"),
      plot.title = element_text(hjust = 0.5, face="bold", size=10),
      
      strip.background = element_blank(),
      strip.text = element_text(face="bold"),
      
      legend.title = element_blank(),
      legend.position="top",
      
      aspect.ratio = 1
    )
}

theme_set(theme_custom(base_size = 9))
update_geom_defaults("line", list(size = 0.8))


#------------------
#------------------


# 1. Read info of subjects
# ------------------------
info_subjects_long <- read.csv("resources/info_subjects.csv") %>%
  mutate(Bedtime_Work = Bedtime_Work %>% parse_hms(), Sleeptime_Work =  Sleeptime_Work %>% parse_hms(), 
         Wakeuptime_Work = Wakeuptime_Work %>% parse_hms(), Bedtime_Free = Bedtime_Free %>% parse_hms(),
         Sleeptime_Free = Sleeptime_Free %>% parse_hms(), Wakeuptime_Free = Wakeuptime_Free %>% parse_hms(),
         
         age = 2011 - Birth_year,
         Sleepduration_Free = Wakeuptime_Free - Sleeptime_Free,
         Sleepduration_Free = ifelse(Sleepduration_Free < 0, 24*3600 + Sleepduration_Free, Sleepduration_Free),
         Sleepduration_Work = Wakeuptime_Work - Sleeptime_Work,
         Sleepduration_Work = ifelse(Sleepduration_Work < 0, 24*3600 + Sleepduration_Work, Sleepduration_Work),
         
         MSF = (Sleeptime_Free + 0.5*(Sleepduration_Free)), #mid sleep free days (Vetter, Roenneberg, Springr Methods Book S. Brown)
         MSF = ifelse(MSF > 24*3600, MSF-24*3600, MSF) %>% seconds_to_period(),
         
         Sleepduration_avg = (Sleepduration_Work * 5 + Sleepduration_Free * 2) / 7,
         MSF_sc = MSF %>% period_to_seconds() - 0.5*(Sleepduration_Free - Sleepduration_avg), #sleep debt-corrected MSF (Vetter Springr)
         MSF_sc = MSF_sc %>% seconds_to_period(),
         MSF_sc = round(MSF_sc %>% time_length(), 0) %>% as_hms,
         MSF_sc = round_hms(as_hms(MSF_sc), 60)) %>%
  dplyr::rename(c("subject"="Subject", "sex"="Sex")) 

info_subjects <- info_subjects_long %>% dplyr::select(subject, sex, Light_condition, age, MSF_sc)
if (!file.exists('resources/info_subjects_short.csv')){
  write.csv(info_subjects, 'resources/info_subjects_short.csv')
}


# 2. Read image files
# -------------------
if (!file.exists("data/raw_MA_files.RDS")){ 
  files <- list.files("/extra/Skin Data/raw data/", full.names = TRUE)
  m <- regexpr("(D|E)\\d+_P\\d+", files, perl = TRUE)
  names <- regmatches(files, m)
  images <- read.maimages(files = files, source = "agilent", green.only = TRUE, names = names, 
                          other.columns = "gIsWellAboveBG")
  
  #annotation of probes
  images$genes$Symbol <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "SYMBOL")
  images$genes$ENSEMBL <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "ENSEMBL")
  images$genes$EntrezID <- mapIds(hgug4112a.db, images$genes$ProbeName, keytype = "PROBEID", column = "ENTREZID")
  
  saveRDS(images, file = "data/raw_MA_files.RDS",compress = "gzip")
} else {
  images <- readRDS("data/raw_MA_files.RDS")
}


# 3. QC and normalization on raw data
# -----------------------------------
y <- limma::backgroundCorrect(images, method = "normexp") 
y <- normalizeBetweenArrays(y, method = "quantile")


# 4. Filter out controls, non-annotated probes and lowly expressed genes from raw data
# ------------------------------------------------------------------------------------
Control <- y$genes$ControlType==1L #Control <- y$genes$ControlType==1 would also work
NoID    <- is.na(y$genes$ENSEMBL)
IsExpr  <- rowSums(y$other$gIsWellAboveBG>0) >= 77 

y0 <- y[!Control & !NoID & IsExpr, ]  # Data (expressed, identified, not controls) with gene annotation
y0$genes <- y0$genes[, c("ProbeName", "Symbol", "ENSEMBL", "EntrezID")] 

yave <- avereps(y0, y0$genes[, "ENSEMBL"])  # Averaging probes mapping to the same gene
rownames(yave$E) <- yave$genes$ProbeName


# 5. PCA of raw data -> What separates first? Are there outliers?
# ---------------------------------------------------------------
metad = data.frame(feature=colnames(yave$E), row.names = colnames(yave$E))
metad$layer = ifelse(grepl("D", row.names(metad)), "dermis", "epidermis")
p <- pca(yave$E, metadata=metad, removeVar = 0.1)
pca_biplot <- biplot(p, showLoadings = FALSE, 
       labSize = 3, pointSize = 3, sizeLoadingsNames = 3, max.overlaps = 10,
       colby = 'layer', colkey = c('dermis' = '#1B9E77', 'epidermis' = '#D95F02'),
       legendPosition = 'right')
pca_pairsplot <- pairsplot(p, colby = 'layer', colkey = c('dermis' = '#1B9E77', 'epidermis' = '#D95F02'))

if (!file.exists(paste0("figures/preanalysis_PCA_biplot.pdf"))){ 
  pca_biplot %>% ggsave(paste0("figures/preanalysis_PCA_biplot.pdf"), .) 
  pca_pairsplot %>% ggsave(paste0("figures/preanalysis_PCA_pairsplot.pdf"), .) 
} 
outliers = "E32_P109"  
  

# 6. Save raw data (annotated, normalized, etc)
# ---------------------------------------------
if (!file.exists("results/rawdata.rds")){ 
  saveRDS(yave, file = "results/rawdata.rds")
} 
  

# 7. Extract sample details from column names
# -------------------------------------------
if (!file.exists("results/experiment.rds")){ 
  experiment <- data.frame(tissue = character(), time = integer(), subject = character()) %>%
    {strcapture("(\\w)(\\d+)_(\\w+)", colnames(y0$E), ., perl = TRUE)}
  saveRDS(experiment, "results/experiment.rds")
} 


# 8. Differential expression analysis
# ------------------------------------
experiment <- readRDS("results/experiment.rds") %>% # read sample details from column names
  full_join(info_subjects) %>% # we're going to correct wall time (sampling time) to internal time
  dplyr::mutate(MSF_sc_dec = lubridate::hms(MSF_sc)) %>% 
  dplyr::mutate(
    MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + lubridate::second(MSF_sc_dec) / 360),2),
    diff_to_refsubj = MSF_sc_dec - median(MSF_sc_dec),
    internal_time_ref = time - diff_to_refsubj) %>%
  mutate(internal_time = time - MSF_sc_dec)

# Remove outliers in yave
ind <- which(colnames(yave) == outliers)    
yave <- yave[, -ind] ; yave_E <- yave[,-(ind-77)]
experiment <- experiment[-ind,]
nrow(experiment) == ncol(yave) 
wts = NULL #since we remove PCA outliers, weights is set to NULL

# Prepare sample details for design matrix
tissue  <- factor(experiment$tissue)
time    <- experiment$internal_time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design and contrast matrices
design <- model.matrix(~ 0 + tissue + tissue:inphase + tissue:outphase) %>% #H0: rhythms are different across tissues
  as.data.frame() %>% dplyr::rename(c("tissueD_inphase"="tissueD:inphase", "tissueE_inphase"="tissueE:inphase",
                                      "tissueD_outphase"="tissueD:outphase", "tissueE_outphase"="tissueE:outphase")) %>% 
  as.matrix()
contr_matrix <- makeContrasts(DvsE = tissueD - tissueE, levels=design)

# duplicate Correlations
dupcor <- duplicateCorrelation(yave, design, block=subject)

# fit yave to linear model and compute contrast to check for DE
fit <- limma::lmFit(yave, design, weights = wts, block=subject, correlation=dupcor$consensus) #fit yave to the linear model
fit_contrast <- limma::contrasts.fit(fit, contr_matrix) #compute contrast
fit_contrast <- limma::eBayes(fit_contrast, trend = TRUE, robust = TRUE) #Bayes statistics of differential expression
limma::volcanoplot(fit_contrast, names=fit$genes$Symbol, main="Dermis vs. Epidermis")

#results_DE <- limma::topTable(fit_contrast, number = Inf, sort.by = "none") #differential expression analysis results
#results_DE$DE_dermis <- ifelse( results_DE$logFC > 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
#results_DE$DE_epidermis <- ifelse( results_DE$logFC < 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
#results_DE$diff_expressed <- ifelse(results_DE$DE_dermis == TRUE, "dermis",
#                                    ifelse(results_DE$DE_epidermis == TRUE, "epidermis", "ns") )
#
#volplot <- ggplot(data=results_DE, aes(x=logFC, y=-log10(adj.P.Val), col=diff_expressed, label=Symbol)) +
#  geom_point(alpha=0.5) + 
#  theme_minimal() +
#  #geom_text_repel(data = results_DE %>% filter(logFC > 0 & adj.P.Val < 0.05) %>% arrange(desc(logFC)) %>% head(25),
#  #                aes(x=logFC, y=-log10(adj.P.Val), label=Symbol),
#  #                max.overlaps=Inf, box.padding=1, size=3.5, point.padding=.5,
#  #                segment.color="grey70", color="black", parse=TRUE) +
#  #geom_text_repel(data = results_DE %>% filter(logFC < 0 & adj.P.Val < 0.05) %>% arrange(logFC) %>% head(25),
#  #                aes(x=logFC, y=-log10(adj.P.Val), label=Symbol),
#  #                max.overlaps=Inf, box.padding=1, size=3.5, point.padding=.5,
#  #                segment.color="grey70", color="black", parse=TRUE) +
#  scale_color_manual(values=c("#1B9E77", "#D95F02", "black")) +
#  geom_hline(yintercept=-log10(0.05), col="grey") + theme_custom() +
#  ggtitle("Differential expression dermis vs. epidermis") + 
#  annotate(geom='text', label=paste0('  ', dim(filter(results_DE, DE_dermis))[1], ' DE genes in D\n',
#                                     '  ', dim(filter(results_DE, DE_epidermis))[1], ' DE genes in E\n',
#                                     '  ', dim(filter(results_DE, diff_expressed=="ns"))[1], ' genes n.s.'), 
#           x=-Inf, y=Inf, hjust=0, vjust=1, color="black") + labs(x="log2 fold change")
#
#volplot %>% ggsave('figures/DE_analysis.png', ., width = 11, height = 3)

tfit <- treat(fit_contrast, lfc=1) #require logFCs to be above a minimum value of 1
dt <- decideTests(tfit)
summary(dt)

results_DE <- limma::topTable(tfit, number = Inf, sort.by = "none")
results_DE$DE_dermis <- ifelse( results_DE$logFC > 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
results_DE$DE_epidermis <- ifelse( results_DE$logFC < 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
results_DE$diff_expressed <- ifelse(results_DE$DE_dermis == TRUE, "dermis",
                                    ifelse(results_DE$DE_epidermis == TRUE, "epidermis", "ns") )
plotMD(tfit, column=1, status=dt[,1], main="Differential expression dermis vs. epidermis\n(min logFC=1)",
                 hl.col=c("#D95F02", "#1B9E77"))

top25_dermis <- results_DE %>% filter(DE_dermis==TRUE & adj.P.Val < 0.05) %>% arrange(desc(logFC)) %>% head(25)
top25_epidermis <- results_DE %>% filter(DE_epidermis==TRUE & adj.P.Val < 0.05) %>% arrange(desc(logFC)) %>% head(25)

# enrichment of DE expressed genes
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ReactomePA))

# GO
gD <- enrichGO(yave[yave$genes$Symbol %in% filter(results_DE, DE_dermis)$Symbol,]$genes$EntrezID,
               universe = yave$genes$EntrezID, ont="BP",
               'org.Hs.eg.db', pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")

gE <- enrichGO(yave[yave$genes$Symbol %in% filter(results_DE, DE_epidermis)$Symbol,]$genes$EntrezID,
               universe = yave$genes$EntrezID, ont="BP",
               'org.Hs.eg.db', pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")

g <- gE %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="epidermis") %>% as.data.frame() %>%
  rbind(gD %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="dermis") %>% as.data.frame()) %>%
  mutate(qvalue = scales::scientific(qvalue, digits = 3),
         pvalue = scales::scientific(pvalue, digits = 3),
         p.adjust = scales::scientific(p.adjust, digits = 3),
         P.DE=as.numeric(pvalue)) %>% 
  dplyr::select(-ID, -junk, -Count) %>% rename(c("Term"="Description"))

GO_DE <- ggplot(g, aes(x=-log10(P.DE), y=tidytext::reorder_within(Term, -log10(P.DE), tissue), color=tissue, size=hits)) + 
  geom_point() +  
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(2.0,5)) + 
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  guides(color = FALSE) + 
  theme_custom() + tidytext::scale_y_reordered() +
  theme(aspect.ratio=1.8, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 
GO_DE %>% ggsave('figures/GO_DE.png', ., width = 11, height = 5.5)

# Reactome
rD <- enrichPathway(gene = filter(results_DE, DE_dermis)$EntrezID, 
                    universe = yave$genes$EntrezID, minGSSize = 20, pvalueCutoff = 0.05, qvalueCutoff = 1.0)
rE <- enrichPathway(gene = filter(results_DE, DE_epidermis)$EntrezID, 
                    universe = yave$genes$EntrezID, minGSSize = 20, pvalueCutoff = 0.05, qvalueCutoff = 1.0)
df_rD <- setReadable(rD, 'org.Hs.eg.db', 'ENTREZID') %>% as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/") %>%
  mutate(hits=DE*100/N) 
df_rE <- setReadable(rE, 'org.Hs.eg.db', 'ENTREZID') %>% as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/") %>%
  mutate(hits=DE*100/N) 

Reac_DE1 <- ggplot(df_rD %>% arrange(desc(-log10(p.adjust))) %>% head(20), 
                 aes(x=-log10(p.adjust), y=reorder(Description, -log10(p.adjust)), size=DE)) + 
  geom_point( color='#1B9E77') +  
  expand_limits(x=c(0,2.0)) + 
  labs(x=bquote(~-log[10]*' adj.'~italic('p')~'value'), y="Pathway", size="Percentage\nof hits") + 
  guides(color = "none") +
  theme_custom() + tidytext::scale_y_reordered() +
  theme(aspect.ratio=2.2, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 
Reac_DE2 <- ggplot(df_rE, aes(x=-log10(p.adjust), y=reorder(Description, -log10(p.adjust)), size=DE)) + 
  geom_point( color='#D95F02') +  
  expand_limits(x=c(0,2.0)) + 
  labs(x=bquote(~-log[10]*' adj.'~italic('p')~'value'), y="Pathway", size="Percentage\nof hits") + 
  guides(color = "none") +
  theme_custom() + tidytext::scale_y_reordered() +
  theme(aspect.ratio=.3, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 
Reac_DE <- plot_grid(Reac_DE1, Reac_DE2, nrow=2, rel_heights = c(1,1))
Reac_DE %>% ggsave('figures/Reac_DE.png', ., width = 11, height = 7.5)
