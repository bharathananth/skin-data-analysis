# open Project 'skin-data-analysis.Rproj'

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hgug4112a.db))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(PCAtools))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(viridis))


dir.create("figures",showWarnings = FALSE)

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
info_subjects_long <- read.csv("data/info_subjects.csv") %>%
  dplyr::mutate(Bedtime_Work = Bedtime_Work %>% parse_hms(), Sleeptime_Work =  Sleeptime_Work %>% parse_hms(), 
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
if (!file.exists('data/info_subjects_short.csv')){
  write.csv(info_subjects, 'data/info_subjects_short.csv')
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
Control <- y$genes$ControlType==1L 
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
       labSize = 3, pointSize = 1, sizeLoadingsNames = 3, max.overlaps = 10,
       colby = 'layer', colkey = c('dermis' = '#1B9E77', 'epidermis' = '#D95F02'),
       legendPosition = 'right')
pca_pairsplot <- pairsplot(p, colby = 'layer', colkey = c('dermis' = '#1B9E77', 'epidermis' = '#D95F02'))
suppfig1A = pca_biplot + theme_custom() + theme(aspect.ratio=.7) +
  guides(color = guide_legend(override.aes = list(size=3)))

outliers = "E32_P109"  
  

# 6. Save raw data (annotated, normalized, etc)
# ---------------------------------------------
if (!file.exists("data/rawdata.rds")){ 
  saveRDS(yave, file = "data/rawdata.rds")
} 
  

# 7. Extract sample details from column names
# -------------------------------------------
if (!file.exists("data/experiment.rds")){ 
  experiment <- data.frame(layer = character(), time = integer(), subject = character()) %>%
    {strcapture("(\\w)(\\d+)_(\\w+)", colnames(y0$E), ., perl = TRUE)}
  saveRDS(experiment, "data/experiment.rds")
} 


# 8. Differential expression analysis
# ------------------------------------
experiment <- readRDS("data/experiment.rds") %>% # read sample details from column names
  full_join(info_subjects) %>% # we're going to correct wall time (sampling time) to internal time
  dplyr::mutate(MSF_sc_dec = lubridate::hms(MSF_sc)) %>% 
  dplyr::mutate(
    MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + lubridate::second(MSF_sc_dec) / 360),2),
    diff_to_refsubj = MSF_sc_dec - median(MSF_sc_dec),
    internal_time_ref = time - diff_to_refsubj) %>%
  dplyr::mutate(internal_time = time - MSF_sc_dec)

# Remove outliers in yave
ind <- which(colnames(yave) == outliers)    
yave <- yave[, -ind] ; yave_E <- yave[,-(ind-77)]
experiment <- experiment[-ind,]
nrow(experiment) == ncol(yave) 
wts = NULL #since we remove PCA outliers, weights is set to NULL

# Prepare sample details for design matrix
layer  <- factor(experiment$layer)
time    <- experiment$internal_time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design and contrast matrices
design <- model.matrix(~ 0 + layer + layer:inphase + layer:outphase) %>% #H0: rhythms are different across layers
  as.data.frame() %>% dplyr::rename(c("layerD_inphase"="layerD:inphase", "layerE_inphase"="layerE:inphase",
                                      "layerD_outphase"="layerD:outphase", "layerE_outphase"="layerE:outphase")) %>% 
  as.matrix()
contr_matrix <- makeContrasts(DvsE = layerD - layerE, levels=design)

# duplicate Correlations
dupcor <- duplicateCorrelation(yave, design, block=subject)

# fit yave to linear model and compute contrast to check for DE
fit <- limma::lmFit(yave, design, weights = wts, block=subject, correlation=dupcor$consensus) #fit yave to the linear model
fit_contrast <- limma::contrasts.fit(fit, contr_matrix) #compute contrast
fit_contrast <- limma::eBayes(fit_contrast, trend = TRUE, robust = TRUE) #Bayes statistics of differential expression
limma::volcanoplot(fit_contrast, names=fit$genes$Symbol, main="Dermis vs. Epidermis")

logfc_cutoff = 1
tfit <- treat(fit_contrast, lfc=logfc_cutoff) #require logFCs to be above a minimum value of 1
dt <- decideTests(tfit)
summary(dt)

results_DE <- limma::topTable(tfit, number = Inf, sort.by = "none")
results_DE$DE_dermis <- ifelse( results_DE$logFC > 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
results_DE$DE_epidermis <- ifelse( results_DE$logFC < 0 & results_DE$adj.P.Val < 0.05, TRUE, FALSE )
results_DE$diff_expressed <- ifelse(results_DE$DE_dermis == TRUE, "dermis",
                                    ifelse(results_DE$DE_epidermis == TRUE, "epidermis", "ns") )

tfit %<>% as.data.frame() %>% 
  mutate(DE_binary=ifelse(abs(coefficients) > logfc_cutoff, TRUE, FALSE),
         DE=ifelse(coefficients > logfc_cutoff, "dermis", ifelse(coefficients < -logfc_cutoff, "epidermis", "ns")))
suppfig1B <- ggplot(tfit) + 
  geom_point(data = tfit %>% filter(!DE_binary), aes(x=Amean, y=coefficients), size=.1, color="black") +
  geom_point(data = tfit %>% filter(DE_binary), aes(x=Amean, y=coefficients, color=DE), size=.7) +
  labs(x="Average log-expression", y="log-fold change") +
  theme_custom() + theme(aspect.ratio=.7) +
  guides(color = guide_legend(override.aes = list(size=3)))

# Supplementary Table 2: Lists of top differentially expressed genes between dermis and epidermis
rank_dermis <- results_DE %>% dplyr::filter(DE_dermis==TRUE & adj.P.Val<0.05) %>% arrange(desc(logFC)) %>% 
  dplyr::select(-ENSEMBL, -EntrezID, -t, -P.Value, -DE_dermis, -DE_epidermis)
rank_epidermis <- results_DE %>% dplyr::filter(DE_epidermis==TRUE & adj.P.Val<0.05) %>% arrange(desc(abs(logFC))) %>% 
  dplyr::select(-ENSEMBL, -EntrezID, -t, -P.Value, -DE_dermis, -DE_epidermis)
sheets <- list("differential expr. dermis" = format.data.frame(rank_dermis, digits=3), 
               "differential expr. epidermis" = format.data.frame(rank_epidermis, digits=3))
if (!file.exists("figures/supp_table2.xlsx")){
  openxlsx::write.xlsx(sheets, file = "figures/supp_table2.xlsx")
}

# Enrichment of DE expressed genes (Reactome pathways)
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

df_rDE <- df_rD %>% mutate(layer="dermis") %>% arrange(desc(-log10(p.adjust))) %>% head(5) %>% 
  full_join(df_rE %>% mutate(layer="epidermis") %>% arrange(desc(-log10(p.adjust))) %>% head(5))

suppfig1C <- ggplot(df_rDE, aes(x=-log10(p.adjust), y=reorder(Description, -log10(p.adjust)), size=DE, color=layer)) + 
  geom_point() +  
  facet_wrap(~layer, scales='free_y') + expand_limits(x=c(2.0,5)) + 
  labs(x=bquote(~-log[10]*' adj.'~italic('p')~'value'), y="Pathway", size="Percentage\nof hits") + 
  guides(color = "none") + 
  theme_custom() + tidytext::scale_y_reordered() +
  theme(aspect.ratio=.5, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 

sfig1_1 <- plot_grid(suppfig1A, NULL, suppfig1B, ncol=3, rel_widths=c(1,0.1,1), labels=c("A", "", "B"))
sfig1_2 <- plot_grid(suppfig1C, labels=c("C"))
sfig1 = plot_grid(sfig1_1, NULL, sfig1_2, nrow=3, rel_heights = c(1,0.1,.5))
sfig1 %>% ggsave('figures/suppfig1.pdf', ., width = 11, height = 5.5)


# 9. Melatonin and cortisol profiles
# ----------------------------------
melatonin <- read.csv("./data/hormones/melatonin.csv") %>%
  tidyr::gather(subject, value, -Zeit) %>%
  dplyr::mutate(mutate(across('subject', str_replace, 'X', 'P'))) %>%
  dplyr::rename(c("Time"="Zeit")) %>% 
  dplyr::filter(subject %in% experiment$subject) %>%
  dplyr::mutate(Time = rep(c(8,12,16,20,24,28,32), 11),
                hormone = "melatonin")
cortisol <- read.csv("./data/hormones/cortisol.csv") %>%
  tidyr::gather(subject, value, -Zeit) %>%
  dplyr::mutate(mutate(across('subject', str_replace, 'X', 'P'))) %>%
  dplyr::rename(c("Time"="Zeit")) %>% 
  dplyr::filter(subject %in% experiment$subject) %>%
  dplyr::mutate(Time = rep(c(8,12,16,20,24,28,32), 11),
                hormone = "cortisol")

suppfig5A <- ggplot(melatonin, aes(x=Time, y=value)) + geom_line() + geom_point(size=1) +
  facet_wrap(~subject, ncol=3, scales="free") + 
  theme_custom() + 
  scale_x_continuous(breaks=c(8,16,24,32), labels=c("8:00", "16:00", "00:00", "8:00")) + 
  scale_y_continuous(limits = c(0,50)) +
  labs(x="wall time (h)", y="Melatonin (pg/mL)") +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))
suppfig5B <- ggplot(cortisol, aes(x=Time, y=value)) + geom_line() + geom_point(size=1) +
  facet_wrap(~subject, ncol=3, scales="free") + 
  theme_custom() + 
  scale_x_continuous(breaks=c(8,16,24,32), labels=c("8:00", "16:00", "00:00", "8:00")) + 
  scale_y_continuous(limits = c(0,1.35), breaks=c(0, 0.3, 0.6, 0.9, 1.2)) +
  labs(x="wall time (h)", y="Cortisol (µg/dL)") +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))

suppfig5 <- plot_grid(suppfig5A, NULL, suppfig5B, ncol=3, labels=c("A","", "B"), rel_widths=c(1,0.05,1))
suppfig5 %>% ggsave('./figures/suppfig5.pdf', ., width = 11, height = 8.5)

##########
##########