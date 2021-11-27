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

setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis")

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

# Find expressed genes in each skin layer separately
yD <- y[,c(1:77)]
yE <- y[,c(78:154)]

# Remove lowly expressed genes 
IsExpr  <- rowSums(yD$other$gIsWellAboveBG>0) >= 39
y0D <- yD[!Control & !NoID & IsExpr, ]  # Data (expressed, identified, not controls) with gene annotation
y0D$genes <- y0D$genes[, c("ProbeName", "Symbol", "ENSEMBL", "EntrezID")] 
yaveD <- avereps(y0D, y0D$genes[, "ENSEMBL"]) 

IsExpr  <- rowSums(yE$other$gIsWellAboveBG>0) >= 39
y0E <- yE[!Control & !NoID & IsExpr, ]  # Data (expressed, identified, not controls) with gene annotation
y0E$genes <- y0E$genes[, c("ProbeName", "Symbol", "ENSEMBL", "EntrezID")] 
yaveE <- avereps(y0E, y0E$genes[, "ENSEMBL"]) 

yaveDE <- list(expr_D = yaveD$genes$ENSEMBL, 
               expr_E = yaveE$genes$ENSEMBL,
               expr_whole = yave$genes$ENSEMBL)

# Check overlap between expressed genes in D vs E or whole skin
fig0_1 <- ggvenn(yaveDE[c(1,2)], fill_color = c("#1B9E77", "#D95F02"), text_size = 3,
                 stroke_size = 0.25, set_name_size = 4) + ggtitle("overlap expressed genes dermis vs epidermis")
fig0_2 <- ggvenn(yaveDE[c(1,3)], fill_color = c("#1B9E77", "grey"), text_size = 3,
                 stroke_size = 0.25, set_name_size = 4) + ggtitle("overlap expressed genes dermis vs whole")
fig0_3 <- ggvenn(yaveDE[c(2,3)], fill_color = c("#D95F02", "grey"), text_size = 3,
                 stroke_size = 0.25, set_name_size = 4) + ggtitle("overlap expressed genes epidermis vs whole")
fig0 <- plot_grid(fig0_1, fig0_2, fig0_3, align='v', nrow=3)


# 5. PCA of raw data -> What separates first? Are there outliers?
# ---------------------------------------------------------------
scale_PCA = TRUE
for_filename <- ifelse(scale_PCA==TRUE, "T", "F")

yave.pca <- prcomp(yave$E, center = TRUE, scale. = scale_PCA)
#summary(yave.pca)

df_out_r <- as.data.frame(yave.pca$rotation)
df_out_r$tissue <- ifelse(grepl("D", row.names(df_out_r)), "dermis", "epidermis")
df_out_r$feature <- row.names(df_out_r)

plot_PCA_all <- ggplot(df_out_r, aes(x=PC1, y=PC2, label=feature, color=tissue)) +
  geom_point() + theme_bw() + geom_text(size=3) + ggtitle("PCA all data")

if (!file.exists(paste0("figures/final/preanalysis_PCA_scale=", for_filename,".pdf"))){ 
  for_filename <- ifelse(scale_PCA==TRUE, "T", "F")
  plot_PCA_all %>% ggsave(paste0("figures/preanalysis_PCA_scale=", for_filename,".pdf"), .) 
} 
outliers = "E32_P109"  
  

# 6. Save raw data (annotated, normalized, etc)
# ---------------------------------------------
if (!file.exists("visualize/data/rawdata.rds")){ 
  saveRDS(yave, file = "visualize/data/rawdata.rds")
} 
if (!file.exists("visualize/data/rawdata_dermis.rds")){ 
  saveRDS(yaveD, file = "visualize/data/rawdata_dermis.rds")
  saveRDS(yaveE, file = "visualize/data/rawdata_epidermis.rds")
} 
  


# 7. Extract sample details from column names
# -------------------------------------------
if (!file.exists("visualize/data/experiment.rds")){ 
  experiment <- data.frame(tissue = character(), time = integer(), subject = character()) %>%
    {strcapture("(\\w)(\\d+)_(\\w+)", colnames(y0$E), ., perl = TRUE)}
  saveRDS(experiment, "visualize/data/experiment.rds")
} 