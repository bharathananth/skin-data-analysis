library(limma)
library(magrittr)
library(hgug4112a.db)
library(ggrepel)
library(ggforce)
library(statmod)
library(GO.db)
library(tibble)
library(massiR)
library(tidyr)
library(dplyr) 
library(tidyr)
library(ggplot2)
library(cowplot)
library(hms)
#library(round_hms)
library(lubridate)
library(gridExtra)

#--------------------------------
#--------------------------------

# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 

#--------------------------------
#--------------------------------

# 2. READ FILES
# -------------
# Read info of subjects
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
info_subjects <- info_subjects_long %>% select(subject, sex, Light_condition, age, MSF_sc)

# Read image files
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

# QC and normalization on raw data
if (!file.exists("visualize/data/rawdata.rds")){
  y <- limma::backgroundCorrect(images, method = "normexp") 
  y <- normalizeBetweenArrays(y, method = "quantile")
}

# Filter out controls, non-annotated probes and lowly expressed genes
if (!file.exists("visualize/data/rawdata.rds")){ 
  Control <- y$genes$ControlType==1L #Control <- y$genes$ControlType==1 would also work
  NoID    <- is.na(y$genes$ENSEMBL)
  IsExpr  <- rowSums(y$other$gIsWellAboveBG>0) >= 77 
  
  y0 <- y[!Control & !NoID & IsExpr, ]  # Data (expressed, identified, not controls) with gene annotation
  y0$genes <- y0$genes[, c("ProbeName", "Symbol", "ENSEMBL", "EntrezID")] 
  
  yave <- avereps(y0, y0$genes[, "ENSEMBL"])  # Averaging probes mapping to the same gene
  rownames(yave$E) <- yave$genes$ProbeName
}

# Save raw data (annotated, normalized, etc)
if (!file.exists("visualize/data/rawdata.rds")){ 
  saveRDS(yave, file = "visualize/data/rawdata.rds")
} else {
  yave <- readRDS("visualize/data/rawdata.rds")
}

# Extract sample details from column names
if (!file.exists("visualize/data/rawdata.rds")){ 
  experiment <- data.frame(tissue = character(), time = integer(), subject = character()) %>%
    {strcapture("(\\w)(\\d+)_(\\w+)", colnames(y0$E), ., perl = TRUE)}
  
  tissue  <- factor(experiment$tissue)
  time    <- experiment$time
  subject <- factor(experiment$subject)
  saveRDS(experiment, "visualize/data/experiment.rds")
} else{
  experiment <- readRDS("visualize/data/experiment.rds")
  experiment %<>% full_join(info_subjects)
  
  tissue  <- factor(experiment$tissue)
  time    <- experiment$time
  subject <- factor(experiment$subject)
  sex <- factor(experiment$sex)
  age <- experiment$age 
}


  
#--------------------------------
#--------------------------------

# 3. DO LINEAR MODEL FITS AND SAVE RESULTS OF FITS
inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

if (!file.exists("visualize/data/results_fig1.rds")){ #Fig1: assuming rhythms are the same across subjects
  # design matrix
  design <- model.matrix(~ 0 + subject + tissue + tissue:inphase + tissue:outphase)
  
  # weights + fits
  wts <- limma::arrayWeights(yave, model.matrix(~tissue + subject)) # model.matrix: simpler model used as suggested in userguide. (avoid?)
  fit <- limma::lmFit(yave, design)#, weights = wts) 
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE) 
  
  # save results of fits
  rhy_indices <- which(grepl("phase",colnames(design)))
  results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
    set_colnames(gsub("\\.","_", colnames(.))) 
  saveRDS(results, file = "visualize/data/results_fig1.rds")
}

#--------------------------------
#--------------------------------

# 4. READ RESULTS AND PLAY WITH THEM
# ----------------------------------
results <- readRDS("visualize/data/results_fig1.rds")
results %<>% 
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)

#--------------------------------
#--------------------------------

# 5. FIGURES
# ----------

# Fig1B: Demographics -- how does our date look like? ~ Chronotype, age, sex...
lims <- as.POSIXct(strptime(c("1970-01-01 02:30:00","1970-01-01 07:15:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
fig1B <- ggplot(experiment) +
  geom_point(aes(as.POSIXct(x = MSF_sc, format = "%H:%M"), age, shape=sex), size=4) + 
  scale_y_continuous(breaks=seq(18,32,2)) +
  scale_x_datetime(date_breaks = "1 hours", date_labels = "%H:%M", limits=lims) +
  theme_bw() + theme(aspect.ratio=1) + 
  xlab("Mid sleep time") + ylab("Age") 

# Fig1C: How does number of rhy genes in D and E change with FDR
results_dermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val, A_D) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_epidermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val, A_E) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_org <- rbind(results_epidermis, results_dermis)
no_dermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="dermis"))[1] # #of genes passing Amp threshold in dermis
no_epidermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="epidermis"))[1] #same but in epidermis

fig1C <- ggplot(data = results_org %>% mutate(len = ifelse(tissue=="dermis", no_dermis, no_epidermis)) %>% 
                  filter(amp_value > amp_cutoff)) + 
  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
  stat_ecdf(aes(x = ifelse(amp_value > amp_cutoff, adj_P_Val), len=len, y = ..y..*len, color=tissue), geom = "step") +  
  scale_x_continuous(breaks = seq(0.00, 0.05, by=0.01)) +
  #facet_wrap(~tissue) + 
  coord_cartesian(xlim = c(-0.0002, 0.053), ylim = c(-0.0002,1500), expand=FALSE) + theme_bw() +
  xlab("False discovery rate") + ylab("number of rhythmic genes") 

#fig1C_1 <- ggplot(data = dplyr::filter(results, A_D > amp_cutoff) %>% #pmax(A_D, A_E)
#                             dplyr::mutate(len = length(adj_P_Val))) + 
#  stat_ecdf(aes(x = adj_P_Val, len=len, y = ..y..*len), geom = "step") +  
#  coord_cartesian(xlim = c(0.001, 0.05), ylim = c(0,1300)) + theme_bw() +
#  xlab("False discovery rate") + ylab("number of rhythmic genes in dermis") 
#fig1C_2 <- ggplot(data = dplyr::filter(results, A_E > amp_cutoff) %>% 
#                    dplyr::mutate(len = length(adj_P_Val))) + 
#  stat_ecdf(aes(x = adj_P_Val, len=len, y = ..y..*len), geom = "step") +  
#  coord_cartesian(xlim = c(0.001, 0.05), ylim = c(0,1300)) + theme_bw() +
#  xlab("False discovery rate") + ylab("number of rhythmic genes in epidermis") 

results_org$tissue <- factor(results_org$tissue, levels = c("epidermis", "dermis"))
suppfig1 <- ggplot(data = results_org %>% filter(amp_value > amp_cutoff & adj_P_Val < fdr_cutoff)) + 
  geom_vline(aes(xintercept=amp_cutoff), color="black", linetype="dashed") +
  #geom_histogram(aes(amp_value, color=tissue), position="identity", fill="white", bins=50, alpha=0.5) + 
  geom_histogram(aes(amp_value, fill=tissue), color="black", bins=50, size=0.2) + 
  theme_bw() + theme(aspect.ratio=1) + xlab("Amplitude") + ylab("number of genes\n(passing amplitude and fdr cutoff)")


# Fig 1D: What happens with the core clock genes in D and E?
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
results_amp <- results %>% 
  filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(Symbol, adj_P_Val, A_D, A_E) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_phase <- results %>% 
  filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(Symbol, adj_P_Val, phaseD, phaseE) %>%
  gather(tissue, phase_value, -adj_P_Val, -Symbol) %>%
  mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))
results_org <- full_join(results_amp, results_phase) %>% 
  filter(Symbol %in% clock_genes)

fig1D <- ggplot(results_org, aes(x=phase_value, y=amp_value, color=tissue)) +
  geom_point() +
  coord_polar() + theme_bw() +
  scale_x_continuous(limits=c(-12, +12), breaks = c(0, 6, 12, -6) ) +
  facet_wrap(~Symbol) +
  geom_segment(aes(x=phase_value, y=0, xend=phase_value, yend=amp_value)) +
  xlab("Phase (h)") + ylab("Amplitude")
  # Amplitude of NPAS2 in dermis is < than amp_cutoff

#fig1D_2 <- ggplot(results_org, aes(x=phase_value, y=amp_value, color=tissue, label=Symbol)) +
#  geom_point() +
#  geom_text_repel() +
#  coord_polar() + theme_bw() +
#  scale_x_continuous(limits=c(-12, +12), breaks = c(0, 6, 12, -6) ) +
#  geom_segment(aes(x=phase_value, y=0, xend=phase_value, yend=amp_value, color=tissue)) 

#--------------------------------
#--------------------------------

# ARRANGE PLOTS IN GRID
# ---------------------
right_grid <- plot_grid(fig1D, ncol=1, labels="D")
left_grid <- plot_grid(fig1B, fig1C, ncol= 1, nrow=2, labels=c("B", "C", ""), rel_heights=c(1,0.5,0.01))

fig1 <- plot_grid(left_grid, NULL, right_grid, align='v', nrow=1, rel_widths = c(1.5,0.1,2.25))

# SAVE FIGURES
# ------------
if (!file.exists("figures/fig1.pdf")){ 
  fig1 %>% ggsave("figures/fig1.pdf", ., width = 11,height = 8.5) 
}
if (!file.exists("figures/suppfig1.pdf")){ 
  suppfig1 %>% ggsave("figures/suppfig1.pdf", ., width=5, height=5) 
}



