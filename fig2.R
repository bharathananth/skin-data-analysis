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
library(tidyverse)
library(gridExtra)
library(mgsub)
library(corrplot)
library(ggvenn)
library(GGally)
library(Hmisc)
suppressPackageStartupMessages(library(clusterProfiler))
library(msigdbr)
setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis")
scale_colour_discrete <- function(...) {
  scale_colour_brewer(..., palette="Dark2")
}

m_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
m_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

#--------------------------------
#--------------------------------

# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 
removal_outliers <- "PCA" #PCA or weights
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R
save_fitdata <- FALSE

#substract_mean_spearman <- TRUE

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
  saveRDS(experiment, "visualize/data/experiment.rds")
} else{
  experiment <- readRDS("visualize/data/experiment.rds")
  experiment %<>% full_join(info_subjects)
}
experiment %<>%
  dplyr::mutate(MSF_sc_dec = hms(MSF_sc)) %>%
  dplyr::mutate(#option1
                MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + second(MSF_sc_dec) / 360),2),
                internal_time1 = time + MSF_sc_dec,
                #option2
                diff_to_refsubj = MSF_sc_dec - median(MSF_sc_dec),
                internal_time = time + diff_to_refsubj)

# Remove outliers
if (removal_outliers == "PCA"){ #option 1
  ind <- which(colnames(yave) == PCA_outliers)    
  yave <- yave[, -ind] 
  experiment <- experiment[-ind,]
  nrow(experiment) == ncol(yave)   #input-check
  dim(yave)
  #since we are removing PCA-outliers, weights from lmFit (performed later) are set to NULL
  wts <- NULL 
} else if(removal_outliers == "weights"){ #option2
  dim(yave)
  wts <- limma::arrayWeights(yave, model.matrix(~experiment$tissue + experiment$subject))  #model.mtx: simpler model used as suggested in userguide (avoid?
} else{
  print("`removal_outliers`` must be either 'PCA' or 'weights' -> check spelling")
}

#--------------------------------
#--------------------------------

# 3. DO LINEAR MODEL FITS AND SAVE RESULTS OF FITS

# Prepare sample details for design matrix
tissue  <- factor(experiment$tissue)
time    <- experiment$internal_time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design matrix
design <- model.matrix(~ 0 + subject + tissue + tissue:inphase + tissue:outphase) #H0: rhythms are different across tissues

# weights + fits
fit <- limma::lmFit(yave, design, weights = wts) #!!!!
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# save results of fits
if (save_fitdata == TRUE){
  rhy_indices <- which(grepl("phase",colnames(design)))
  results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
    set_colnames(gsub("\\.","_", colnames(.))) 
  saveRDS(results, file = "visualize/data/results_fig2_internaltime.rds")
}  

#--------------------------------
#--------------------------------

# 4. READ RESULTS AND PLAY WITH THEM
# ----------------------------------
results <- readRDS("visualize/data/results_fig2_internaltime.rds")
results %<>% 
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)
# phase = 0 means that the respective gene peaks at 8AM (time of first sampling)

#--------------------------------
#--------------------------------

# 5. FIGURES
# ----------

# Fig2C: How does number of rhy genes in D and E change with FDR
results_dermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val, A_D, ProbeName) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_epidermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val, A_E, ProbeName) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_org <- rbind(results_epidermis, results_dermis)
no_dermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="dermis"))[1] # #of genes passing Amp threshold in dermis
no_epidermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="epidermis"))[1] #same but in epidermis

fig2C <- ggplot(data = results_org %>% mutate(len = ifelse(tissue=="dermis", no_dermis, no_epidermis)) %>% 
                  filter(amp_value > amp_cutoff)) + 
  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
  stat_ecdf(aes(x = ifelse(amp_value > amp_cutoff, adj_P_Val), len=len, y = ..y..*len, color=tissue), geom = "step") +  
  scale_x_continuous(breaks = seq(0.00, 0.05, by=0.01)) + 
  #facet_wrap(~tissue) + 
  coord_cartesian(xlim = c(-0.0002, 0.053), ylim = c(-0.0002,1500), expand=FALSE)  +
  theme_bw() + theme(aspect.ratio=1) + 
  xlab("False discovery rate\n") + ylab("# of rhythmic genes") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position="top")


# Fig 1D: What happens with the core clock genes in D and E?
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
results_amp <- results %>% 
  #filter(adj_P_Val < fdr_cutoff) %>%
  #filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(Symbol, adj_P_Val, A_D, A_E, AveExpr, ProbeName) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol, -AveExpr, -ProbeName) %>%
  filter(amp_value > amp_cutoff) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_phase <- results %>% 
  #filter(adj_P_Val < fdr_cutoff) %>%
  filter(pmax(A_E, A_D) > amp_cutoff) %>%
  #filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(Symbol, adj_P_Val, phaseD, phaseE, ProbeName) %>%
  gather(tissue, phase_value, -adj_P_Val, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))
results_passAmpcutoff <- inner_join(results_amp, results_phase) #full_join if we want all clock genes, independently of amp > cutoff
  # results_passAmpcutoff contains genes with amp>cutoff, but no filtering with fdr has been done yet

lims <- as.POSIXct(strptime(c("1970-01-01 00:00:00","1970-01-01 23:59:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
exps <- as.POSIXct(strptime(c("1969-12-31 23:50:00","1970-01-02 00:10:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
fig2D <- ggplot(results_passAmpcutoff %>% filter(Symbol %in% clock_genes & adj_P_Val < fdr_cutoff) %>%
                  mutate(phase_clock1 = phase_value + 8) %>% 
                  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
                  mutate(phase_clock = format(.POSIXct(3600*phase_clock1, "UTC"), "%Y-%m-%d %H:%M:%S")) %>%
                  mutate(time_of_day=hms::hms(second(phase_clock),minute(phase_clock),hour(phase_clock))) %>%
                  mutate(Symbol_it = paste0("italic('", Symbol, "')")),
                aes(x=phase_clock1, y=amp_value)) +
                #aes(x=as.POSIXct(x = time_of_day, format = "%H:%M"), y=amp_value)) +
  geom_point(aes(color=tissue)) +
  coord_polar(start = 0, direction=1) + theme_bw() +
  #scale_x_datetime(date_breaks = "6 hours", date_labels = "%H:%M", limits = lims) +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 6, 12, 18), labels = c(0, 6, 12, 18)) +
  scale_y_continuous(breaks = seq(0.2,0.8,0.2)) +
  facet_wrap(~tissue, ncol=1, nrow=2) +
  #geom_segment(aes(x=as.POSIXct(x = time_of_day, format = "%H:%M"), y=0, 
  #                 xend=as.POSIXct(x = time_of_day, format = "%H:%M"), yend=amp_value, color=tissue)) +
  geom_segment(aes(x=phase_clock1, y=0, 
                   xend=phase_clock1, yend=amp_value, color=tissue)) +
  geom_text_repel(aes(label=Symbol_it), max.overlaps=Inf, box.padding=1, size=2.5, point.padding=.5,
                  segment.color="grey70", color="grey50", parse=TRUE) +
  xlab("Phase around clock") + ylab("Amplitude") + guides(color = FALSE) +
  theme(axis.line.y = element_line(colour = "black"),#, inherit.blank = TRUE),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        panel.spacing.y = unit(1.5, "lines"),
        plot.margin=unit(c(1,2,1.5,2),"cm"))#,panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 
  # Amplitude of NPAS2 in dermis is < than amp_cutoff


# Fig 1E: Heat map of circadian rhythmic genes: z scores, acropase-ordered
rhy_results <- results_passAmpcutoff %>% filter(adj_P_Val < fdr_cutoff)

#toplot <- yave$E %>% transform(ProbeName = yave$genes$ProbeName, 
#                                  Symbol = yave$genes$Symbol) %>% as_tibble() %>%
#  tidyr::gather(junk, value, -ProbeName, -Symbol) %>%
#  tidyr::separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
#  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
#  inner_join(rhy_results %>% mutate(tissue=ifelse(tissue=="dermis", "D", "E")) ) %>%
#  #dplyr::filter(Symbol %in% rhy_results$Symbol) %>%
#  dplyr::group_by(ProbeName, Symbol, tissue, time) %>%
#  dplyr::summarise(value = mean(value)) %>% 
#  dplyr::mutate(z.score = (value - mean(value)) / sd(value))
#
toplot <- yave$E %>% transform(ProbeName = yave$genes$ProbeName, 
                     Symbol = yave$genes$Symbol) %>% as_tibble() %>%
  tidyr::gather(junk, value, -ProbeName, -Symbol) %>%
  tidyr::separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
  inner_join(rhy_results %>% mutate(tissue=ifelse(tissue=="dermis", "D", "E")) ) %>%
  group_by(ProbeName, Symbol, tissue, subject) %>%
  mutate(z.score=(value-mean(value)) / sd(value)) %>% as.data.frame() %>%
  dplyr::group_by(time, Symbol, tissue) %>%
  summarise(z.score=mean(z.score))

toplot %<>% as.data.frame() %>% mutate(tissue=ifelse(tissue=="D", "dermis", "epidermis")) %>%
  full_join(rhy_results %>% select(Symbol, phase_value, tissue)) %>% 
  rename(c("phase"="phase_value")) %>% arrange(phase) 
toplot$Symbol_ord <- factor(toplot$Symbol, levels = rev(unique(toplot$Symbol)), ordered=TRUE)

fig2E <- ggplot(toplot %>% arrange(phase), 
                aes(x=time, y=Symbol_ord)) +
  geom_tile(aes(fill=z.score)) + #ylab("time (h)") +
  labs(x="time (h)", y="", fill=expression(italic('z')~'score')) + 
  guides(color = FALSE) +
  facet_wrap(~tissue, scales="free_y", ncol=1, nrow=2) + theme_bw() +
  #ylab("") + #scale_fill_gradient2() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        panel.spacing.y = unit(1.5, "lines")) + 
  scale_fill_distiller(palette = "RdBu", limits=c(-1.75,1.75), breaks = c(1.5, 0, -1.5)) + 
  scale_x_continuous(breaks=time %>% unique)  + theme(aspect.ratio=1) 

## Fig 1F: PSEA in D and E
## in a terminal, execute the .jar PSEA file: > java -jar PSEA-master/PSEA1.1_VectorGraphics.jar
#fdr_cutoff_PSEA <- 0.1
#if (!file.exists("visualize/data/phases_fig1_D.csv")){ 
#  results_passAmpcutoff %>% filter(tissue=="dermis" & adj_P_Val < fdr_cutoff_PSEA) %>% 
#    select(Symbol, phase_value) %>% mutate(phase_value=round(phase_value, 0) %>% as.numeric()) %>%
#    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_D.txt", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
#}
#
#if (!file.exists("visualize/data/phases_fig1_E.csv")){ 
#  results_passAmpcutoff %>% filter(tissue=="epidermis" & adj_P_Val < fdr_cutoff_PSEA) %>% 
#    select(Symbol, phase_value) %>% mutate(phase_value=round(phase_value, 0) %>% as.numeric()) %>%
#    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_E.txt", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
#}
## Significant enriched time-dependent pathways of circadian genes identified in dermis and epidermis determined with the PSEA
## Java app (-12, +12, 5, 10000, 0.05). A less significant cutoff (fdr<0.1) was used to select circadian genes for performing 
## PSEA analysis (q<0.05). Gene sets were downloaded from the Molecular Signatures database (MSigDB) 
## all C2 (KEGG gene sets) or C5 (GO:BP). Sets containing fewer than five circadian transcripts were excluded from the analysis. 
## The Kuiper test was used to identify circadian gene sets by comparing the  acrophases of all circadian transcripts 
## (rounded to the full hour). Domain from -12 to +12 since that is our range of phases.
## READ RUKEIA PAPER & WU2020 SUPP S1C&D
## PSEA results saved in /visualize/data/PSEA/dermis or /epidermis
#
## PSEA KEGG
#q_cutoff <- 0.25# 0.25 for KEGG, 0.05 for GO:BP
#PSEA_string <- "KEGG" #KEGG or GOBP
#
#PSEA_d_K <- read.csv(paste0("visualize/data/PSEA/dermis/results_", PSEA_string, ".txt"), sep='\t') %>% 
#  filter(str_detect(Set.ID, PSEA_string) & Set.N >= 5)
#PSEA_d_K <- PSEA_d_K[PSEA_d_K[,4] <  q_cutoff, ] %>% mutate(tissue = "dermis")
#PSEA_e_K <- read.csv(paste0("visualize/data/PSEA/epidermis/results_", PSEA_string, ".txt"), sep='\t')%>% 
#  filter(str_detect(Set.ID, PSEA_string) & Set.N >= 5)
#PSEA_e_K <- PSEA_e_K[PSEA_e_K[,4] <  q_cutoff, ] %>% mutate(tissue = "epidermis")
#PSEA_K <- full_join(PSEA_d_K %>% select(Set.ID, Set.N, Vector.average.value, tissue), 
#                    PSEA_e_K %>% select(Set.ID, Set.N, Vector.average.value, tissue)) %>%
#  mutate(term = gsub("_", " ", Set.ID)) %>% mutate(term=gsub("KEGG", "", term)) %>% mutate(term = term %>% tolower()) %>%
#  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
#                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))
#
#fig1F_KEGG <- ggplot(PSEA_K %>%
#                       mutate(phase_clock1 = Vector.average.value + 8) %>% 
#                       mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)), 
#                     aes(x=phase_clock1, y=term)) +
#  geom_point(aes(size=Set.N, color=tissue)) + 
#  facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
#  #facet_wrap(~tissue, scales="free_y", space="free_y", ncol=1, nrow=2) +
#  theme_bw() + xlab("phase (h)\n") + ylab(ifelse(PSEA_string=="GOBP", "GO:BP term", "KEGG pathway")) + 
#  guides(color = FALSE) + labs(size='# of genes') + # theme(aspect.ratio=1) + 
#  theme(axis.line = element_line(colour = "black"),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        strip.background = element_blank(),
#        #strip.text = element_text(face="bold"))
#        strip.text = element_blank(),
#        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
#  ggtitle("PSEA analysis\n") + expand_limits(y=0) + 
#  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
#  scale_size(breaks=c(5,10,15))
#
#
## PSEA GOBP
#q_cutoff <- 0.05# 0.25 for KEGG, 0.05 for GO:BP
#PSEA_string <- "GOBP" #KEGG or GOBP
#PSEA_d_GO <- read.csv(paste0("visualize/data/PSEA/dermis/results_", PSEA_string, ".txt"), sep='\t') %>% 
#  filter(str_detect(Set.ID, PSEA_string) & Set.N >= 5)
#PSEA_d_GO <- PSEA_d_GO[PSEA_d_GO[,4] <  q_cutoff, ] %>% mutate(tissue = "dermis")
#PSEA_e_GO <- read.csv(paste0("visualize/data/PSEA/epidermis/results_", PSEA_string, ".txt"), sep='\t')%>% 
#  filter(str_detect(Set.ID, PSEA_string) & Set.N >= 5)
#PSEA_e_GO <- PSEA_e_GO[PSEA_e_GO[,4] <  q_cutoff, ] %>% mutate(tissue = "epidermis")
#
#PSEA_GO <- full_join(PSEA_d_GO %>% select(Set.ID, Set.N, Vector.average.value, tissue), 
#                     PSEA_e_GO %>% select(Set.ID, Set.N, Vector.average.value, tissue)) %>%
#  mutate(term = gsub("_", " ", Set.ID)) %>% mutate(term=gsub("GOBP", "", term)) %>% mutate(term = term %>% tolower()) %>%
#  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
#                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))
#
#fig1F_GOBP <- ggplot(PSEA_GO %>%
#                       mutate(phase_clock1 = Vector.average.value + 8) %>% 
#                       mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)), 
#                     aes(x=phase_clock1, y=term)) +
#  geom_point(aes(size=Set.N, color=tissue)) +   
#  facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
#  #facet_wrap(~tissue, scales="free_y", space="free_y", ncol=1, nrow=2) +
#  theme_bw() + xlab("phase (h)\n") + ylab(ifelse(PSEA_string=="GOBP", "GO:BP term", "KEGG pathway")) + 
#  guides(color = FALSE) + labs(size='# of genes') + #theme(aspect.ratio=4/3) + 
#  theme(axis.line = element_line(colour = "black"),
#        panel.border = element_blank(),
#        panel.background = element_blank(),
#        strip.background = element_blank(),
#        #strip.text = element_text(face="bold"))
#        strip.text = element_blank(),
#        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
#  ggtitle("PSEA analysis\n") +  expand_limits(y=0) + 
#  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) 


# Figure 1G: 
# GO analysis (this is complementary to the PSEA, since PSEA takes no background and GO does)
gD <- goana(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
gE <- goana(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
#universe is the 'background' universe for the analysis, if none provided, then all EntrezID will be used
topGO(gD %>% filter(Ont=="BP"), n=20, truncate.term = "42") # super significant terms have p values in the order of 10e-8
topGO(gE %>% filter(Ont=="BP"), n=20, truncate.term = "42") 

g <- gE %>% filter(Ont=="BP") %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N, tissue="epidermis") %>% as.data.frame() %>%
  rbind(gD %>% filter(Ont=="BP") %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N, tissue="dermis") %>% as.data.frame())

fig2G_GOBP <- ggplot(g, aes(x=hits, y=Term, size=P.DE, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="GO:BP term", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("GO:BP analysis")

# KEGG analysis
kD <- kegga(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
kE <- kegga(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
topKEGG(kD, n=20, truncate.path = "42")
topKEGG(kE, n=20, truncate.path = "42")

k <- kE %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N, tissue="epidermis") %>% as.data.frame() %>%
  rbind(kD %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N, tissue="dermis") %>% as.data.frame())

fig2G_KEGG <- ggplot(k, aes(x=hits, y=Pathway, size=P.DE, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="KEGG pathway", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("KEGG pathway analysis")

#################################
#################################

# Supp Fig7A: Venn Diagram of number of rhythmic genes
yrhy <- list(dermis = yave$genes$ProbeName[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol], 
             epidermis = yave$genes$ProbeName[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol])
suppfig7A <- ggvenn(yrhy, fill_color = c("#1B9E77", "#D95F02"), 
                    stroke_size = 0.5, set_name_size = 4)

# Suppfig6B: How many genes (that pass FDR and amp cutoff) have a certain amplitude? Histogram of this
rhy_results$tissue <- factor(rhy_results$tissue, levels = c("epidermis", "dermis"))
suppfig7B <- ggplot(data = rhy_results) + 
  geom_vline(aes(xintercept=amp_cutoff), color="black", linetype="dashed") +
  #geom_histogram(aes(amp_value, color=tissue), position="identity", fill="white", bins=50, alpha=0.5) + 
  geom_histogram(aes(amp_value, fill=tissue), color="black", bins=50, size=0.2) + 
  theme_bw() + 
  theme(aspect.ratio=1,        
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom") + 
  xlab("Amplitude") + ylab("number of genes\n(passing amplitude and fdr cutoff)") +
  scale_fill_brewer(palette="Dark2", direction=-1) + ggtitle("\n")

# Supp Fig7C: Heatmaps of Spearman correlation between each pair of clock genes for dermis and epidermis? (~1B El Athman)
timeseries_cg <- yave[which(yave$genes$Symbol %in% clock_genes),]$E %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  rename(c("ProbeName"="rowname")) %>% gather(tissuetime_subj, value, -ProbeName) %>%
  inner_join(yave$genes %>% as.data.frame() %>% select(ProbeName, Symbol)) %>% select(-ProbeName) %>%
  tidyr::separate(tissuetime_subj, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1)

timeseries_cg_raw <- timeseries_cg
timeseries_cg %<>%
  group_by(tissue, subject, Symbol) %>%
  mutate(value=value-mean(value)) %>% as.data.frame()

#if(substract_mean_spearman == TRUE){
#  #means <- timeseries_cg %>% group_by(tissue, subject, Symbol) %>% dplyr::summarise(mean_value = mean(value))
#  timeseries_cg %<>%
#    group_by(tissue, subject, Symbol) %>%
#    mutate(value=value-mean(value)) %>% as.data.frame()
#}

data_corrmat <- timeseries_cg %>% #mutate(gene_tissue = paste(Symbol, tissue, "_")) %>% 
  spread(Symbol, value)
data_corrmat_raw <- timeseries_cg_raw %>% #mutate(gene_tissue = paste(Symbol, tissue, "_")) %>% 
  spread(Symbol, value)

corrmat_D <- data_corrmat %>% filter(tissue=="D") %>% mutate(subjtime = paste(subject, time, sep="_")) %>%
  select(-tissue, -time, -subject) %>% tibble::column_to_rownames('subjtime') %>% cor(., method="spearman") %>% as.matrix()
corrmat_E <- data_corrmat %>% filter(tissue=="E") %>% mutate(subjtime = paste(subject, time, sep="_")) %>%
  select(-tissue, -time, -subject) %>% tibble::column_to_rownames('subjtime') %>% cor(., method="spearman") %>% as.matrix() 

corrmat <- corrmat_D %>% as.data.frame() %>% mutate(tissue="\ndermis") %>% tibble::rownames_to_column() %>% 
  gather(key, value, -tissue,-rowname) %>%
  full_join(corrmat_E %>% as.data.frame() %>% mutate(tissue="\nepidermis") %>% tibble::rownames_to_column() %>% 
              gather(key, value, -tissue,-rowname))
corrmat$key <- factor(corrmat$key, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "CRY1", 
                                            "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP"))
corrmat$rowname<- factor(corrmat$rowname, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "CRY1", 
                                                   "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP"))

my.lines <- data.frame(x=c(.5,6.5), y=c(6.5,.5), 
                       xend=c(13.5,6.5), yend=c(6.5,13.5))
suppfig7C <- ggplot(data = corrmat, aes(x=key, y=rowname, fill=value)) + facet_wrap(~tissue) +
  geom_tile() + theme_bw() + xlab("") + ylab("") + labs(fill=expression("Spearman's"~rho)) +
  scale_fill_distiller(palette = "RdBu", limits=c(-1.05,1.05), breaks = c(1, 0.5, 0, -0.5, -1)) + 
  #geom_hline(yintercept = "CRY1"+1) +
  theme(aspect.ratio=1,
        axis.text.x = element_text(angle = 90, vjust = 0.5, face="italic", hjust=1),
        axis.text.y = element_text(face="italic"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        legend.title = element_text(vjust = 0.8),
        legend.position="bottom") + 
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F) 


# Supp Fig7D: Scatterplot of FC Amp vs acrophase and histograms (~1E El Athman)
toplot <- rhy_results %>% filter(tissue=="dermis") %>% mutate(FC_amp = 2^(2*amp_value)) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))
suppfig7D_1 <- ggplot(toplot, aes(x=phase_clock1, y=amp_value, color=tissue)) +
  geom_point(alpha=0.3) +
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=phase_clock1, y=amp_value), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.5, point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + scale_y_continuous(trans='log2') + 
  #facet_wrap(~tissue, nrow=2, ncol=1, scales="free") +
  ylab("amplitude") + xlab("phase (h)") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold", size=10)) + guides(color = FALSE) + 
  ggtitle("\ndermis") #+ ylim(2^(2*amp_cutoff),4.5)
suppfig7D_1 <- ggExtra::ggMarginal(suppfig7D_1, type = 'histogram', fill="#1B9E77", color="white")

toplot <- rhy_results %>% filter(tissue=="epidermis") %>% mutate(FC_amp = 2^(2*amp_value)) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))
suppfig7D_2 <- ggplot(toplot, aes(x=phase_clock1, y=amp_value, color=tissue)) +
  geom_point(alpha=0.3, color="#D95F02") +
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=phase_clock1, y=amp_value), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.5, point.padding=.5, size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + scale_y_continuous(trans='log2')+
  #facet_wrap(~tissue, nrow=2, ncol=1, scales="free") +
  ylab("amplitude") + xlab("phase (h)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold", size=10)) + guides(color = FALSE) + 
  ggtitle("\nepidermis") #+ ylim(2^(2*amp_cutoff),4.5)
suppfig7D_2 <- ggExtra::ggMarginal(suppfig7D_2, type = 'histogram', color="white", fill="#D95F02")

suppfig7D <- plot_grid(suppfig7D_1, suppfig7D_2, ncol=1, nrow=2)

#############
#############

# SUPP FIG 4&5: Comparison of GO/KEGG (respectively) done by different means
# (limma, clusterProfiler::enrichGO/KEGG, clusterProfiler with msigdbr)

#GO and KEGG with the clusterProfiler package
go_D_1 <- enrichGO(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID, org.Hs.eg.db, 
                   ont = "BP", universe = yave$genes$EntrezID, pvalueCutoff = 1.0, 
                   qvalueCutoff = 1.0)
go_E_1 <- enrichGO(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID, org.Hs.eg.db, 
                  ont = "BP", universe = yave$genes$EntrezID, pvalueCutoff = 1.0, 
                   qvalueCutoff = 1.0)
kegg_D_1 <- enrichKEGG(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID, 
                       organism="hsa", universe = yave$genes$EntrezID, keyType = "ncbi-geneid", pvalueCutoff = 1.0, 
                       qvalueCutoff = 1.0, minGSSize = 5)
kegg_E_1 <- enrichKEGG(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID, 
                       organism="hsa", universe = yave$genes$EntrezID, keyType = "ncbi-geneid", pvalueCutoff = 1.0, 
                       qvalueCutoff = 1.0, minGSSize = 5)

#head(go_D_1, n=20) %>%
#  mutate(Description = str_wrap(Description, width=20),
#         qvalue = scales::scientific(qvalue, digits = 3),
#         GeneRatio = GeneRatio) %>%
#  dplyr::select(Description, qvalue, GeneRatio) %>%
#  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE)
#head(go_E_1, n=20) %>%
#  mutate(Description = str_wrap(Description, width=20),
#         qvalue = scales::scientific(qvalue, digits = 3),
#         GeneRatio = GeneRatio) %>%
#  dplyr::select(Description, qvalue, GeneRatio) %>%
#  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE)
#head(kegg_D_1, n=20) %>%
#  mutate(Description = str_wrap(Description, width=20),
#         qvalue = scales::scientific(qvalue, digits = 3),
#         GeneRatio = GeneRatio) %>%
#  dplyr::select(Description, qvalue, GeneRatio) %>%
#  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE)
#head(kegg_E_1, n=20) %>%
#  mutate(Description = str_wrap(Description, width=20),
#         qvalue = scales::scientific(qvalue, digits = 3),
#         GeneRatio = GeneRatio) %>%
#  dplyr::select(Description, qvalue, GeneRatio) %>%
#  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) #similar results to limma::goana/kegga analysis

kegg_1 <- head(kegg_E_1, n=20) %>%
  mutate(Description = str_wrap(Description, width=20),
         #qvalue = scales::scientific(qvalue, digits = 3),
         GeneRatio = GeneRatio) %>%
  arrange(as.numeric(pvalue)) %>%
  dplyr::select(Description, pvalue, GeneRatio) %>%
  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
  as.data.frame() %>% mutate(hits=Count*100/N, tissue="epidermis") %>% 
  rbind(head(kegg_D_1, n=20) %>%
          mutate(Description = str_wrap(Description, width=20),
                 #qvalue = scales::scientific(qvalue, digits = 3),
                 GeneRatio = GeneRatio) %>%
          dplyr::select(Description, pvalue, GeneRatio)%>%
          arrange(as.numeric(pvalue)) %>%
          tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
          as.data.frame() %>% mutate(hits=Count*100/N, tissue="dermis")) %>%
  mutate(Description = gsub("\n", " ", Description))

go_1 <- head(go_E_1, n=20) %>%
  mutate(Description = str_wrap(Description, width=20),
         #qvalue = scales::scientific(qvalue, digits = 3),
         GeneRatio = GeneRatio) %>%
  arrange(as.numeric(pvalue)) %>%
  dplyr::select(Description, pvalue, GeneRatio) %>%
  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
  as.data.frame() %>% mutate(hits=Count*100/N, tissue="epidermis") %>% 
  rbind(head(go_D_1, n=20) %>%
          mutate(Description = str_wrap(Description, width=20),
                 #qvalue = scales::scientific(qvalue, digits = 3),
                 GeneRatio = GeneRatio) %>%
          dplyr::select(Description, pvalue, GeneRatio)%>%
          arrange(as.numeric(pvalue)) %>%
          tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
          as.data.frame() %>% mutate(hits=Count*100/N, tissue="dermis")) %>%
  mutate(Description = gsub("\n", " ", Description))

suppfig8A <- fig2G_GOBP + ggtitle("GO:BP Analysis with limma")
suppfig8B <- ggplot(go_1, aes(x=hits, y=Description, size=pvalue, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="GO:BP Term", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("GO:BP Analysis with clusterProfiler")

suppfig9A <- fig2G_KEGG + ggtitle("KEGG pathway analysis with limma")
suppfig9B <- ggplot(kegg_1, aes(x=hits, y=Description, size=pvalue, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="KEGG pathway", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("KEGG pathway analysis with clusterProfiler")

#GO and KEGG with the clusterProfiler package and msigdb thing
go_D_2 <- enricher(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$Symbol, 
                   universe = yave$genes$Symbol, 
                   TERM2GENE = m_c5 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame(), 
                   qvalueCutoff = 1.0, pAdjustMethod = "none")
go_E_2 <- enricher(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$Symbol, 
                   universe = yave$genes$Symbol, 
                   TERM2GENE = m_c5 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame(), 
                   qvalueCutoff = 1.0, pAdjustMethod = "none")
kegg_D_2 <- enricher(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$Symbol, 
                     universe = yave$genes$Symbol, 
                     TERM2GENE = m_c2 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame(), 
                     qvalueCutoff = 1.0, pAdjustMethod = "none")
kegg_E_2 <- enricher(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$Symbol, 
                     universe = yave$genes$Symbol, 
                     TERM2GENE = m_c2 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame(), 
                     qvalueCutoff = 1.0, pAdjustMethod = "none")

go_2 <- head(go_E_2, n=20) %>%
  mutate(Description = str_wrap(Description, width=20),
         #qvalue = scales::scientific(qvalue, digits = 3),
         GeneRatio = GeneRatio) %>%
  arrange(as.numeric(pvalue)) %>%
  dplyr::select(Description, pvalue, GeneRatio) %>%
  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
  as.data.frame() %>% mutate(hits=Count*100/N, tissue="epidermis") %>% 
  rbind(head(go_D_1, n=20) %>%
          mutate(Description = str_wrap(Description, width=20),
                 #qvalue = scales::scientific(qvalue, digits = 3),
                 GeneRatio = GeneRatio) %>%
          dplyr::select(Description, pvalue, GeneRatio)%>%
          arrange(as.numeric(pvalue)) %>%
          tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
          as.data.frame() %>% mutate(hits=Count*100/N, tissue="dermis")) %>%
  mutate(Description = gsub("\n", " ", Description)) %>% mutate(term = gsub("_", " ", Description)) %>% 
  mutate(term=gsub("GOBP", "", term)) %>%
  mutate(term = term %>% tolower()) %>%
  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))
kegg_2 <- head(kegg_E_2, n=20) %>%
  mutate(Description = str_wrap(Description, width=20),
         #qvalue = scales::scientific(qvalue, digits = 3),
         GeneRatio = GeneRatio) %>%
  arrange(as.numeric(pvalue)) %>%
  dplyr::select(Description, pvalue, GeneRatio) %>%
  tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
  as.data.frame() %>% mutate(hits=Count*100/N, tissue="epidermis") %>% 
  rbind(head(kegg_D_1, n=20) %>%
          mutate(Description = str_wrap(Description, width=20),
                 #qvalue = scales::scientific(qvalue, digits = 3),
                 GeneRatio = GeneRatio) %>%
          dplyr::select(Description, pvalue, GeneRatio)%>%
          arrange(as.numeric(pvalue)) %>%
          tidyr::separate(GeneRatio, c("Count","N"), sep = "/", convert = TRUE) %>% 
          as.data.frame() %>% mutate(hits=Count*100/N, tissue="dermis")) %>%
  mutate(Description = gsub("\n", " ", Description)) %>% mutate(term = gsub("_", " ", Description)) %>% 
  mutate(term=gsub("KEGG", "", term)) %>%
  mutate(term = term %>% tolower()) %>%
  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))

suppfig8C <- ggplot(go_2, aes(x=hits, y=term, size=pvalue, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="GO:BP Term", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("GO:BP Analysis with clusterProfiler and MSigDB C5 GO:BP")

suppfig9C <- ggplot(kegg_2, aes(x=hits, y=term, size=pvalue, color=tissue)) + geom_point() + scale_size(trans='reverse') +
  #facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0,105), y=0) + 
  labs(x="Hits (%)", y="KEGG pathway", size=expression(italic('p')~'value')) + guides(color = FALSE) +
  theme_bw() + 
  theme(aspect.ratio=2.3,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle("KEGG pathway analysis with clusterProfiler and MSigDB C2 CP:KEGG")

#################################
#################################

#AK's suggestion: how are amplitudes/phases different in Dermis vs Epidermis?
df = rhy_results[which(rhy_results$Symbol %in% rhy_results[duplicated(rhy_results$Symbol),"Symbol"]),] %>% 
  arrange(desc(amp_value))
df_amp   = df %>% dplyr::select(Symbol, tissue, amp_value) %>% spread(tissue, amp_value) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% arrange(desc(dermis))
df_phase = df %>% dplyr::select(Symbol, tissue, phase_value) %>% spread(tissue, phase_value) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))

suppfig10A <- ggplot(df_amp, aes(x=dermis, y=epidermis)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="grey47", alpha=0.5) +
  geom_point(data = filter(df_amp, Symbol %in% clock_genes), 
             aes(x=dermis, y=epidermis), color="red") +
  geom_text_repel(data = filter(df_amp, Symbol %in% clock_genes), 
                  aes(x=dermis, y=epidermis, label=Symbol_it), color="red", parse=TRUE) +
  coord_fixed() + theme_bw() + xlab("amplitude dermis") + ylab("amplitude epidermis") +
  theme(aspect.ratio=1.0,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) 
print(paste0(which(df_amp$epidermis > df_amp$dermis) %>% length(), "/", dim(df_amp)[1], 
             " genes (rhythmic in both tissues) with higher amplitude in epidermis than dermis"))

suppfig10B <- ggplot(df_phase, aes(x=dermis, y=epidermis)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="grey47", alpha=0.5) +
  geom_point(data = filter(df_phase, Symbol %in% clock_genes), 
             aes(x=dermis, y=epidermis), color="red") +
  geom_text_repel(data = filter(df_phase, Symbol %in% clock_genes), 
                  aes(x=dermis, y=epidermis, label=Symbol_it), color="red", parse=TRUE) +
  coord_fixed() + theme_bw() + xlab("acrophase dermis (h)") + ylab("acrophase epidermis (h)") +
  theme(aspect.ratio=1.0,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  scale_x_continuous(breaks=c(-12,-6,0,6,12), labels=c(-12,-6,0,6,12)) +
  scale_y_continuous(breaks=c(-12,-6,0,6,12), labels=c(-12,-6,0,6,12))      

#################################
#################################

# ARRANGE PLOTS IN GRID
# ---------------------

# Fig1
#fig2F <- fig2F_KEGG
fig2G <- fig2G_GOBP

#fig1_top <- plot_grid(NULL, fig1B, fig1C, ncol=3, nrow=1, labels=c("A", "B", "C"), rel_widths = c(1.35,0.8,0.83))
#fig1_center <- plot_grid(fig1D, suppfig1D, fig1F, ncol=3, nrow=1, labels=c("D", "E", "F"), rel_widths = c(1,1.2,2.5))
#fig1_bottom <- plot_grid(fig1G, labels=c("G"))
#fig1 <- plot_grid(fig1_top, fig1_center, fig1_bottom, align='v', nrow=3, rel_heights = c(1, 1.2, 1.6))


#fig1_top1 <- plot_grid(NULL, NULL, fig1B, ncol=3, nrow=1, labels=c("A", "", "B"), rel_widths = c(1.35,0.4,0.9))
#fig1_top2 <- plot_grid(fig1C,NULL, fig1D + facet_wrap(~tissue, ncol=2, nrow=1), ncol=3, labels=c("C", "", "D"), 
#                       rel_widths = c(0.8,0.01,1.6))
#fig1_bottom1 <- plot_grid(suppfig1D, fig1F, ncol=2, nrow=1, labels=c("E", "F"), rel_widths = c(2,2.55))
#fig1_bottom2 <- plot_grid(fig1G, labels=c("G"))
#fig1 <- plot_grid(fig1_top1, fig1_top2, fig1_bottom1, fig1_bottom2, align='v', nrow=4, rel_heights = c(1,1, 1.2, 1.6))

fig2_top2 <- plot_grid(fig2C, NULL, ncol=2, labels=c("C", "D"), 
                       rel_widths = c(2,2.55))
fig2_bottom1 <- plot_grid(fig2D, suppfig7D, ncol=2, nrow=1, labels=c("E", "F"), rel_widths = c(2,2.55))
fig2_bottom2 <- plot_grid(fig2G, labels=c("G"))
fig2 <- plot_grid(fig2_top2, fig2_bottom1, NULL, fig2_bottom2, align='v', nrow=4, 
                  rel_heights = c(1.2, 1.8, 0.02, 1.6))

# Supp Fig1
suppfig7_top <- plot_grid(suppfig7A, NULL, suppfig7B, nrow=1, ncol=3, rel_widths = c(1,0.05,1),labels=c("A", "", "B")) #rel_widths = c(1,1.75), 
suppfig7_bottom <- plot_grid(suppfig7C, NULL, 
                             fig2E + facet_wrap(~tissue, ncol = 2) + theme(legend.position="bottom",
                                                                           legend.title = element_text(vjust = 0.8)), 
                             nrow=1, ncol=3, rel_widths = c(1,0.05,1), labels = c("C", "", "D"))
suppfig7 <- plot_grid(suppfig7_top, NULL, suppfig7_bottom, nrow=3, ncol=1, rel_heights = c(1,0.01,1))


# Supp Fig2:
suppfig8 <- plot_grid(suppfig8A, suppfig8B, suppfig8C, nrow=3, ncol=1, labels=c("A", "B", "C"))
suppfig9 <- plot_grid(suppfig9A, suppfig9B, suppfig9C, nrow=3, ncol=1, labels=c("A", "B", "C"))
suppfig10 <- plot_grid(suppfig10A, NULL, suppfig10B, nrow=1, ncol=3, labels=c("A", "", "B"), rel_widths = c(1,0.05,1))

# SAVE FIGURES
# ------------
if (!file.exists("figures/fig2.pdf")){ 
  fig2 %>% ggsave("figures/fig2_GOBP.pdf", ., width = 11,height = 15.5) 
}
if (!file.exists("figures/suppfig7.pdf")){ 
  suppfig7 %>% ggsave("figures/suppfig7.pdf", ., width=11, height=9) 
}
if (!file.exists("figures/suppfig8.pdf")){ 
  suppfig8 %>% ggsave("figures/suppfig8.pdf", ., width = 11,height = 15.5) 
}
if (!file.exists("figures/suppfig9.pdf")){ 
  suppfig9 %>% ggsave("figures/suppfig9.pdf", ., width = 11,height = 15.5) 
}
if (!file.exists("figures/suppfig10.pdf")){ 
  suppfig10 %>% ggsave("figures/suppfig10.pdf", ., width = 11,height = 5.5) 
}
############
############


results_fig1 <- readRDS("visualize/data/results_fig1.rds")
results_fig1 %<>% 
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)
rhy_some_fig1 <- results_fig1 %>% filter(pmax(A_D, A_E) > amp_cutoff & adj_P_Val < fdr_cutoff)
rhy_both_fig1 <- rhy_some_fig1 %>% filter(A_D > amp_cutoff & A_E > amp_cutoff)

rhy_some_fig2 <- results %>% filter(pmax(A_D, A_E) > amp_cutoff & adj_P_Val < fdr_cutoff)
rhy_both_fig2 <- rhy_some_fig2 %>% filter(A_D > amp_cutoff & A_E > amp_cutoff)


rhy <- list(dermis_fig1 = filter(rhy_some_fig1, A_D > amp_cutoff)$Symbol,
            dermis_fig2 = filter(rhy_some_fig2, A_D > amp_cutoff)$Symbol,
            epidermis_fig1 = filter(rhy_some_fig1, A_E > amp_cutoff)$Symbol,
            epidermis_fig2 = filter(rhy_some_fig2, A_E > amp_cutoff)$Symbol)

# Overlaps between rhythmic genes found in analysis from Fig1 and 2  
ggvenn(rhy[1:2], fill_color = c("#00798c", "#edae49"), 
       stroke_size = 0.5, set_name_size = 4) + ggtitle("Rhythmic genes in Fig1 and Fig2 in dermis: overlaps")
dermis_genes_fig2_notin_fig1 <- rhy$dermis_fig2[which(!(rhy$dermis_fig2 %in% rhy$dermis_fig1))]

ggvenn(rhy[3:4], fill_color = c("#00798c", "#edae49"), 
       stroke_size = 0.5, set_name_size = 4) + ggtitle("Rhythmic genes in Fig1 and Fig2 in epidermis: overlaps")
epidermis_genes_fig2_notin_fig1 <- rhy$epidermis_fig2[which(!(rhy$epidermis_fig2 %in% rhy$epidermis_fig1))]

# Correlation of amplitudes of rhythmic genes found in analysis from Fig1 and 2
rhy_some_figs12 <- inner_join(rhy_some_fig1 %>% dplyr::select(Symbol, A_D, A_E, phaseD, phaseE) %>%
                                dplyr::rename(c("A_D_fig1"="A_D", "A_E_fig1"="A_E", "phaseD_fig1"="phaseD", "phaseE_fig1"="phaseE")),
                              rhy_some_fig2 %>% dplyr::select(Symbol, A_D, A_E, phaseD, phaseE) %>%
                                dplyr::rename(c("A_D_fig2"="A_D", "A_E_fig2"="A_E", "phaseD_fig2"="phaseD", "phaseE_fig2"="phaseE")))
rhy_some_figs12 <- full_join(rhy_some_figs12 %>% select(Symbol, contains("A_")) %>% 
                               gather(key, amp_value, -Symbol) %>%
                               tidyr::separate(key, c("junk", "tissue", "fig"), sep = "_", convert = TRUE) %>% select(-junk) %>%
                               spread(fig, amp_value) %>%
                               dplyr::rename(c("amp_fig1"="fig1", "amp_fig2"="fig2")) %>%
                               mutate(tissue = ifelse(tissue=="D", "dermis", "epidermis")),
                             
                             rhy_some_figs12 %>% select(Symbol, contains("phase")) %>% 
                               gather(key, phase_value, -Symbol) %>%
                               tidyr::separate(key, c("tissue", "fig"), sep = "_", convert = TRUE) %>% 
                               spread(fig, phase_value) %>%
                               dplyr::rename(c("phase_fig1"="fig1", "phase_fig2"="fig2")) %>%
                               mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis")))

ggplot( rhy_some_figs12, aes(x=amp_fig1, y=amp_fig2, color=tissue) ) + geom_point(alpha=0.5) +
  coord_fixed() + theme_bw() + 
  xlab("amplitude estimated from analysis fig 1 (a.u.)") + ylab("amplitude estimated from analysis fig 2 (a.u.)") +
  theme(aspect.ratio=1.0,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) + 
  ggtitle('Amplitude correlation analysis figs1 and 2\n(genes that are rhythmic in at least one tissue)')

ggplot( rhy_some_figs12, aes(x=phase_fig1, y=phase_fig2, color=tissue) ) + geom_point(alpha=0.5) +
  coord_fixed() + theme_bw() + 
  xlab("acrophase estimated from analysis fig 1 (h)") + ylab("acrophase estimated from analysis fig 2 (h)") +
  theme(aspect.ratio=1.0,
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) + 
  ggtitle('Acrophase correlation analysis figs1 and 2\n(genes that are rhythmic in at least one tissue)')
