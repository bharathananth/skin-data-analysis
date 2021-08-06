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
library(lubridate)
library(gridExtra)
library(data.table)
library(tibble)
library(circular)
library(scales)

setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis/")
source("funcs.R")

#--------------------------------
#--------------------------------

# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 
removal_outliers <- "PCA" #PCA or weights
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R
save_fitdata <- FALSE
n_cutoff <- 11#8

#--------------------------------
#--------------------------------

# 2. READ FILES
# -------------
# Read info of subjects, calculate mid sleeping time
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
  rename(c("subject"="Subject", "sex"="Sex"))
info_subjects <- info_subjects_long %>% select(subject, sex, Light_condition, age, MSF_sc)

# Raw data (annotated, normalized, etc)
yave <- readRDS("visualize/data/rawdata.rds")

# Extract sample details from column names
experiment <- readRDS("visualize/data/experiment.rds")
experiment %<>% full_join(info_subjects)

#--------------------------------
#--------------------------------

# 3. DEALING WITH OUTLIERS: REMOVING MANUALLY THE OUTLIERS SEEN IN PCA (OPTION1) OR WEIGHING SAMPLES WITH ARRAYWEIGHTS (OPTION2)
# ------------------------------------------------------------------------------------------------------------------------------
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

# Prepare sample details for design matrix
# ----------------------------------------
tissue  <- factor(experiment$tissue)
time    <- experiment$time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + b*sin(wt) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

#--------------------------------
#--------------------------------

# 4. CREATE DESIGN MATRIX AND FIT DATA
# ------------------------------------
# design matrix
design <- model.matrix(~ 0 + subject + tissue + subject:tissue:inphase + subject:tissue:outphase) #H0: rhythms are different across subjects and tissues

# weights + fits
fit <- limma::lmFit(yave, design, weights = wts) #!!!!
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

if (save_fitdata == TRUE){
  # save results of fits -> topTable is done in each subject individually
  for (i in unique(subject)){
    for (j in unique(tissue)){
      rhy_indices <- which(grepl(paste0(i,".*tissue",j,".*phase"), colnames(design)))
      results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
        set_colnames(gsub("\\.","_", colnames(.))) %>%
        rename(!! paste0("P_Value_", i) := P_Value) %>%
        rename(!! paste0("adj_P_Val_", i) := adj_P_Val) 
      saveRDS(results, file = paste0("visualize/data/results_fig2_", i, "_", j, ".rds"))
    }
  }
}

#--------------------------------
#--------------------------------

# 5. READ RESULTS AND CALCULATE AMPLITUDES, PHASES 
# ------------------------------------------------
i <- "P115" # from unique(subject)
j <- "D"
resultsi <- readRDS(paste0("visualize/data/results_fig2_", i, "_", j, ".rds")) %>% select(-AveExpr, -F)

results <- matrix(0, nrow=11578) %>% as.data.frame()
for (i in unique(subject)){
  for (j in unique(tissue)){
    # read results of each patient (TopTable was done in each patient individually, see notebook 26.03.2021)
    results_i <- readRDS(paste0("visualize/data/results_fig2_", i, "_", j, ".rds")) %>% select(-AveExpr, -F) %>%
      rename(!! paste0("P_Value_",i,"_",j) := paste0("P_Value_",i)) %>%
      rename(!! paste0("adj_P_Val_",i,"_",j) := paste0("adj_P_Val_",i))
    # calculate amplitudes and phases from fits
    results_i[,paste0("A_", j, "_", i)] <- sqrt(results_i[,5]^2 + results_i[,6]^2)
    results_i[,paste0("phase", j, "_", i)] <- atan2(results_i[,6], results_i[,5])*12/pi
    results_i[,paste0("rhythm_", i, "_",j)] <- complex(modulus=results_i[,9], argument=results_i[,10]*pi/12) #argument in rad!!
  
    # stack all columns together
    if (i == "P100" & j == "D"){
      results <- cbind(results, results_i) %>% select(-V1)
    } else{
      results <- full_join(results, results_i)
    }
  }
}

#--------------------------------
#--------------------------------

# 6. DATA FRAME OF COMPLEX NUMBERS (~RHYTHM) AND CALCULATION OF RSS
# -----------------------------------------------------------------

# Data frame of complex numbers
results_complex <- results %>% select(ProbeName, Symbol, contains("rhythm"))
temp <- results_complex %>% gather(key, rhythm_complex, -ProbeName, -Symbol) %>%
  separate(key, c("junk","subject","tissue"), sep = "_", convert = TRUE) %>% select(-junk) 

# Determine fitted values and residuals
all(fitted.values(fit) == fitted.values(fit2)) #input check that fitted values are the same in fit and fit2
all(residuals(fit, yave) == residuals(fit2, yave)) #input check that residuals (expression_value - fitted_value) are the same in fit and fit2

fitted_values <- fitted.values(fit) %>% as.data.frame()
colnames(fitted_values) <- colnames(yave$E)
resid_values  <- residuals(fit, yave) %>% as.data.frame()

rss <- resid_values %>% rownames_to_column() %>% gather(key, value, -rowname) %>%
  separate(key, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
  dplyr::group_by(rowname, tissue, subject) %>%
  dplyr::summarise(rss = sum(value^2)) %>% as.data.frame() %>% 
  mutate(rse = sqrt(rss/fit$sigma[1])) #different values to the fit$sigma (probably because for fit$sigma what was counted is the RSS in D&E)
names(rss)[1] <- "ProbeName"

dim(rss)[1] == dim(yave)[1]*length(unique(tissue))*length(unique(subject)) #input check: as many RSS as total_no_genes*11subjects*2tissues

#--------------------------------
#--------------------------------

# 7. DATA FRAME OF WEIGHTED RHYTHMS
# ---------------------------------

# Calculate weighted rhythms (modulus_rhy / rss)
rhythms_weighted <- temp %>% full_join(rss) %>% #rhythm weighed == rhythm (in complex form) divided by the rss (serves as "weight")
  mutate(rhythm_weighted = rhythm_complex / rss, # I choose rss, easier to justify
         modulus_rhy = Mod(rhythm_complex),
         argument_rhy = Arg(rhythm_complex), #same value as argument_rhywt :)
         modulus_rhywt = Mod(rhythm_weighted)) 
# NOTES!!!
# 1. Some genes have high modulus (A) and high RSS -> end up with low-ish modulus(weighted_rhy)
# 2. Some genes have initially low modulus (A) and not too bad RSS -> end up with high modulus(weighted_rhy)

#--------------------------------
#--------------------------------

# 8. DETERMINATION OF RHYTHMIC GENES ACROSS SUBJECTS
# --------------------------------------------------
# 2 step filtering strategy to end up with rhythmic genes
rhythms_weighted_filter <- rhythms_weighted %>% 
  filter(modulus_rhywt >= amp_cutoff) %>% #filter1: by weighted modulus (this way genes with very high RSS are removed)
  filter(modulus_rhy >= amp_cutoff) #filter2: still we might have the NOTE2 scenario (genes with originally very low modulus)
rhythms_weighted_filter$Symbol %>% unique %>% length #we have a potpourri of basically all genes..

counts <- rhythms_weighted_filter %>% group_by(tissue) %>% count(Symbol) %>% as.data.frame()
rhythms_weighted_filter <- full_join(counts, rhythms_weighted_filter)

# Plots of rhythmic genes across subjects/tissues
plot_bar1 <- ggplot(rhythms_weighted_filter %>% select(Symbol, subject, tissue) ) + 
  geom_bar(aes(x=subject, fill=tissue)) + facet_grid(~tissue) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Number of genes with |rhy_weighted| > amp_cutoff && |rhy| > amp_cutoff")
plot_bar2 <- ggplot(rhythms_weighted_filter %>% select(Symbol, subject, tissue) ) + 
  geom_bar(aes(x=Symbol, fill=tissue)) + facet_grid(~tissue) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Number of patients in which a gene is rhythmic (|rhy_wt| && |rhy| > amp_cutoff)")
plot_bar3 <- ggplot(rhythms_weighted_filter %>% select(Symbol, n, tissue)) +
  geom_bar(aes(x=n, fill=tissue)) + facet_grid(~tissue) + theme_bw() +
  ggtitle("Number of genes that are rhythmic (|rhy_wt| && |rhy| > amp_cutoff)\nin accumulating number of subjects")

# Keep genes that are rhythmic in at least n_threshold subjects
rhythms_weighted_filter %<>% filter(n >= n_cutoff) 
print(paste0("total # of rhy_genes in D in at least ", n_cutoff, " subjects: ",
             rhythms_weighted_filter %>% filter(tissue=="D") %$% Symbol %>% unique %>% length))
print(paste0("total # of rhy_genes in E in at least ", n_cutoff, " subjects: ",
             rhythms_weighted_filter %>% filter(tissue=="E") %$% Symbol %>% unique %>% length))

# Calculate sum of UNweighted rhythmicity indeces of a gene across subjects (where gene is rhy): 
# high values ~ gene in phase across subjs
# (Weighted scores just to filter rhythmic genes)
sum_rhythms <- rhythms_weighted_filter %>% 
  dplyr::group_by(ProbeName, Symbol, tissue) %>%
  dplyr::summarise(sum_rhy = sum(rhythm_complex)) %>% as.data.frame() %>%
  mutate(modulus_sum = Mod(sum_rhy))

rhythms_weighted_filter %<>% full_join(sum_rhythms) %>% 
  mutate(modulus_sum_norm = modulus_sum/n) %>% #to distinguish if high sum is because gene is rhythmic in all 11 patients vs just 8
  arrange(desc(modulus_sum_norm)) 

# Top genes with highest modulus_sum or modulus_sum_norm in each tissue
rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum)) %$% Symbol %>% unique %>% head
rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %$% Symbol %>% unique %>% head
rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum)) %$% Symbol %>% unique %>% head
rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %$% Symbol %>% unique %>% head

# NOTES!!!
# 3. Regarding the sum of modula to see if a gene is in phase across subj -> if we get a value of high sum, 
#    - how to distinguish if the gene is really in phase across subjs or just has a couple of subjs with high A?

rhythms_weighted_filter %>% select(tissue, Symbol, n, rss, modulus_rhy, modulus_rhywt, modulus_sum) %>% head

# Some plots 
plot_gene("PER2", "D") #gene with min rss
plot_gene("TLR5", "E") #gene with max rss


# separate matrices of D&E
rhythmicity.D <- list(
  rhythm = rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %>% 
    select(ProbeName, Symbol, tissue, subject, rhythm_complex) %>%
    mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
    spread(tissuesubject, rhythm_complex) %>% tibble::column_to_rownames("ProbeName") %>%
    mutate(sum = rowSums(.[2:12], na.rm = TRUE), mod_sum=Mod(sum), mod_sum_norm=mod_sum/rowSums(!is.na(.[2:12]))),
  rss = rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %>% 
    select(ProbeName, Symbol, tissue, subject, rss) %>%
    mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
    spread(tissuesubject, rss) %>% tibble::column_to_rownames("ProbeName"),
  rse = rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %>% 
    select(ProbeName, Symbol, tissue, subject, rse) %>%
    mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
    spread(tissuesubject, rse) %>% tibble::column_to_rownames("ProbeName"),
  rhythm_wtd = rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %>% 
    select(ProbeName, Symbol, tissue, subject, rhythm_weighted) %>%
    mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
    spread(tissuesubject, rhythm_weighted) %>% tibble::column_to_rownames("ProbeName")
)

rhythmicity.E <- list(
 rhythm = rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %>% 
   select(ProbeName, Symbol, tissue, subject, rhythm_complex) %>%
   mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
   spread(tissuesubject, rhythm_complex) %>% tibble::column_to_rownames("ProbeName")%>%
   mutate(sum = rowSums(.[2:12], na.rm = TRUE), mod_sum=Mod(sum), mod_sum_norm=mod_sum/rowSums(!is.na(.[2:12]))),
 rss = rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %>% 
   select(ProbeName, Symbol, tissue, subject, rss) %>%
   mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
   spread(tissuesubject, rss) %>% tibble::column_to_rownames("ProbeName"),
 rse = rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %>%
   select(ProbeName, Symbol, tissue, subject, rse) %>%
   mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
   spread(tissuesubject, rse) %>% tibble::column_to_rownames("ProbeName"),
 rhythm_wtd = rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %>% 
   select(ProbeName, Symbol, tissue, subject, rhythm_weighted) %>%
   mutate(tissuesubject=paste0(subject, "_", tissue)) %>% select(-tissue, -subject) %>%
   spread(tissuesubject, rhythm_weighted) %>% tibble::column_to_rownames("ProbeName")
)

top_phasicgenes_D <- rhythmicity.D$rhythm %>% arrange(desc(mod_sum_norm)) %$% Symbol %>% head
top_phasicgenes_E <- rhythmicity.E$rhythm %>% arrange(desc(mod_sum_norm)) %$% Symbol %>% head

all(rhythmicity.D$rhythm[,2:12] / rhythmicity.D$rss[,2:12] == rhythmicity.D$rhythm_wtd[,2:12], na.rm=TRUE) #input checks
all(rhythmicity.E$rhythm[,2:12] / rhythmicity.E$rss[,2:12] == rhythmicity.E$rhythm_wtd[,2:12], na.rm=TRUE)
all(top_phasicgenes_D == rhythms_weighted_filter %>% filter(tissue=="D") %>% arrange(desc(modulus_sum_norm)) %$% Symbol %>% unique %>% head)
all(top_phasicgenes_E == rhythms_weighted_filter %>% filter(tissue=="E") %>% arrange(desc(modulus_sum_norm)) %$% Symbol %>% unique %>% head)

plot_gene(top_phasicgenes_D[1], "D")
plot_gene(top_phasicgenes_E[1], "E")  #ZBTB16 is rhy in 8subj with max modulus_sum


#########
#########

# 9. CALCULATE DELTA_PHASE BETWEEN EACH SUBJECT AND THE MEAN PHASE ACROSS SUBJECTS (for rhy genes)
# --------------------------------------------------------------------------------
ref_subject <- info_subjects %>% filter(MSF_sc == median(MSF_sc)) %$% subject [1] %>% as.character() #ref subject is the one with mid sleep time in middle (median)
idx_refsubj <- which(subject %>% unique == ref_subject)

# Phases of the rhythmic genes (calculated from the complex numbers, thus in RAD) -> coincides with phases from `results` (in h)
phases_D <- Arg(rhythmicity.D$rhythm[,2:12] %>% as.matrix) %>% as.data.frame() %>% rownames_to_column() %>% 
  full_join(rhythmicity.D$rhythm %>% rownames_to_column() %>% select(rowname, Symbol)) %>% 
  gather(key, value, -Symbol, -rowname) %>% mutate(value=ifelse(value<0, value+2*pi, value)) %>%
  spread(key, value) %>% tibble::column_to_rownames("rowname")
round(phases_D[1,2:12]*12/pi, 4) == 
  round(results[which(results$Symbol=="OVGP1"), grepl("phaseD",colnames(results))] + 24, 4) #input check

phases_E <- Arg(rhythmicity.E$rhythm[,2:12] %>% as.matrix) %>% as.data.frame() %>% rownames_to_column() %>% 
  full_join(rhythmicity.E$rhythm %>% rownames_to_column() %>% select(rowname, Symbol)) %>% 
  gather(key, value, -Symbol, -rowname) %>% mutate(value=ifelse(value<0, value+2*pi, value)) %>%
  spread(key, value) %>% tibble::column_to_rownames("rowname")
round(phases_E[1,2:12]*12/pi, 4) == 
  round(results[which(results$Symbol=="OVGP1"), grepl("phaseE",colnames(results))] + 24, 4) #input check

any(phases_D[,2:12] < 0 | phases_D[,2:12] > 2*pi, na.rm=TRUE) #input check, make sure all phases (rad) are between 0 and 2pi
any(phases_E[,2:12] < 0 | phases_E[,2:12] > 2*pi, na.rm=TRUE)

#phases_D[,2:12] <- phases_D[,2:12]/(2*pi)
#phases_E[,2:12] <- phases_E[,2:12]/(2*pi)

# ----------
# Trying to merge everything
phases <- full_join(phases_D, phases_E) %>% gather(key, value, -Symbol) %>% 
  separate(key, c("subject","tissue"), sep = "_", convert = TRUE) %>%
  mutate(value=round(value,6))

summary_phases_D <- phases %>% filter(tissue=="D") %>% drop_na() %>%
  group_by(Symbol) %>% dplyr::summarise(mean_D = mean(circular(value), na.rm=TRUE),#
                                        median_D = median(circular(value), na.rm=TRUE),
                                        sd_D = sd(circular(value), na.rm=TRUE)) %>% as.data.frame() %>%
  mutate(mean_D = ifelse(mean_D < 0, mean_D+2*pi, mean_D),
         median_D = ifelse(median_D < 0, median_D+2*pi, median_D))
summary_phases_E <- phases %>% filter(tissue=="E") %>% drop_na() %>%
  group_by(Symbol) %>% dplyr::summarise(mean_E = mean(circular(value), na.rm=TRUE),#
                                        median_E = median(circular(value), na.rm=TRUE),
                                        sd_E = sd(circular(value), na.rm=TRUE)) %>% as.data.frame() %>%
  mutate(mean_E = ifelse(mean_E < 0, mean_E+2*pi, mean_E),
         median_E = ifelse(median_E < 0, median_E+2*pi, median_E))

summary_phases <- summary_phases_D %>% full_join( summary_phases_E ) %>% gather(key, value, -Symbol) %>%
  separate(key, c("key","tissue"), sep = "_", convert = TRUE) %>%
  mutate(value=round(value,6)) 
summary_phases %<>% full_join( phases[which(phases$value %in% summary_phases$value),] %>% drop_na() )

plot_phases <- ggplot() + 
  #geom_boxplot(data = phases_D %>% gather(key, value, -Symbol), aes(x=Symbol, y=value)) + #problematic with circular distribs
  geom_point(data = phases, aes(x=Symbol, y=value, color=subject)) + 
  geom_point(data=summary_phases %>% filter(key=="median"), aes(x=Symbol, y=value), 
             color="black", shape=1, size=4) +
  geom_text_repel(data=summary_phases %>% filter(key=="median"), aes(x=Symbol, y=value, label=subject)) +
  facet_wrap(~tissue, ncol=1, nrow=2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Phases of rhythmic genes") + xlab("") + ylab("phase [rad]")

plot_subjects_inmedianphase <- ggplot(data = summary_phases %>% filter(key=="median") %>% drop_na() ) +
  geom_bar(aes(x=subject, fill=tissue)) + facet_wrap(~tissue) + ylab("number of times subject's phase is median phase") +
  theme_bw()

# trying in polar:
rhy_D <- rhythmicity.D$rhythm %>% select(-contains("sum")) %>% gather(key, value, -Symbol)
rhy_E <- rhythmicity.E$rhythm %>% select(-contains("sum")) %>% gather(key, value, -Symbol)
rhy <- rbind(rhy_D, rhy_E) %>% separate(key, c("subject","tissue"), sep = "_", convert = TRUE)

summary_rhy <- rhy %>% # OJO PORQUE MEDIAN AMP MIGHT NOT CORRESPOND TO GENE WITH MEDIAN_PHASE
  group_by(Symbol, tissue) %>% 
  dplyr::summarise(mean_phi   = mean.circular(circular(Arg(value)), na.rm=TRUE),
                   median_phi = median.circular(circular(Arg(value)), na.rm=TRUE),
                   sd_phi     = sd.circular(circular(Arg(value)), na.rm=TRUE)) %>% as.data.frame()

pi_scales <- math_format(.x * pi, format = function(x) x / pi)
plot_phases_polar <- ggplot(data = rhy) +
  geom_point(aes(x=Arg(value), y=1, color=tissue), size=0.5, alpha=0.5, shape=4) + 
  geom_segment(aes(x=Arg(value), y=0, xend=Arg(value), yend=1, color=tissue), alpha=0.5, linetype="dotted") +
  facet_wrap(~Symbol, ncol=7, nrow=8) + #facet_wrap_paginate(~Symbol, ncol=7, nrow = 2, page=i) +
  #geom_point(aes(x=ph_diff),bins = 50, fill=NA, color="grey20") + 
  theme_bw() +
  xlab("Phase_Epidermis - Phase_Dermis") + coord_polar(start = pi/2, direction=-1) +
  geom_point(data=summary_rhy, aes(x=median_phi, y=1, color=tissue)) +
  geom_segment(data=summary_rhy, aes(x=median_phi, y=0, xend=median_phi, yend=1, color=tissue)) +
  scale_x_continuous(labels = pi_scales, breaks = c(0, pi/2, pi, -pi/2)) + #scale_y_continuous(breaks = c("", "", "", ""))
  xlab("Phase (h)") + ylab("") + theme(axis.text.y=element_blank(),
                                       axis.ticks.y=element_blank())

plot_phases_polar %>% ggsave(sprintf("figures/rhythmic_genes/rhy_genes.pdf"),.,width = 7,height = 10) 



#ggplot() + geom_point(data = phases_E %>% gather(key, value, -Symbol), aes(x=Symbol, y=value, color=key)) + 
#  geom_point(data=summary_phases_E, aes(x=Symbol, y=mean), color="black", shape=8, size=4) +
#  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Phases of rhythmic genes in Epidermis") +
#  xlab("") + ylab("phase [rad]")

# Calculate delta phases (in rad) with respect to the median (?) phase across subjects (mean delta_phase == 0!)
# deltaphase < 0 => the ref subject is behind || deltaphase > 0 => ref subject is ahead
deltaphases_D <- phases_D %>% rownames_to_column() %>% gather(key, phase, -Symbol, -rowname) %>% full_join(summary_phases_D) %>%
  mutate(deltaphase = phase - value) %>% select(-phase, -value) %>% 
  mutate(deltaphase = ifelse(deltaphase > pi, deltaphase - 2*pi, deltaphase)) %>% #this line and the one below make sure all phases are between -pi&pi
  mutate(deltaphase = ifelse(deltaphase < -pi, deltaphase + 2*pi, deltaphase)) %>%
  spread(key, deltaphase) %>% 
  column_to_rownames("rowname")

deltaphases_E <- phases_E %>% rownames_to_column() %>% gather(key, phase, -Symbol, -rowname) %>% full_join(summary_phases_E) %>%
  mutate(deltaphase = phase - value) %>% select(-phase, -value) %>% 
  mutate(deltaphase = ifelse(deltaphase > pi, deltaphase - 2*pi, deltaphase)) %>% #this line and the one below make sure all phases are between -pi&pi
  mutate(deltaphase = ifelse(deltaphase < -pi, deltaphase + 2*pi, deltaphase)) %>%
  spread(key, deltaphase) %>% 
  column_to_rownames("rowname")

# This chunk below works for the phase of a reference subject instead of median of phases, but what if a gene is not rhy in the ref?
#deltaphases_D <- phases_D[,1:11] - phases_D[,idx_refsubj]
#deltaphases_D <- deltaphases_D %>% rownames_to_column() %>% 
#  full_join(rhythmicity.D$rhythm %>% rownames_to_column() %>% select(rowname, Symbol)) %>% 
#  gather(key, value, -Symbol, -rowname) %>% 
#  mutate(value = ifelse(value > pi, value - 2*pi, value)) %>% #this line and the one below make sure all phases are between -pi&pi
#  mutate(value = ifelse(value < -pi, value + 2*pi, value)) %>%
#  spread(key, value)  %>% tibble::column_to_rownames("rowname") 
#  
#deltaphases_E <- phases_E[,1:11] - phases_E[,idx_refsubj]
#deltaphases_E <- deltaphases_E %>% rownames_to_column() %>% 
#  full_join(rhythmicity.E$rhythm %>% rownames_to_column() %>% select(rowname, Symbol)) %>% 
#  gather(key, value, -Symbol, -rowname) %>% 
#  mutate(value = ifelse(value > pi, value - 2*pi, value)) %>% #this line and the one below make sure all phases are between -pi&pi
#  mutate(value = ifelse(value < -pi, value + 2*pi, value)) %>%
#  spread(key, value)  %>% tibble::column_to_rownames("rowname") 

# Delta_phase of top_phasicgenes should be quite low
deltaphases_D %>% filter(Symbol %in% top_phasicgenes_D)

# Calculate median delta_phase of rhythmic genes in each subject
summary_deltaphases_D_acrosssubjs <- deltaphases_D %>% gather(key, value, -Symbol) %>%
  group_by(Symbol) %>% dplyr::summarise(value = median(value, na.rm=TRUE)) %>% as.data.frame() %>%
  inner_join(deltaphases_D %>% gather(key, value, -Symbol)) %>%
  separate(key, c("subject","tissue"), sep = "_", convert = TRUE) %>% select(-tissue) %>%
  mutate(topphasic_gene = ifelse(Symbol %in% top_phasicgenes_D, TRUE, FALSE))
summary_deltaphases_E_acrosssubjs <- deltaphases_E %>% gather(key, value, -Symbol) %>%
  group_by(Symbol) %>% dplyr::summarise(value = median(value, na.rm=TRUE))  %>% as.data.frame() %>% 
  inner_join(deltaphases_E %>% gather(key, value, -Symbol)) %>%
  separate(key, c("subject","tissue"), sep = "_", convert = TRUE) %>% select(-tissue) %>%
  mutate(topphasic_gene = ifelse(Symbol %in% top_phasicgenes_E, TRUE, FALSE))
# Delta_phase of top_phasicgenes is quite low

summary_deltaphases_D <- deltaphases_D %>% gather(key, value, -Symbol) %>%
  group_by(key) %>% dplyr::summarise(value = median(value, na.rm=TRUE)) %>% as.data.frame()
summary_deltaphases_E <- deltaphases_E %>% gather(key, value, -Symbol) %>%
  group_by(key) %>% dplyr::summarise(value = median(value, na.rm=TRUE)) %>% as.data.frame()

ggplot() + geom_hline(yintercept = 0, linetype="dashed") +
  geom_point(data = deltaphases_D %>% gather(key, value, -Symbol), aes(x=Symbol, y=value, color=key)) + 
  geom_point(data=summary_deltaphases_D_acrosssubjs, aes(x=Symbol, y=value, label=subject, shape=topphasic_gene), 
             size = 3, show.legend = FALSE) + scale_shape_manual(values=c(0,2)) +
  geom_text_repel(data=summary_deltaphases_D_acrosssubjs, aes(x=Symbol, y=value, label=subject), color=NA) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Delta_phases of rhythmic genes in Dermis") +
  xlab("") + ylab("delta phase [rad]\n(phase - median_phase)")
ggplot() + geom_hline(yintercept = 0, linetype="dashed") +
  geom_point(data = deltaphases_E %>% gather(key, value, -Symbol), aes(x=Symbol, y=value, color=key)) + 
  geom_point(data=summary_deltaphases_E_acrosssubjs, aes(x=Symbol, y=value, label=subject, shape=topphasic_gene), 
             size = 3, show.legend = FALSE) + scale_shape_manual(values=c(0,2)) +
  geom_text_repel(data=summary_deltaphases_E_acrosssubjs, aes(x=Symbol, y=value, label=subject), color=NA) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Delta_phases of rhythmic genes in Epidermis") +
  xlab("") + ylab("delta phase [rad]\n(phase - median_phase)")

####################
####################

# 10. PLOT MEDIAN DELTA_PHASE VS. MEAN SLEEP TIME
# --------------------------------------------
summary_deltaphases <- rbind(summary_deltaphases_D, summary_deltaphases_E) %>% as.data.frame() %>%
  separate(key, c("subject","tissue"), sep = "_", convert = TRUE) 
#unlisted_info <- unlist(strsplit(summary_rhygenes$key, "_"))
#summary_rhygenes$tissue <- unlisted_info[seq(3, length(unlisted_info), 4)]
#summary_rhygenes$subject <- unlisted_info[seq(4, length(unlisted_info), 4)]

summary_deltaphases <- summary_deltaphases %>% full_join(info_subjects)

lims <- as.POSIXct(strptime(c("1970-01-01 02:30:00","1970-01-01 07:15:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
fig2A <- ggplot(summary_deltaphases, aes(as.POSIXct(x = MSF_sc, format = "%H:%M"), y=value, 
                                         color=tissue, shape=sex, label=subject)) +
  geom_point(size=3, alpha=0.5) +
  geom_text_repel() +
  scale_x_datetime(date_breaks = "1 hours", date_labels = "%H:%M", limits=lims) +
  xlab("Mid sleep time") + 
  ylab(paste0("median delta_phase of rhythmic genes (h)\n(delta_phase = phase-mean_phase)",
              "\n(proxy for chronotype)")) + theme_bw() + theme(aspect.ratio=1) + 
  ggtitle(paste0("Rhythm filtration: ",
                 "\n|rhy| & |rhy_wtd| > amp_cutoff & n_cutoff=", n_cutoff, 
                 "\nrhygenes_D=", dim(rhythmicity.D$rhythm)[1], ", rhygenes_E=", dim(rhythmicity.E$rhythm)[1]))

if (!file.exists("figures/fig2A.pdf")){ 
  fig2A %>% ggsave("figures/fig2A.pdf", ., width = 5, height = 5) 
}




##########
##########

# Clock genes that are rhythmic (mod_rhy & mod_rhywt > amp_cutoff) in at least 8 subjects:
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
clockgenes_weighted_filter <- rhythms_weighted_filter[which(rhythms_weighted_filter$Symbol %in% clock_genes),]

print(paste0("Rhythmic clock genes in D in at least ", n_cutoff, " subjects: ",
             list(clockgenes_weighted_filter %>% filter(tissue=="D") %$% Symbol %>% unique %>% sort)))
print(paste0("Rhythmic clock genes in E in at least ", n_cutoff, " subjects: ",
             list(clockgenes_weighted_filter %>% filter(tissue=="E") %$% Symbol %>% unique %>% sort)))

#plot_gene("ARNTL", "D") %>% ggsave("figures/clockgenes/ARNTL_D.png", ., width = 8, height = 5)
#plot_gene("ARNTL2","D") %>% ggsave("figures/clockgenes/ARNTL2_D.png", ., width = 8, height = 5)
#plot_gene("CLOCK", "D") %>% ggsave("figures/clockgenes/CLOCK_D.png", ., width = 8, height = 5)
#plot_gene("CRY1",  "D") %>% ggsave("figures/clockgenes/CRY1_D.png", ., width = 8, height = 5)
#plot_gene("CRY2",  "D") %>% ggsave("figures/clockgenes/CRY2_D.png", ., width = 8, height = 5)
#plot_gene("DBP",   "D") %>% ggsave("figures/clockgenes/DBP_D.png", ., width = 8, height = 5)
#plot_gene("NPAS2", "D") %>% ggsave("figures/clockgenes/NPAS2_D.png", ., width = 8, height = 5)
#plot_gene("NR1D1", "D") %>% ggsave("figures/clockgenes/NR1D1_D.png", ., width = 8, height = 5)
#plot_gene("NR1D2", "D") %>% ggsave("figures/clockgenes/NR1D2_D.png", ., width = 8, height = 5)
#plot_gene("PER1",  "D") %>% ggsave("figures/clockgenes/PER1_D.png", ., width = 8, height = 5)
#plot_gene("PER2",  "D") %>% ggsave("figures/clockgenes/PER2_D.png", ., width = 8, height = 5)
#plot_gene("PER3",  "D") %>% ggsave("figures/clockgenes/PER3_D.png", ., width = 8, height = 5)
#plot_gene("RORA",  "D") %>% ggsave("figures/clockgenes/RORA_D.png", ., width = 8, height = 5)
#
#plot_gene("ARNTL", "E") %>% ggsave("figures/clockgenes/ARNTL_E.png", ., width = 8, height = 5)
#plot_gene("ARNTL2","E") %>% ggsave("figures/clockgenes/ARNTL2_E.png", ., width = 8, height = 5)
#plot_gene("CLOCK", "E") %>% ggsave("figures/clockgenes/CLOCK_E.png", ., width = 8, height = 5)
#plot_gene("CRY1",  "E") %>% ggsave("figures/clockgenes/CRY1_E.png", ., width = 8, height = 5)
#plot_gene("CRY2",  "E") %>% ggsave("figures/clockgenes/CRY2_E.png", ., width = 8, height = 5)
#plot_gene("DBP",   "E") %>% ggsave("figures/clockgenes/DBP_E.png", ., width = 8, height = 5)
#plot_gene("NPAS2", "E") %>% ggsave("figures/clockgenes/NPAS2_E.png", ., width = 8, height = 5)
#plot_gene("NR1D1", "E") %>% ggsave("figures/clockgenes/NR1D1_E.png", ., width = 8, height = 5)
#plot_gene("NR1D2", "E") %>% ggsave("figures/clockgenes/NR1D2_E.png", ., width = 8, height = 5)
#plot_gene("PER1",  "E") %>% ggsave("figures/clockgenes/PER1_E.png", ., width = 8, height = 5)
#plot_gene("PER2",  "E") %>% ggsave("figures/clockgenes/PER2_E.png", ., width = 8, height = 5)
#plot_gene("PER3",  "E") %>% ggsave("figures/clockgenes/PER3_E.png", ., width = 8, height = 5)
#plot_gene("RORA",  "E") %>% ggsave("figures/clockgenes/RORA_E.png", ., width = 8, height = 5)

#--------------------------------
#--------------------------------


