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
  rename(c("subject"="Subject", "sex"="Sex"))
info_subjects <- info_subjects_long %>% select(subject, sex, Light_condition, age, MSF_sc)

# Raw data (annotated, normalized, etc)
yave <- readRDS("visualize/data/rawdata.rds")


# Extract sample details from column names
experiment <- readRDS("visualize/data/experiment.rds")
experiment %<>% full_join(info_subjects)

tissue  <- factor(experiment$tissue)
time    <- experiment$time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

#--------------------------------
#--------------------------------

# 3. DO LINEAR MODEL FITS AND SAVE RESULTS OF FITS
inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

if (!file.exists("visualize/data/results_fig2_P100.rds")){ #Fig1: assuming rhythms are the same across subjects
  # design matrix
  #design <- model.matrix(~ 0 + subject + tissue + subject:inphase + subject:outphase) #H0: rhythms are different across subjects but are same in tissue
  design <- model.matrix(~ 0 + subject + tissue + subject:tissue:inphase + subject:tissue:outphase) #H0: rhythms are different across subjects and tissues
  
  # weights + fits
  wts <- limma::arrayWeights(yave, model.matrix(~tissue + subject)) # model.matrix: simpler model used as suggested in userguide. (avoid?)
  fit <- limma::lmFit(yave, design)#, weights = wts) 
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE) 
  
  # save results of fits
  for (i in unique(subject)){
    rhy_indices <- which(grepl(paste0(i,".*phase"), colnames(design)))
    results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
      set_colnames(gsub("\\.","_", colnames(.))) %>%
      rename(!! paste0("P_Value_", i) := P_Value) %>%
      rename(!! paste0("adj_P_Val_", i) := adj_P_Val) 
    saveRDS(results, file = paste0("visualize/data/results_fig2_", i, ".rds"))
  }
}

#--------------------------------
#--------------------------------

# 4. READ RESULTS AND PLAY WITH THEM
# ----------------------------------
i <- "P100" # from unique(subject)
resultsi <- readRDS(paste0("visualize/data/results_fig2_", i, ".rds")) %>% select(-AveExpr, -F)

results <- matrix(0, nrow= 11578) %>% as.data.frame()
for (i in unique(subject)){
  # read results of each patient (TopTable was done in each patient individually, see notebook 26.03.2021)
  results_i <- readRDS(paste0("visualize/data/results_fig2_", i, ".rds")) %>% select(-AveExpr, -F) 
  # calculate amplitudes and phases from fits
  results_i$A_D <- sqrt(results_i[,5]^2 + results_i[,7]^2)
  results_i$A_E <- sqrt(results_i[,6]^2 + results_i[,8]^2) 
  results_i$phaseD <- atan2(results_i[,7], results_i[,5])*12/pi
  results_i$phaseE <- atan2(results_i[,8], results_i[,6])*12/pi
  results_i %<>% 
    rename(!! paste0("A_D_", i) := A_D) %>%
    rename(!! paste0("A_E_", i) := A_E) %>%
    rename(!! paste0("phaseD_", i) := phaseD) %>%
    rename(!! paste0("phaseE_", i) := phaseE)
  
  if (i == "P100"){
    results <- cbind(results, results_i) %>% select(-V1)
  } else{
    results <- full_join(results, results_i)
  }
}

## VOY POR AQUI!!! lo de abajo es antiguo

# 100 most rhythmic genes across subjects
top100_rhygenes <- results[1:100,] %>%
  mutate(A_min = pmin(A_P100, A_P102, A_P103, A_P106, A_P107, A_P108, A_P109, A_P111, A_P113, A_P114, A_P115)) 
top100_rhygenes_noAmp_cutoff <- top100_rhygenes %>% filter(A_min < amp_cutoff) #genes in which AT LEAST ONE amplitude across subjects is < cutoff
#top100_rhygenes_org <- top100_rhygenes %>%
#  dplyr::select(Symbol, matches("A_P"), matches("phase_P"))

# amplitude distributions of 100 most rhythmic genes in different subjects
top100_rhygenes_extended <- top100_rhygenes %>%
  dplyr::select(Symbol, matches("A_P"), adj_P_Val) %>%
  filter(adj_P_Val < fdr_cutoff) %>% #although from the top100_rhygenes, all have fdr < cutoff
  gather(subject, amp_value, -Symbol, -adj_P_Val) %>%
  #filter(amp_value > amp_cutoff) %>% #to plot only genes whose amplitude surpass the amp_threshold (OJO because if I filter based on A_min I might lose genes that have high Amp in one suject but not in other)
  mutate(subject = gsub("A_", "", subject)) %>%
  left_join(experiment %>% dplyr::select(sex, subject)) %>% distinct #un poco chapuza, no se por que acabo con 15400 rows

plot_Amplitudes_100TopRhyGenes <- ggplot(data = top100_rhygenes_extended) +
  geom_vline(aes(xintercept=amp_cutoff), color="black", linetype="dashed") +
  geom_histogram(aes(amp_value, fill=sex), color="black", bins=30, size=0.1, alpha=0.5) + 
  facet_wrap(~subject) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw() + theme(aspect.ratio=1) + 
  xlab("Amplitude") + ylab("number of genes\n(from 100 most rhythmic genes, all passed fdr cutoff)")

# calculate phase differences between each subject and a reference subject
temp <- top100_rhygenes %>%
  mutate(phase_diff_P100 = phase_P100 - phase_P100, #P100 is the median in terms of MSF_sc, thus I choose it as ref. subject
         phase_diff_P102 = phase_P102 - phase_P100,
         phase_diff_P103 = phase_P103 - phase_P100,
         phase_diff_P106 = phase_P106 - phase_P100,
         phase_diff_P107 = phase_P107 - phase_P100,
         phase_diff_P108 = phase_P108 - phase_P100,
         phase_diff_P109 = phase_P109 - phase_P100,
         phase_diff_P111 = phase_P111 - phase_P100,
         phase_diff_P113 = phase_P113 - phase_P100,
         phase_diff_P114 = phase_P114 - phase_P100,
         phase_diff_P115 = phase_P115 - phase_P100) %>%
  dplyr::select(Symbol, matches("phase_diff_P"), adj_P_Val) %>%
  filter(adj_P_Val < fdr_cutoff) %>% #although from the top100_rhygenes, all have fdr < cutoff
  gather(subject, phase_diff_value, -Symbol, -adj_P_Val) %>%
  #filter(amp_value > amp_cutoff) %>% #to plot only genes whose amplitude surpass the amp_threshold (OJO because if I filter based on A_min I might lose genes that have high Amp in one suject but not in other)
  mutate(subject = gsub("phase_diff_", "", subject)) %>%
  left_join(experiment %>% dplyr::select(sex, subject)) %>% distinct

top100_rhygenes_extended %<>% full_join(temp)

# determine mean phase_diff of each subject (proxy for mean chronotype) -> convenio: owls have positive POE (see Springr book chapter)
mean_chronotype <- top100_rhygenes_extended %>% 
  group_by(subject) %>%
  summarise(mean_chronotype = mean(phase_diff_value)) %>% as.data.frame()
  


# Fig2C: How does number of rhy genes (in all patients) change with FDR
fig2C <- ggplot(data = dplyr::filter(results, pmax(A_P100, A_P102, A_P103, A_P106, A_P107, A_P108, A_P109, A_P111, A_P113, 
                                                   A_P114, A_P115) > amp_cutoff) %>% 
                  dplyr::mutate(len = length(adj_P_Val))) + 
  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
  stat_ecdf(aes(x = adj_P_Val, len=len, y = ..y..*len), geom = "step") + 
  scale_x_continuous(breaks = seq(0.00, 0.05, by=0.01)) +
  coord_cartesian(xlim = c(0.001, 0.053), ylim = c(0,2000)) + theme_bw() +
  xlab("False discovery rate") + ylab("number of rhythmic genes") 
