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

#--------------------------------
#--------------------------------

# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 
removal_outliers <- "weights" #PCA or weights
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R
save_fitdata <- TRUE

#--------------------------------
#--------------------------------

# 2. READ FILES
# -------------
# Read info of subjects, calculate mid sleeping time
setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis/")
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
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

#--------------------------------
#--------------------------------

# 4. CREATE DESIGN MATRIX AND FIT DATA
# ------------------------------------
if (save_fitdata == TRUE){
#if (!file.exists("visualize/data/results_fig2_P100_D.rds")){
  # design matrix
  design <- model.matrix(~ 0 + subject + tissue + subject:tissue:inphase + subject:tissue:outphase) #H0: rhythms are different across subjects and tissues
  
  # weights + fits
  fit <- limma::lmFit(yave, design, weights = wts) #!!!!
  fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  
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

# 6.1. CHECK MAXIMUM P_VALUES OF EACH GENE ACROSS SUBJECTS. CHECK RHYTHMIC GENES
# ------------------------------------------------------------------------------
# Separate results of subjects-dermis & subject-epidermis
# calculate first and second largest p values of each gene
# do multiple testing at end (although topTable already did it) to remove false positives

# A. Dermis
results_D <- results %>% select(ProbeName, Symbol, ENSEMBL, EntrezID, contains("_D"), contains("tissueD"), contains("phaseD"),
                                -contains("adj_P_Val"))

results_D$FirstLargest_pVal <- apply(select(results_D, starts_with("P_Val")), 1, 
                                     function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[1])
results_D$SecondLargest_pVal <- apply(select(results_D, starts_with("P_Val")), 1, 
                                 function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[2])
results_D$ThirdLargest_pVal <- apply(select(results_D, starts_with("P_Val")), 1, 
                                      function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[3])
results_D$FourthLargest_pVal <- apply(select(results_D, starts_with("P_Val")), 1, 
                                     function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[4])
results_D$FifthLargest_pVal <- apply(select(results_D, starts_with("P_Val")), 1, 
                                      function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[5])
results_D$FirstLargest_adjpVal <- p.adjust(results_D$FirstLargest_pVal, method ="BH")
results_D$SecondLargest_adjpVal <- p.adjust(results_D$SecondLargest_pVal, method ="BH")
results_D$ThirdLargest_adjpVal <- p.adjust(results_D$ThirdLargest_pVal, method ="BH")
results_D$FourthLargest_adjpVal <- p.adjust(results_D$FourthLargest_pVal, method ="BH")
results_D$FifthLargest_adjpVal <- p.adjust(results_D$FifthLargest_pVal, method ="BH")

## manual calculation of adjusted p values -> wiki: pval*total_no_tests/rank
results_D$rank_1pval <-  rank(results_D$FirstLargest_pVal)
results_D$rank_2pval <-  rank(results_D$SecondLargest_pVal)
results_D$rank_3pval <-  rank(results_D$ThirdLargest_pVal)
results_D$rank_4pval <-  rank(results_D$FourthLargest_pVal)
results_D$rank_5pval <-  rank(results_D$FifthLargest_pVal)
results_D$FirstLargest_adjpVal_manual  <- (results_D$FirstLargest_pVal  * dim(results_D)[1]) / results_D$rank_1pval
results_D$SecondLargest_adjpVal_manual <- (results_D$SecondLargest_pVal * dim(results_D)[1]) / results_D$rank_2pval
results_D$ThirdLargest_adjpVal_manual  <- (results_D$ThirdLargest_pVal  * dim(results_D)[1]) / results_D$rank_3pval
results_D$FourthLargest_adjpVal_manual <- (results_D$FourthLargest_pVal * dim(results_D)[1]) / results_D$rank_4pval
results_D$FifthLargest_adjpVal_manual  <- (results_D$FifthLargest_pVal  * dim(results_D)[1]) / results_D$rank_5pval

results_D$FirstLargest_adjpVal_manual[results_D$FirstLargest_adjpVal_manual > 1] <- 1
results_D$SecondLargest_adjpVal_manual[results_D$SecondLargest_adjpVal_manual > 1] <- 1
results_D$ThirdLargest_adjpVal_manual[results_D$ThirdLargest_adjpVal_manual > 1] <- 1
results_D$FourthLargest_adjpVal_manual[results_D$FourthLargest_adjpVal_manual > 1] <- 1
results_D$FifthLargest_adjpVal_manual[results_D$FifthLargest_adjpVal_manual > 1] <- 1

#input check of my manual adj p value calculation
ggplot(results_D) + geom_abline(aes(slope=1, intercept=0), color='grey') +
  geom_hline(yintercept=fdr_cutoff, color="red", linetype="dashed") +
  geom_vline(xintercept=fdr_cutoff, color="red", linetype="dashed") +
  geom_point(aes(FifthLargest_adjpVal, FifthLargest_adjpVal_manual)) + #xlim(0.95,1) + ylim(0.95,1)
  theme_bw() + theme(aspect.ratio=1) + xlab("p.adjust adj p_values") + ylab("manual adj p_values") + ggtitle("Dermis") 

# B. Epidermis
results_E <- results %>% select(ProbeName, Symbol, ENSEMBL, EntrezID, contains("_E"), contains("tissueE"), contains("phaseE"),
                                -contains("adj_P_Val")) 

results_E$FirstLargest_pVal <- apply(select(results_E, starts_with("P_Val")), 1, 
                                     function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[1]) #equivalent to a gene being rhy in all subjects
results_E$SecondLargest_pVal <- apply(select(results_E, starts_with("P_Val")), 1,
                                 function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[2]) #equivalent to a gene being rhy in all subjects except 1
results_E$ThirdLargest_pVal <- apply(select(results_E, starts_with("P_Val")), 1,
                                      function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[3]) #equivalent to a gene being rhy in all subj except 2
results_E$FourthLargest_pVal <- apply(select(results_E, starts_with("P_Val")), 1, 
                                      function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[4])
results_E$FifthLargest_pVal <- apply(select(results_E, starts_with("P_Val")), 1, 
                                      function(x) c(x[is.na(x)], sort(x, decreasing=TRUE))[5])
results_E$FirstLargest_adjpVal <- p.adjust(results_E$FirstLargest_pVal, method ="BH")
results_E$SecondLargest_adjpVal <- p.adjust(results_E$SecondLargest_pVal, method ="BH")
results_E$ThirdLargest_adjpVal <- p.adjust(results_E$ThirdLargest_pVal, method ="BH")
results_E$FourthLargest_adjpVal <- p.adjust(results_E$FourthLargest_pVal, method ="BH")
results_E$FifthLargest_adjpVal <- p.adjust(results_E$FifthLargest_pVal, method ="BH")

## manual calculation of adjusted p values -> wiki: pval*total_no_tests/rank
results_E$rank_1pval <-  rank(results_E$FirstLargest_pVal)
results_E$rank_2pval <-  rank(results_E$SecondLargest_pVal)
results_E$rank_3pval <-  rank(results_E$ThirdLargest_pVal)
results_E$rank_4pval <-  rank(results_E$FourthLargest_pVal)
results_E$rank_5pval <-  rank(results_E$FifthLargest_pVal)
results_E$FirstLargest_adjpVal_manual  <- (results_E$FirstLargest_pVal  * dim(results_E)[1]) / results_E$rank_1pval
results_E$SecondLargest_adjpVal_manual <- (results_E$SecondLargest_pVal * dim(results_E)[1]) / results_E$rank_2pval
results_E$ThirdLargest_adjpVal_manual  <- (results_E$ThirdLargest_pVal  * dim(results_E)[1]) / results_E$rank_3pval
results_E$FourthLargest_adjpVal_manual <- (results_E$FourthLargest_pVal * dim(results_E)[1]) / results_E$rank_4pval
results_E$FifthLargest_adjpVal_manual  <- (results_E$FifthLargest_pVal  * dim(results_E)[1]) / results_E$rank_5pval

results_E$FirstLargest_adjpVal_manual[results_E$FirstLargest_adjpVal_manual > 1] <- 1
results_E$SecondLargest_adjpVal_manual[results_E$SecondLargest_adjpVal_manual > 1] <- 1
results_E$ThirdLargest_adjpVal_manual[results_E$ThirdLargest_adjpVal_manual > 1] <- 1
results_E$FourthLargest_adjpVal_manual[results_E$FourthLargest_adjpVal_manual > 1] <- 1
results_E$FifthLargest_adjpVal_manual[results_E$FifthLargest_adjpVal_manual > 1] <- 1

#input check of my manual adj p value calculation
ggplot(results_E) + geom_abline(aes(slope=1, intercept=0), color='grey') +
  geom_hline(yintercept=fdr_cutoff, color="red", linetype="dashed") +
  geom_vline(xintercept=fdr_cutoff, color="red", linetype="dashed") +
  geom_point(aes(FifthLargest_adjpVal, FifthLargest_adjpVal_manual)) + #xlim(0.95,1) + ylim(0.95,1)
  theme_bw() + theme(aspect.ratio=1) + xlab("p.adjust adj p_values") + ylab("manual adj p_values") + ggtitle("Epidermis") 

# check rhythmic genes in each tissue
pval_filter <- "FifthLargest_adjpVal"

rhygenes_D <- results_D[ which(results_D[, pval_filter] < fdr_cutoff) ,] %>% #results_D %>% filter(pval_filter < fdr_cutoff) %>%
  mutate(Minimum_A_D = pmin(!!! select(., starts_with("A_D"))))
rhygenes_E <- results_E[ which(results_E[, pval_filter] < fdr_cutoff) ,] %>%
  mutate(Minimum_A_E = pmin(!!! select(., starts_with("A_E"))))
print(paste0('# of rhygenes in D = ', dim(rhygenes_D)[1]))
print(paste0('# of rhygenes in E = ', dim(rhygenes_E)[1]))


# 6.2. SORT GENES BY LARGEST P VALUE. CHECK RHYTHMIC GENES
# --------------------------------------------------------
results_D %<>% arrange(!!! pval_filter) 
results_E %<>% arrange(!!! pval_filter) 

#very few rhythmic genes if I want genes that are rhythmic in all subjects (except a few, 
#depending on the pval_filter) when I compare the xth adjusted p value to fdr_cutoff!

#--------------------------------
#--------------------------------

# 7. FOR THE RHYTHMIC GENES, CALCULATE PHASE_DIFF BETWEEN EACH SUBJECT AND A REFERENCE SUBJECT
# -------------------------------------------------------------------------------------------
ref_subject <- info_subjects %>% filter(MSF_sc == median(MSF_sc)) %$% subject [1] %>% as.character() #ref subject is the one with mid sleep time in middle (median)

temp_D <- rhygenes_D %>% select(contains("phaseD")) 
names(temp_D) <- gsub(x = names(temp_D), pattern = "phase", replacement = "delta_phase_")  
temp_D[1:ncol(temp_D)] <- temp_D[1:ncol(temp_D)]-temp_D[,1]
rhygenes_D %<>% cbind(temp_D)

temp_E <- rhygenes_E %>% select(contains("phaseE")) 
names(temp_E) <- gsub(x = names(temp_E), pattern = "phase", replacement = "delta_phase_")  
temp_E[1:ncol(temp_E)] <- temp_E[1:ncol(temp_E)]-temp_E[,1]
rhygenes_E %<>% cbind(temp_E)

# make sure delta_phase values go from -12 to +12
rhygenes_D[,49:59][rhygenes_D[,49:59] < -12] <- rhygenes_D[,49:59][rhygenes_D[,49:59] < -12] + 24
rhygenes_D[,49:59][rhygenes_D[,49:59] >  12] <- rhygenes_D[,49:59][rhygenes_D[,49:59] > 12] - 24
rhygenes_E[,49:59][rhygenes_E[,49:59] < -12] <- rhygenes_E[,49:59][rhygenes_E[,49:59] < -12] + 24
rhygenes_E[,49:59][rhygenes_E[,49:59] >  12] <- rhygenes_E[,49:59][rhygenes_E[,49:59] > 12] - 24 #before: 63:ncol(rhygenes_E)

rhygenes_D %>% select(contains("phase"))
rhygenes_E %>% select(contains("phase"))

# Calculate mean delta_phase of rhythmic genes in each subject
summary_rhygenes_D <- rhygenes_D %>% select(Symbol, contains("delta_phase")) %>% 
  gather(key, value, -Symbol) %>%
  group_by(key) %>% dplyr::summarise(value = mean(value)) 
summary_rhygenes_E <- rhygenes_E %>% select(Symbol, contains("delta_phase")) %>% 
  gather(key, value, -Symbol) %>%
  group_by(key) %>% dplyr::summarise(value = median(value)) 

#--------------------------------
#--------------------------------

# 8. PLOT MEAN DELTA_PHASE VS. MEAN SLEEP TIME
# --------------------------------------------
summary_rhygenes <- rbind(summary_rhygenes_D, summary_rhygenes_E) %>% as.data.frame() %>%
  mutate(tissue = ifelse(grepl("_D_", key), "dermis", "epidermis"))
unlisted_info <- unlist(strsplit(summary_rhygenes$key, "_"))
summary_rhygenes$tissue <- unlisted_info[seq(3, length(unlisted_info), 4)]
summary_rhygenes$subject <- unlisted_info[seq(4, length(unlisted_info), 4)]

summary_rhygenes <- summary_rhygenes %>% full_join(info_subjects)

lims <- as.POSIXct(strptime(c("1970-01-01 02:30:00","1970-01-01 07:15:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
fig2A <- ggplot(summary_rhygenes, aes(as.POSIXct(x = MSF_sc, format = "%H:%M"), y=value, color=tissue, shape=sex, label=subject)) +
  geom_point(size=3, alpha=0.8) +
  geom_text_repel() +
  scale_x_datetime(date_breaks = "1 hours", date_labels = "%H:%M", limits=lims) +
  xlab("Mid sleep time") + 
  ylab(paste0("median delta_phase of rhythmic genes (h)\n(with respect to ref. subject ", ref_subject, ")",
              "\n(proxy for chronotype)")) +
  theme_bw() + theme(aspect.ratio=1) + ggtitle(paste(paste0("Rhythm filtration: ", pval_filter),
                                                     paste0("rhygenes_D=", dim(rhygenes_D)[1], 
                                                            ", rhygenes_E=", dim(rhygenes_E)[1]), sep="\n"))
if (!file.exists("figures/fig2A.pdf")){ 
  fig2A %>% ggsave("figures/fig2A.pdf", ., width = 5, height = 5) 
}
#--------------------------------
#--------------------------------


## amplitude distributions of 100 most rhythmic genes in different subjects
#top100_rhygenes_extended <- top100_rhygenes %>%
#  dplyr::select(Symbol, matches("A_"), matches("P_Val"))# %>%
#  filter(adj_P_Val < fdr_cutoff) %>% #although from the top100_rhygenes, all have fdr < cutoff
#  gather(subject, amp_value, -Symbol, -adj_P_Val) %>%
#  #filter(amp_value > amp_cutoff) %>% #to plot only genes whose amplitude surpass the amp_threshold (OJO because if I filter based on A_min I might lose genes that have high Amp in one suject but not in other)
#  mutate(subject = gsub("A_", "", subject)) %>%
#  left_join(experiment %>% dplyr::select(sex, subject)) %>% distinct #un poco chapuza, no se por que acabo con 15400 rows
#
#plot_Amplitudes_100TopRhyGenes <- ggplot(data = top100_rhygenes_extended) +
#  geom_vline(aes(xintercept=amp_cutoff), color="black", linetype="dashed") +
#  geom_histogram(aes(amp_value, fill=sex), color="black", bins=30, size=0.1, alpha=0.5) + 
#  facet_wrap(~subject) +
#  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
#  theme_bw() + theme(aspect.ratio=1) + 
#  xlab("Amplitude") + ylab("number of genes\n(from 100 most rhythmic genes, all passed fdr cutoff)")
#
## calculate phase differences between each subject and a reference subject
#temp <- top100_rhygenes %>%
#  mutate(phase_diff_P100 = phase_P100 - phase_P100, #P100 is the median in terms of MSF_sc, thus I choose it as ref. subject
#         phase_diff_P102 = phase_P102 - phase_P100,
#         phase_diff_P103 = phase_P103 - phase_P100,
#         phase_diff_P106 = phase_P106 - phase_P100,
#         phase_diff_P107 = phase_P107 - phase_P100,
#         phase_diff_P108 = phase_P108 - phase_P100,
#         phase_diff_P109 = phase_P109 - phase_P100,
#         phase_diff_P111 = phase_P111 - phase_P100,
#         phase_diff_P113 = phase_P113 - phase_P100,
#         phase_diff_P114 = phase_P114 - phase_P100,
#         phase_diff_P115 = phase_P115 - phase_P100) %>%
#  dplyr::select(Symbol, matches("phase_diff_P"), adj_P_Val) %>%
#  filter(adj_P_Val < fdr_cutoff) %>% #although from the top100_rhygenes, all have fdr < cutoff
#  gather(subject, phase_diff_value, -Symbol, -adj_P_Val) %>%
#  #filter(amp_value > amp_cutoff) %>% #to plot only genes whose amplitude surpass the amp_threshold (OJO because if I filter based on A_min I might lose genes that have high Amp in one suject but not in other)
#  mutate(subject = gsub("phase_diff_", "", subject)) %>%
#  left_join(experiment %>% dplyr::select(sex, subject)) %>% distinct
#
#top100_rhygenes_extended %<>% full_join(temp)
#
## determine mean phase_diff of each subject (proxy for mean chronotype) -> convenio: owls have positive POE (see Springr book chapter)
#mean_chronotype <- top100_rhygenes_extended %>% 
#  group_by(subject) %>%
#  summarise(mean_chronotype = mean(phase_diff_value)) %>% as.data.frame()
#  
#
#
## Fig2C: How does number of rhy genes (in all patients) change with FDR
#fig2C <- ggplot(data = dplyr::filter(results, pmax(A_P100, A_P102, A_P103, A_P106, A_P107, A_P108, A_P109, A_P111, A_P113, 
#                                                   A_P114, A_P115) > amp_cutoff) %>% 
#                  dplyr::mutate(len = length(adj_P_Val))) + 
#  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
#  stat_ecdf(aes(x = adj_P_Val, len=len, y = ..y..*len), geom = "step") + 
#  scale_x_continuous(breaks = seq(0.00, 0.05, by=0.01)) +
#  coord_cartesian(xlim = c(0.001, 0.053), ylim = c(0,2000)) + theme_bw() +
#  xlab("False discovery rate") + ylab("number of rhythmic genes") 
#