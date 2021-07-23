library('doParallel')
library('dplyr')
library('ggplot2')
library('tidyr')
library('zeitzeiger') #https://github.com/hugheylab/zeitzeiger

library(limma)
setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis")
scale_colour_discrete <- function(...) {
  scale_colour_brewer(..., palette="Dark2")
}

# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 
removal_outliers <- "PCA" #PCA or weights
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R
save_fitdata <- FALSE

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

# variance of each time-series (for each gene in each subj+tissue)
variances <- yave$E %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  rename(c("ProbeName"="rowname")) %>% inner_join(results %>% select(ProbeName, Symbol)) %>%
  select(-ProbeName) %>% 
  gather(key, value, -Symbol) %>% separate(key, c("tissuetime", "subject"), sep="_", convert=TRUE)  %>%
  separate(tissuetime, into = c("tissue", "time"), "(?<=[A-Z])(?=[0-9])") %>%
  mutate(time = time %>% as.numeric()) %>% dplyr::group_by(Symbol, tissue, subject) %>%
  dplyr::summarise(sd = sd(value)) %>% as.data.frame() %>% mutate(variance=sd^2)
mean_vars <- variances %>% dplyr::group_by(Symbol, tissue) %>% dplyr::summarise(mean_var = mean(variance)) %>% as.data.frame()

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
# phase = 0 means that the respective gene peaks at 8AM (time of first sampling)

results_amp <- results %>% 
  #filter(adj_P_Val < fdr_cutoff) %>%
  #filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(ProbeName, Symbol, adj_P_Val, A_D, A_E, AveExpr) %>%
  gather(tissue, amp_value, -adj_P_Val, -Symbol, -AveExpr, -ProbeName) %>%
  filter(amp_value > amp_cutoff) %>%
  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
results_phase <- results %>% 
  #filter(adj_P_Val < fdr_cutoff) %>%
  filter(pmax(A_E, A_D) > amp_cutoff) %>%
  #filter(pmax(A_E, A_D) > amp_cutoff & adj_P_Val < fdr_cutoff) %>%
  select(ProbeName, Symbol, adj_P_Val, phaseD, phaseE) %>%
  gather(tissue, phase_value, -adj_P_Val, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))
results_passAmpcutoff <- inner_join(results_amp, results_phase)

rhy_results <- results_passAmpcutoff %>% filter(adj_P_Val < fdr_cutoff)


#--------------------------------
#--------------------------------
#--------------------------------
#--------------------------------
#--------------------------------


## 4. CREATE DESIGN MATRIX AND FIT DATA (different from Fig1.R)
## ------------------------------------
## design matrix
#design <- model.matrix(~ 0 + subject + tissue + subject:tissue:inphase + subject:tissue:outphase) #H0: rhythms are different across subjects and tissues
#
## weights + fits
#fit <- limma::lmFit(yave, design, weights = wts) #!!!!
#fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
#
#if (save_fitdata == TRUE){
#  # save results of fits -> topTable is done in each subject individually
#  for (i in unique(subject)){
#    for (j in unique(tissue)){
#      rhy_indices <- which(grepl(paste0(i,".*tissue",j,".*phase"), colnames(design)))
#      results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
#        set_colnames(gsub("\\.","_", colnames(.))) %>%
#        rename(!! paste0("P_Value_", i) := P_Value) %>%
#        rename(!! paste0("adj_P_Val_", i) := adj_P_Val) 
#      saveRDS(results, file = paste0("visualize/data/results_fig2_", i, "_", j, ".rds"))
#    }
#  }
#}
#
#

#--------------------------------
#--------------------------------

# 6. ZEITZEIGER (https://zeitzeiger.hugheylab.org/articles/introduction.html#load-the-necessary-packages-1)
# -------------

# Dataframe preparation for Zeitzeiger: 
# cols = genes, rows = all observations (across all subjects) -> separately for each tissue
xD <- yave$E %>% as.data.frame() %>% select(contains("D")) %>% 
  #tibble::rownames_to_column("ProbeName") %>% 
  #inner_join(rhy_results %>% filter(tissue=="dermis")) %>% 
  #select(-ProbeName, -adj_P_Val, -AveExpr, -tissue, -amp_value, -phase_value) %>% tibble::column_to_rownames("Symbol") %>%
  t()
rownames_ord <- paste0(sapply(strsplit(rownames(xD), split="_"), "[", 2), '_',
                       sapply(strsplit(rownames(xD), split="_"), "[", 1)) 
rownames(xD) <- rownames_ord
xD <- xD[order(rownames(xD)), ] %>% as.matrix() #ZeitZeiger takes input in this form


xE <- yave$E %>% as.data.frame() %>% select(contains("E")) %>% 
  #tibble::rownames_to_column("ProbeName") %>% 
  #inner_join(rhy_results %>% filter(tissue=="epidermis")) %>% 
  #select(-ProbeName, -adj_P_Val, -AveExpr, -tissue, -amp_value, -phase_value) %>% tibble::column_to_rownames("Symbol") %>%
  t()
rownames_ord <- paste0(sapply(strsplit(rownames(xE), split="_"), "[", 2), '_',
                       sapply(strsplit(rownames(xE), split="_"), "[", 1)) 
rownames(xE) <- rownames_ord
xE <- xE[order(rownames(xE)), ]%>% as.matrix()



# Crossvalidation
# ---------------
registerDoParallel(cores = 2)
sumabsv = c(1, 1.5, 2, 3)
nSpc = 1:5
nFolds = length(unique(subject))
nObs_D = dim(xD)[1]; nObs_E = dim(xE)[1]
nFeatures_D = dim(xD)[2]; nFeatures_E = dim(xE)[2]

#foldid_D = rep(1:nFolds, length.out = dim(xD)[1])
#foldid_E = rep(1:nFolds, length.out = dim(xD)[1])[-73] #E32_P109 is outlier so is renoved (quite chapuza)
foldid_D = rep(1:nFolds, each = length(unique(time)))
foldid_E = rep(1:nFolds, each = length(unique(time)))[-49] #E32_P109 is outlier so is renoved (quite chapuza)
#foldid_D = sample(rep(1:nFolds, length.out = dim(xD)[1]))

#time_D = str_sub(sapply(strsplit(rownames(xD), split="_"), "[", 1), start=2) %>% as.numeric()
#time_E = str_sub(sapply(strsplit(rownames(xE), split="_"), "[", 1), start=2) %>% as.numeric()
time_D = str_sub(sapply(strsplit(rownames(xD), split="_"), "[", 2), start=2) %>% as.numeric() 
time_E = str_sub(sapply(strsplit(rownames(xE), split="_"), "[", 2), start=2) %>% as.numeric() 
time_D = time_D/24; time_E = time_E/24

fitResultList_D = zeitzeigerFitCv(xD, time_D, foldid_D)
fitResultList_E = zeitzeigerFitCv(xE, time_E, foldid_E)

spcResultList_D = list()
spcResultList_E = list()
for (ii in 1:length(sumabsv)) {
  spcResultList_D[[ii]] = zeitzeigerSpcCv(fitResultList_D, sumabsv = sumabsv[ii])
  spcResultList_E[[ii]] = zeitzeigerSpcCv(fitResultList_E, sumabsv = sumabsv[ii])}

predResultList_D = list()
predResultList_E = list()
for (ii in 1:length(sumabsv)) {
  predResultList_D[[ii]] = zeitzeigerPredictCv(xD, time_D, foldid_D, spcResultList_D[[ii]], nSpc = nSpc)
  predResultList_E[[ii]] = zeitzeigerPredictCv(xE, time_E, foldid_E, spcResultList_E[[ii]], nSpc = nSpc)}


# Plot the error for each set of parameter values
# -----------------------------------------------
# Before plotting, we need to reorganize the output, making a data.frame with the information for each prediction.
timePredList_D = lapply(predResultList_D, function(a) a$timePred)
timePredList_E = lapply(predResultList_E, function(a) a$timePred)

cvResult_D = data.frame(do.call(rbind, timePredList_D),
                        timeObs = rep(time_D, length(sumabsv)),
                        sumabsv = rep(sumabsv, each = length(time_D)),
                        obs = rep(1:nObs_D, length(sumabsv)),
                        stringsAsFactors = FALSE)
cvResult_E = data.frame(do.call(rbind, timePredList_E),
                        timeObs = rep(time_E, length(sumabsv)),
                        sumabsv = rep(sumabsv, each = length(time_E)),
                        obs = rep(1:nObs_E, length(sumabsv)),
                        stringsAsFactors = FALSE)

cvResultGath_D = gather(cvResult_D, key = nSpc, value = timePred, -obs, -timeObs, -sumabsv)
cvResultGath_D$nSpc = as.integer(sapply(as.character(cvResultGath_D$nSpc),
                                        function(a) substr(a, 2, nchar(a))))
cvResultGath_D$sumabsv = factor(cvResultGath_D$sumabsv)
cvResultGath_D$timeError = getCircDiff(cvResultGath_D$timePred, cvResultGath_D$timeObs)

cvResultGath_E = gather(cvResult_E, key = nSpc, value = timePred, -obs, -timeObs, -sumabsv)
cvResultGath_E$nSpc = as.integer(sapply(as.character(cvResultGath_E$nSpc),
                                        function(a) substr(a, 2, nchar(a))))
cvResultGath_E$sumabsv = factor(cvResultGath_E$sumabsv)
cvResultGath_E$timeError = getCircDiff(cvResultGath_E$timePred, cvResultGath_E$timeObs)

# Now calculate the median absolute error for each set of parameter values
cvResultGathGroup_D = cvResultGath_D %>%
  group_by(sumabsv, nSpc) %>%
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()
cvResultGathGroup_D[,"tissue"] <- "dermis"
cvResultGathGroup_E = cvResultGath_E %>%
  group_by(sumabsv, nSpc) %>%
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()
cvResultGathGroup_E[,"tissue"] <- "epidermis"
cvResultGathGroup <- rbind(cvResultGathGroup_D, cvResultGathGroup_E)

ggplot(cvResultGathGroup) + facet_wrap(~tissue) +
  geom_point(aes(x = nSpc, y = medae, shape = sumabsv, color = sumabsv), size = 2) +
  labs(x = 'Number of SPCs', y = 'Median absolute error') +
  theme_bw() + theme(legend.position = c(0.7, 0.7))


# Train a model on the full dataset
# ---------------------------------
fitResultFinal_D = zeitzeigerFit(xD, time_D)
fitResultFinal_E = zeitzeigerFit(xE, time_E)
spcResultFinal_D = zeitzeigerSpc(fitResultFinal_D$xFitMean, fitResultFinal_D$xFitResid, 
                                 sumabsv = 3) #sumabsv=3 (maybe even 2?) gives least MAE in D
spcResultFinal_E = zeitzeigerSpc(fitResultFinal_E$xFitMean, fitResultFinal_E$xFitResid, 
                                 sumabsv = 4) #sumabsv=4 gives least MAE in E

dfVar_D = data.frame(spc = 1:length(spcResultFinal_D$d), propVar = spcResultFinal_D$d^2 / sum(spcResultFinal_D$d^2))
dfVar_E = data.frame(spc = 1:length(spcResultFinal_E$d), propVar = spcResultFinal_E$d^2 / sum(spcResultFinal_E$d^2))
dfVar_D[,"tissue"] = "dermis"; dfVar_E[,"tissue"] = "epidermis"
dfVar = rbind(dfVar_D, dfVar_E)

ggplot(dfVar) + facet_wrap(~tissue) +
  geom_point(aes(x = spc, y = propVar), size = 2, shape = 1) +
  scale_x_continuous(breaks = seq(1, 10)) +
  labs(x = 'SPC', y = 'Proportion of\nvariance explained') + theme_bw() #only first 3 SPCs explain variance


# Plot the behavior of SPCs over time
# -----------------------------------
zD = xD %*% spcResultFinal_D$v[, 1:3]
zE = xE %*% spcResultFinal_E$v[, 1:3]
colnames(zD) = c('SPC 1', 'SPC 2', 'SPC 3'); colnames(zE) = colnames(zD)

zGath_D = gather(data.frame(zD, obs = 1:nObs_D, Time = time_D, check.names = FALSE),
               key = SPC, value = Abundance, -obs, -Time)
zGath_E = gather(data.frame(zE, obs = 1:nObs_E, Time = time_E, check.names = FALSE),
                 key = SPC, value = Abundance, -obs, -Time)
zGath_D[,"tissue"] = "dermis"; zGath_E[,"tissue"] = "epidermis"
zGath = rbind(zGath_D, zGath_E)

ggplot(zGath) + 
  facet_wrap(tissue~SPC, scales = 'free_y') +
  geom_point(aes(x = Time, y = Abundance, color = tissue), size = 2, shape = 1) + theme_bw()

ggplot(zGath) + 
  facet_grid(SPC~tissue, scales = 'free_y') +
  geom_point(aes(x = Time, y = Abundance, color = tissue), size = 2, shape = 1) + theme_bw()


# Plot coefficients of the features (time-telling genes) for the SPCs
# -------------------------------------------------------------------
vD = data.frame(spcResultFinal_D$v[, 1:3])
vE = data.frame(spcResultFinal_E$v[, 1:3])
colnames(vD) = c('SPC 1', 'SPC 2', 'SPC 3'); colnames(vE) = colnames(vD)
vD = vD[apply(vD, 1, function(r) any(r != 0)), ]; vE = vE[apply(vE, 1, function(r) any(r != 0)), ]
vD[vD == 0] = NA; vE[vE == 0] = NA
vD = vD[do.call(order, vD), ]; vE = vE[do.call(order, vE), ]
vD$feature = rownames(vD); vE$feature = rownames(vE)
vD = inner_join(vD, yave$genes %>% as.data.frame() %>% select(Symbol) %>% mutate(feature = as.character(1:n())))
vE = inner_join(vE, yave$genes %>% as.data.frame() %>% select(Symbol) %>% mutate(feature = as.character(1:n())))

vGath_D = gather(vD, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vD$feature)),
                Symbol = factor(Symbol, levels = rev(vD$Symbol)), 
                tissue = "dermis")
vGath_E = gather(vE, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vE$feature)),
                Symbol = factor(Symbol, levels = rev(vE$Symbol)), 
                tissue = "epidermis")

pdf("figures/fig2_ZZ.pdf")
ggplot(vGath_D) + facet_wrap(~ spc, scales="free") +
  geom_bar(aes(x = Symbol2, y = Coefficient), stat = 'identity', fill="#1B9E77") +
  labs(x = 'Feature') + coord_flip() +
  theme_bw() +  ggtitle("ZeitZeiger on dermis") +
  theme(panel.spacing = unit(1.2, 'lines'),
        axis.text.y = element_text(face="italic", size=6),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold", size=10)) 

ggplot(vGath_E) + facet_wrap(~ spc, scales="free") +
  geom_bar(aes(x = Symbol, y = Coefficient), stat = 'identity', fill="#D95F02") +
  labs(x = 'Feature') + coord_flip() +
  theme_bw() + ggtitle("ZeitZeiger on epidermis") +
  theme(panel.spacing = unit(1.2, 'lines'),
      axis.text.y = element_text(face="italic", size=6),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face="bold"),
      plot.title = element_text(hjust = 0.5, face="bold", size=10)) 
dev.off()
