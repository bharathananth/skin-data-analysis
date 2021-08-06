library(Rcpp)
library(RcppArmadillo)
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
library(ggpubr)
library(plotly) #interactive plot: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html
library(broom)

setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis/")
sourceCpp("circSVD.cpp")
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
results <- matrix(0, nrow=11578) %>% as.data.frame()
for (i in unique(subject)){
  for (j in unique(tissue)){
    # read results of each patient (TopTable was done in each patient individually, see notebook 26.03.2021)  
    results_i <- readRDS(paste0("visualize/data/results_fig2_", i, "_", j, ".rds")) %>% select(-AveExpr) %>%
      rename(!! paste0("F_",i,"_",j) := F) %>%
      rename(!! paste0("P_Value_",i,"_",j) := paste0("P_Value_",i)) %>%
      rename(!! paste0("adj_P_Val_",i,"_",j) := paste0("adj_P_Val_",i))
    # calculate amplitudes and phases from fits
    results_i[,paste0("A_", j, "_", i)] <- sqrt(results_i[,5]^2 + results_i[,6]^2)
    results_i[,paste0("phase", j, "_", i)] <- atan2(results_i[,6], results_i[,5])*12/pi #phase in h
    results_i[,paste0("rhythm_", i, "_",j)] <- complex(modulus=results_i[,10], argument=results_i[,11]*pi/12) #argument in rad!!
    
    # stack all columns together
    if (i == "P100" & j == "D"){
      results <- cbind(results, results_i) %>% select(-V1)
    } else{
      results <- full_join(results, results_i)
    }
  }
}

genes <- results %>% select(ProbeName, Symbol)

#phase_from_htorad <- results %>% select(contains("phaseD")) %>% gather(key, value) %>% mutate(value_rad = value/12*pi) %>% 
#  separate(key,c("junk", "subject"), sep="_", convert=TRUE) %>% select(-junk) %>% head
#phase_from_complex <- results %>% select(contains("rhythm")) %>% select(contains("D")) %>%
#  gather(key, value) %>% mutate(argument=Arg(value), modulus=Mod(value)) %>% 
#  separate(key,c("junk", "subject", "junk2"), sep="_", convert=TRUE) %>% select(-contains("junk")) %>% head
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
  mutate(rhythm_weighted = rhythm_complex / sqrt(rss), # I choose rss, easier to justify
         modulus_rhy = Mod(rhythm_complex),
         argument_rhy = Arg(rhythm_complex), #same value as argument_rhywt :)
         modulus_rhywt = Mod(rhythm_weighted)) 

rhy <- results %>% # How many genes with fdr<cutoff and mod(WEIGHTED, UNWEIGHTED) > amp_cutoff
  select(ProbeName, contains("P_Value"), Symbol) %>% gather(key, P_Value, -ProbeName, -Symbol) %>%
  #select(ProbeName, contains("adj_P")) %>% gather(key, P_Value, -ProbeName) %>%
  separate(key, c("junk1", "junk2","subject","tissue"), sep = "_", convert = TRUE) %>% 
  select(-contains("junk")) %>% 
  full_join(rhythms_weighted %>% select(-Symbol, -rse, -argument_rhy, -contains("modulus"))) %>%
  filter(p.adjust(P_Value, method="BH") < fdr_cutoff & Mod(rhythm_weighted) > amp_cutoff & Mod(rhythm_complex) > amp_cutoff) 
counts <- rhy %>% group_by(tissue) %>% count(tissue) %>% as.data.frame()

plot_Mod_fdr <- ggplot(data = results %>% 
                  #select(ProbeName, contains("adj_P")) %>% gather(key, P_Value, -ProbeName) %>%
                  #separate(key, c("junk1","junk2","junk3","subject","tissue"), sep = "_", convert = TRUE) %>% 
                  select(ProbeName, contains("P_Value")) %>% gather(key, P_Value, -ProbeName) %>%
                  separate(key, c("junk1","junk2","subject","tissue"), sep = "_", convert = TRUE) %>% 
                  select(-contains("junk")) %>% 
                  full_join(rhythms_weighted %>% select(-Symbol, -rse, -argument_rhy, -contains("modulus")))) +
  annotate("rect", xmin = 10^(-9), xmax = fdr_cutoff, ymin = amp_cutoff, ymax = 10^(2), alpha = .2, fill="yellow") +
  #geom_point(aes(x=P_Value, y=Mod(rhythm_weighted), color=tissue, alpha=0.4)) + 
  geom_point(aes(x=p.adjust(P_Value, method ="BH"), y=Mod(rhythm_weighted), color=tissue, alpha=0.4)) + 
  facet_wrap(~tissue) + 
  geom_hline(yintercept = amp_cutoff, linetype="dashed") + geom_vline(xintercept = fdr_cutoff, linetype="dashed") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_bw() + xlab("adjusted p value, BH") +
  ggtitle("Correlation metric vs fdr\nmetric == rhythm_complex/sqrt(rss)\ncareful because yellow box has some Mod(original)<cutoff")

plot_Mod_Fstat <- ggplot(data = results %>% 
                  select(ProbeName, contains("F")) %>% gather(key, F_stat, -ProbeName) %>%
                  separate(key, c("junk1","subject","tissue"), sep = "_", convert = TRUE) %>% 
                  select(-contains("junk")) %>% 
                  full_join(rhythms_weighted %>% select(-Symbol, -rse, -argument_rhy, -contains("modulus")))) +
  geom_point(aes(x=F_stat, y=Mod(rhythm_weighted), color=tissue, alpha=0.4)) + facet_wrap(~tissue) +
  geom_hline(yintercept = amp_cutoff, linetype="dashed") + #geom_vline(xintercept = fdr_cutoff, linetype="dashed")  + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_bw() + xlab("F statistic") +
  ggtitle("Correlation metric vs F stat\nmetric == rhythm_complex/sqrt(rss)\ncareful because some Mod(rhy_wt)>cutoff have Mod(rhy)<cutoff")

#--------------------------------
#--------------------------------

# 8. CIRCSVD
# ----------
sourceCpp("circSVD.cpp")
var_threshold <- .75
# Gene expression matrix X_t as complex numbers (~rhythmicity) NOT weighted by sqrt(rss)
XD_t <- rhythms_weighted %>% filter(tissue=="D" & Symbol %in% filter(mean_vars, tissue=="D" & mean_var > var_threshold)$Symbol) %>% 
  select(ProbeName, rhythm_weighted, subject, tissue) %>% 
  mutate(subjecttissue = paste(subject, tissue, sep="_")) %>% select(-tissue, -subject) %>% 
  spread(subjecttissue, rhythm_weighted) %>% column_to_rownames("ProbeName") 
XE_t <- rhythms_weighted %>% filter(tissue=="E" & Symbol %in% filter(mean_vars, tissue=="E" & mean_var > var_threshold)$Symbol) %>% 
  select(ProbeName, rhythm_weighted, subject, tissue) %>% 
  mutate(subjecttissue = paste(subject, tissue, sep="_")) %>% select(-tissue, -subject) %>% 
  spread(subjecttissue, rhythm_weighted) %>% column_to_rownames("ProbeName") 

# X matrix has cols=genes and rows=subjects
XD <- t(XD_t)
colnames(XD) <- rownames(XD_t)
rownames(XD) <- colnames(XD_t)

XE <- t(XE_t)
colnames(XE) <- rownames(XE_t)
rownames(XE) <- colnames(XE_t)

# compute circSVD
circSVD_D <- circSVD(XD)
circSVD_E <- circSVD(XE)

# QUESTIONS
# 1. who is the 'reference subject' when calculating Us? No one has a phase=0
# 2. which is the 'reference gene' when calculating Vs?

# U -> ~ chronotypes
UD <- circSVD_D$U #Mod of all u elements is the same
UE <- circSVD_E$U
UD_H <- t(Conj(UD))
UD_H * XD[,1]

rownames(UD) <- rownames(XD)
rownames(UE) <- rownames(XE)
colnames(UD) <- "U"
colnames(UE) <- "U"
UD %<>% as.data.frame() %>% tibble::rownames_to_column() %>%
  separate(rowname, c("subject", "tissue"), sep="_", convert=TRUE)
UE %<>% as.data.frame() %>% tibble::rownames_to_column() %>%
  separate(rowname, c("subject", "tissue"), sep="_", convert=TRUE)
U <- full_join(UD, UE) %>% arrange(desc(Mod(U))) %>% mutate(modulus = Mod(U), argument = Arg(U))

# V -> ~ genes that are important for chronotypes
VD <- circSVD_D$V
VE <- circSVD_E$V

rownames(VD) <- colnames(XD)
rownames(VE) <- colnames(XE)
colnames(VD) <- "V"
colnames(VE) <- "V"

VD %<>% as.data.frame() %>% tibble::rownames_to_column() %>% rename(c("ProbeName" = "rowname")) %>%
  mutate(tissue="D") %>% full_join(genes)
VE %<>% as.data.frame() %>% tibble::rownames_to_column() %>% rename(c("ProbeName" = "rowname")) %>%
  mutate(tissue="E") %>% full_join(genes)
V <- full_join(VD, VE) %>% arrange(desc(Mod(V))) %>% mutate(modulus = Mod(V), argument = Arg(V))

# Some plots
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
clock_genes <- genes %>% filter(Symbol %in% clock_genes)

plot_gene_rhyU_polar("PER2", "D")
plot_gene_rhyU_polar("KLF5", "D")


lims <- as.POSIXct(strptime(c("1970-01-01 02:30:00","1970-01-01 07:15:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
df_plot <- full_join(info_subjects, U) %>% mutate(argument=Arg(U), phase_h = argument*12/pi) %>%
  mutate(MSF_sc=as.POSIXct(MSF_sc, format = "%H:%M"))

#ggplot(data = df_plot, aes(as.POSIXct(x = MSF_sc, format = "%H:%M"), y=phase_h)) + 
#  geom_point(aes(shape=sex),size=4, alpha=0.7) + geom_text_repel(box.padding = 0.5, max.overlaps = Inf, aes(label=subject)) +
#  geom_smooth(method='lm') +
#  facet_wrap(~tissue) + theme_bw() +
#  xlab("Mid sleep time") + ylab("estimated chronotype (h)") +  
#  scale_x_datetime(date_breaks = "1 hours", date_labels = "%H:%M", limits=lims) +
#  ggtitle("Correlation between mid sleep time and chronotype")

ggscatter(df_plot, 
          x = "MSF_sc", y = "phase_h", shape = "sex", size = 3,
          color = "tissue", palette = "jco",
          add = "reg.line") + 
  geom_text_repel(box.padding = 1, max.overlaps = Inf, aes(label=subject)) + facet_wrap(~tissue) +
  stat_cor(label.y = 0) +
  stat_regline_equation(label.y = -1.5) +
  xlab("Mid sleep time") + ylab("estimated chronotype (h)") +  
  scale_x_datetime(date_breaks = "1 hours", date_labels = "%H:%M", limits=lims) +
  ggtitle("Correlation between mid sleep time and chronotype (analysis U)") #+
  #scale_y_continuous(limits=c(-7,10), breaks = c(-5, -2.5, 0, 2.5, 5))

chronotype_linreg <- df_plot %>% select(MSF_sc, phase_h, tissue) %>% group_by(tissue) %>%
  do(broom::tidy(lm(phase_h ~ MSF_sc, .))) %>% as.data.frame() 
chronotype_slopes <- chronotype_linreg %>% filter(term != "MSF_sc")

V %>% filter(modulus != 0 & tissue=="D") %>% dim #t=1e6: 11383, t=1e3: 11383, t=1e1: 219, t=50: 4674
V %>% filter(modulus != 0 & tissue=="E") %>% dim
V %>% filter(modulus != 0 & tissue=="D") %>% head
V %>% filter(modulus != 0 & tissue=="E") %>% head

hist(log10(V %>% filter(tissue=="D") %$% modulus))


# GO analysis on top V genes
# --------------------------
top_ngenes_GO <- 100
gD <- goana(yave[ yave$genes$ProbeName %in% V[which(V$tissue=="D" & V$modulus!=0),]$ProbeName ,]$genes$EntrezID, 
#gD <- goana(yave[ yave$genes$ProbeName %in% V[which(V$tissue=="D"),][1:top_ngenes_GO,]$ProbeName ,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
topGO(gD, n=20, truncate.term = "60") %>% as.data.frame() %>% filter(Ont=="BP") # signif. terms have p values ~ 10e-8
gE <- goana(yave[ yave$genes$ProbeName %in% V[which(V$tissue=="E" & V$modulus!=0),]$ProbeName ,]$genes$EntrezID, 
#gE <- goana(yave[ yave$genes$ProbeName %in% V[which(V$tissue=="E"),][1:top_ngenes_GO,]$ProbeName ,]$genes$EntrezID, 
            universe = yave$genes$EntrezID) 
topGO(gE, n=20, truncate.term = "60") %>% as.data.frame() %>% filter(Ont=="BP") # signif. terms have p values ~ 10e-8

gD %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N) %>% as.data.frame() %>% filter(Ont=="BP") %>% 
  ggplot(aes(x=hits, y=Term, colour=P.DE, size=DE)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count") +
  ggtitle(paste0("GO Term analysis on D (top ", top_ngenes_GO, " genes of circSVD analysis)")) +
  theme_bw()
gE %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N) %>% as.data.frame() %>% filter(Ont=="BP") %>% 
  ggplot(aes(x=hits, y=Term, colour=P.DE, size=DE)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count") +
  ggtitle(paste0("GO Term analysis on D (top ", top_ngenes_GO, " genes of circSVD analysis)")) +
  theme_bw()


# Bland Altman plot
# -----------------
slope_D <- chronotype_slopes %>% filter(tissue=="D") %$% estimate
slope_E <- chronotype_slopes %>% filter(tissue=="E") %$% estimate
df.ba <- df_plot %>% select(subject, sex, tissue, phase_h) %>% spread(tissue, phase_h)
df.ba %<>% mutate(D_corr = D-slope_D, E_corr = E-slope_E)
df.ba %<>% rowwise() %>% mutate(avg = mean(c(D_corr, E_corr)), diff = D_corr-E_corr) %>% as.data.frame()

ggplot(df.ba, aes(avg, diff)) + 
  geom_hline(yintercept=mean(df.ba$diff), linetype='dashed', color='grey55') + 
  geom_hline(yintercept=mean(df.ba$diff) + 1.96*sd(df.ba$diff), linetype='dotted', color='grey55') + 
  geom_hline(yintercept=mean(df.ba$diff) - 1.96*sd(df.ba$diff), linetype='dotted', color='grey55') + 
  #1.96 is based on the fact that 95% of area of a normal distrib is within 1.96 sd of the mean
  geom_point(aes(shape=sex), size=4) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, aes(label=subject)) + theme_bw() + coord_fixed() +
  xlab("average between chronotype in D and E") + 
  ylab("difference between chronotype in D and E (D-E)\n(U values have been 'corrected' with slopes)") +
  ggtitle("Bland-Altman plot (mean +- 1.96 sd plotted)")
  

# Determine the number of rhythmic genes with |v| != 0 as a function of the variance
# ----------------------------------------------------------------------------------
var_span <- seq(0.001, 0.50, by=0.001)
df <- data.frame(variance=NULL, total_no_genes_D=NULL, total_no_genes_E=NULL, v_no_genes_D=NULL, v_no_genes_E=NULL)

for (var_i in var_span) {
  XD_ti <- rhythms_weighted %>% filter(tissue=="D" & Symbol %in% filter(mean_vars, tissue=="D" & mean_var > var_i)$Symbol) %>% 
    select(ProbeName, rhythm_weighted, subject, tissue) %>% 
    mutate(subjecttissue = paste(subject, tissue, sep="_")) %>% select(-tissue, -subject) %>% 
    spread(subjecttissue, rhythm_weighted) %>% column_to_rownames("ProbeName") 
  XE_ti <- rhythms_weighted %>% filter(tissue=="E" & Symbol %in% filter(mean_vars, tissue=="E" & mean_var > var_i)$Symbol) %>% 
    select(ProbeName, rhythm_weighted, subject, tissue) %>% 
    mutate(subjecttissue = paste(subject, tissue, sep="_")) %>% select(-tissue, -subject) %>% 
    spread(subjecttissue, rhythm_weighted) %>% column_to_rownames("ProbeName") 
  
  # X matrix has cols=genes and rows=subjects
  XDi <- t(XD_ti)
  colnames(XDi) <- rownames(XD_ti)
  rownames(XDi) <- colnames(XD_ti)
  
  XEi <- t(XE_ti)
  colnames(XEi) <- rownames(XE_ti)
  rownames(XEi) <- colnames(XE_ti)
  
  # compute circSVD
  circSVD_D_i <- circSVD(XDi)
  circSVD_E_i <- circSVD(XEi)
  
  # V -> ~ genes that are important for chronotypes
  VD_i <- circSVD_D_i$V
  VE_i <- circSVD_E_i$V
  
  rownames(VD_i) <- colnames(XDi)
  rownames(VE_i) <- colnames(XEi)
  colnames(VD_i) <- "V"
  colnames(VE_i) <- "V"
  
  VD_i %<>% as.data.frame() %>% tibble::rownames_to_column() %>% rename(c("ProbeName" = "rowname")) %>%
    mutate(tissue="D") %>% full_join(genes)
  VE_i %<>% as.data.frame() %>% tibble::rownames_to_column() %>% rename(c("ProbeName" = "rowname")) %>%
    mutate(tissue="E") %>% full_join(genes)
  Vi <- full_join(VD_i, VE_i) %>% arrange(desc(Mod(V))) %>% mutate(modulus = Mod(V), argument = Arg(V))
  
  df_i <- data.frame(variance = var_i, 
                     total_no_genes_D = dim(Vi %>% filter(tissue=="D"))[1], 
                     total_no_genes_E = dim(Vi %>% filter(tissue=="E"))[1],
                     v_no_genes_D = dim(Vi %>% filter(modulus != 0 & tissue=="D"))[1], 
                     v_no_genes_E = dim(Vi %>% filter(modulus != 0 & tissue=="E"))[1])
  df <- rbind(df, df_i)
  
}
df2 <- df %>% mutate(percentage_no_genes_D = v_no_genes_D / total_no_genes_D * 100,
                     percentage_no_genes_E = v_no_genes_E / total_no_genes_E * 100) %>% gather(key, value, -variance) %>% 
  separate(key, c("filter", "junk1", "junk2", "tissue"), sep="_", convert=TRUE) %>% select(-contains("junk"))
df2$filter_f = factor(df2$filter, levels=c("total", "v", "percentage"))
filter_f.labs <- c("total number of genes", "genes with |v| != 0,\ndouble t = 1e1", "percentage of genes with |v|!=0\nfrom all genes that pass\nvariance filter")
names(filter_f.labs) <- c("total", "v", "percentage")

ggplot(df2, aes(variance, value)) + geom_line(aes(color=tissue)) + 
  theme_bw() + facet_wrap(~filter_f, scales="free", labeller = labeller(filter_f = filter_f.labs)) +
  xlab("Variance of time series") + ylab("number of genes") +  
  ggtitle("Number of genes that pass variance filter and that have |v| != 0")

########

# Heatmap of rhy in ONE subject -> Maybe better heat maps in fig1, where rhy analysis is not done on each subj individually????
subj_toplot <- "P100"
data_toplot <- fitted_values %>% tibble::rownames_to_column() %>% 
  rename(c("ProbeName"="rowname")) %>% inner_join(results %>% select(ProbeName, Symbol)) %>%
  select(-ProbeName) %>% 
  gather(key, value, -Symbol) %>% 
  separate(key, c("tissuetime", "subject"), sep="_", convert=TRUE)  %>%
  separate(tissuetime, into = c("tissue", "time"), "(?<=[A-Z])(?=[0-9])") %>%
  mutate(time = time %>% as.numeric()) %>%
  dplyr::group_by(Symbol, tissue, subject) %>%
  #dplyr::summarise(mean_value = mean(value)) %>% 
  dplyr::mutate(value = value - mean(value)) %>% as.data.frame()
  #filter(subject==subj_toplot)

phases <- results %>% select(contains("phaseD") | contains("phaseE") | Symbol) %>% gather(key, phase, -Symbol) %>% 
  separate(key, c("tissue", "subject"), sep="_", convert=TRUE)  %>%
  mutate(tissue=gsub("phase", "", tissue))
data_toplot %<>% full_join(phases) %>% arrange(phase)

ggplot(data_toplot %>% filter(Symbol %in% rhy$Symbol) %>% filter(subject==subj_toplot), aes(x=time, y=Symbol)) +
  geom_tile(aes(fill=value)) +
  facet_wrap(~tissue) +
  ggtitle(paste0("Expressed transcriptome -- cos fits in ", subj_toplot)) + 
  scale_fill_viridis_c(option="inferno")#scale_fill_distiller(palette = "RdPu")

superheat(data_toplot,
          scale = TRUE,
          left.label.text.size=3,
          bottom.label.text.size=3,
          bottom.label.size = .05,
          row.dendrogram = TRUE )
