suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr) )
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidytext))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(ggh4x))

setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis")

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

# Optional step to run analysis in parallel on multicore machines (here we use 4 threads)
cl <- makeCluster(4)
registerDoParallel(cl)
bpparam <- MulticoreParam(workers = 6)


#################################
#################################


# 1. LOAD GENE EXPRESSION DATA AS WELL AS RHYTHMIC RESULTS (ANALYSIS FIG 1)
# -------------------------------------------------------------------------
PCA_outliers <- "E32_P109" 
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 

# Gene expression data
yave <- readRDS("visualize/data/rawdata.rds")
ind <- which(colnames(yave) == PCA_outliers)    
yave <- yave[, -ind] 

geneExpr <- yave$E
genes <- yave$genes

results <- readRDS("visualize/data/results_populationrhy_internaltime.rds") %>% 
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)

results$rhythmic_in_D <- ifelse(results$diff_rhythmic==FALSE, TRUE,
                                ifelse(results$diff_rhythmic==TRUE & results$A_D > amp_cutoff, TRUE, FALSE)) 
results$rhythmic_in_E <- ifelse(results$diff_rhythmic==FALSE, TRUE, 
                                ifelse(results$diff_rhythmic==TRUE & results$A_E > amp_cutoff, TRUE, FALSE))

some_rhy <- dplyr::filter(results, rhythmic_in_D | rhythmic_in_E) #some_rhy as in fig1 analysis
geneExpr <- geneExpr %>% as.data.frame %>% filter(rownames(.) %in% some_rhy$ProbeName)

# Metadata
info_subjects <- read.csv("resources/info_subjects_short.csv") %>% dplyr::select(-X) #read info of subjects

info_exp <- data.frame(experiment2= geneExpr %>% colnames) %>% mutate(experiment=experiment2) %>%
  tidyr::separate(experiment2, c("tissuetime","subject"), sep = "_", convert=TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) 

experiment <- readRDS("visualize/data/experiment.rds") %>% full_join(info_subjects) # read sample details from column names
experiment <- experiment[-ind,] #remove outlier
info_exp <- info_exp %>% full_join(experiment) %>% 
  mutate(tissue = ifelse(tissue=="D", "dermis", "epidermis"),
         time = as.character(time),
         MSF_sc_dec = hms(MSF_sc),
         MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + second(MSF_sc_dec) / 360),2),
         inphase = cos(2*pi*as.numeric(internal_time)/24),
         outphase = sin(2*pi*as.numeric(internal_time)/24))



# 2. EXECUTE LINEAR MIXED MODEL WITH VARIANCE PARTITION
# -----------------------------------------------------
form <- ~ 1 + inphase + outphase + (1 + inphase + outphase|subject) + (1 + inphase + outphase||tissue)
fitList <- fitVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE, BPPARAM = bpparam) 

# Calculate variance in amplitude and phase from all rhythmic genes -> Error propagation
df_total = data.frame(ProbeName=NULL, Amp=NULL, var_A_layer=NULL, var_A_subject=NULL,
                      phase=NULL, var_phi_layer=NULL, var_phi_subject=NULL, 
                      magn=NULL, var_magn_layer=NULL, var_magn_subject=NULL)

for (i in 1:length(fitList)){
  # fixed effects (coefficients)
  m_i <- fixef(fitList[[i]])[1] %>% as.numeric 
  a_i <- fixef(fitList[[i]])[2] %>% as.numeric 
  b_i <- fixef(fitList[[i]])[3] %>% as.numeric
  
  # amplitude and phase (fixed effects) of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix Sigma (because there are >1 subject/tissue to do fits)
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList[[i]])))
  sigma_i_S <- sigma_i[1:3,1:3] #covariance_subject
  sigma_i_T <- sigma_i[4:6,4:6] #covariance_tissue
  
  # define Jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), # dA/dm, dA/da, dA/db
                   nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), 
                     nrow=1)
  Jac = matrix( c(0, a_i/A_i, b_i/A_i,
                  0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)) * (12/pi), (1/(1+(b_i/a_i)^2)) * (1/a_i) * (12/pi)), 
                nrow = 2, byrow = TRUE) 
  
  # determine variance in amplitude and phase, separating tissue and subject contributions
  var_A_subject <- Jac_A %*% sigma_i_S %*% t(Jac_A)
  var_A_layer <- Jac_A %*% sigma_i_T %*% t(Jac_A)
  var_phi_subject <- Jac_phi %*% sigma_i_S %*% t(Jac_phi)
  var_phi_layer <- Jac_phi %*% sigma_i_T %*% t(Jac_phi)
  
  var_S <- Jac %*% sigma_i_S %*% t(Jac) 
  var_T <- Jac %*% sigma_i_T %*% t(Jac)
  
  df_i <- data.frame(ProbeName=names(fitList)[i], 
                     Amp=A_i,     var_A_layer=var_A_layer,           var_A_subject=var_A_subject,
                     phase=phi_i, var_phi_layer=var_phi_layer,       var_phi_subject=var_phi_subject,
                     magn=m_i,    var_magn_layer=sigma_i_T[1,1], var_magn_subject=sigma_i_S[1,1])
  df_total <- rbind(df_total, df_i)
}

# Save results of linear mixed model + variances of amplitude, phase 
if (!file.exists("visualize/data/variance_rhythmic_parameters_full.csv")){
  write.csv(df_total %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol)),
            "visualize/data/variance_rhythmic_parameters_full.csv")
}
df_total <- read.csv("visualize/data/variance_rhythmic_parameters_full.csv") %>% dplyr::select(-X)
hist(df_total$Amp, breaks=100)
df_total %<>% filter(Amp>.15) # filter out genes with low amp_fit that result in high variability and "mask" variable genes
hist(df_total$Amp, breaks=100)


# 3. ARRANGE DATA FOR VISUALIZATION
# ---------------------------------
# Calculate sd and cv from variances, organize results for plots
variation_A   <- df_total %>% dplyr::select(ProbeName, Amp, var_A_subject, var_A_layer) %>% 
  gather(variable, variance, -ProbeName, -Amp) %>%
  mutate(variable = ifelse(variable=="var_A_subject", "A_S", "A_T"),
         sd = sqrt(variance),
         cv = sd/Amp) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol, EntrezID))
variation_phi <- df_total %>% dplyr::select(ProbeName, phase, var_phi_subject, var_phi_layer) %>%
  gather(variable, variance, -ProbeName, -phase) %>%
  mutate(variable = ifelse(variable=="var_phi_subject", "phi_S", "phi_T"),
         sd = sqrt(variance),
         cv = sd/24) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol, EntrezID))
variation_magn <- df_total %>% dplyr::select(ProbeName, magn, var_magn_subject, var_magn_layer) %>%
  gather(variable, variance, -ProbeName, -magn) %>%
  mutate(variable = ifelse(variable=="var_magn_subject", "magn_S", "magn_T"),
         sd = sqrt(variance),
         cv = sd/magn) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol, EntrezID))

variation_full <- rbind(variation_A %>% dplyr::rename(c("value_fit"="Amp")), 
                        variation_phi %>% dplyr::rename(c("value_fit"="phase")) %>%
                          mutate(value_fit = ifelse(value_fit < 0, value_fit + 24, value_fit))) %>% 
  rbind(variation_magn %>% dplyr::rename(c("value_fit"="magn"))) %>%
  tidyr::separate(variable, c("rhythmic_par","variable"), sep = "_", convert = TRUE) %>%
  mutate(rhythmic_par=ifelse(rhythmic_par == "A", "amplitude", 
                             ifelse(rhythmic_par=="phi", "phase", "magnitude")),
         variable=ifelse(variable == "S", "subject", "tissue")) 

variation_full$rhythmic_par <- factor(variation_full$rhythmic_par, levels=c("magnitude", "amplitude", "phase"))
variation_full$effect <- ifelse(variation_full$variable=="tissue", "layer", "subject")
variation_full$effect <- factor(variation_full$effect, levels=c("layer", "subject"))


# Correlations of fixed effects and variances
DR_genes    <- filter(some_rhy, diff_rhythmic)$Symbol
nonDR_genes <- filter(some_rhy, diff_rhythmic==FALSE)$Symbol
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP",
                 "NFIL3", "TEF", "BHLHE40", "BHLHE41")#, "HLF")

variation_full %<>% mutate(Symbol_it = paste0("italic('", Symbol, "')"))


# Fig2A: Distribution of sds -> How variable are magnitudes, amplitudes and phases across layers/subjects?
fig2A <- ggplot(variation_full %>% dplyr::select(rhythmic_par, sd, Symbol, effect)) +
  geom_violin(aes(y=effect, x=sd, fill=effect, color=effect), alpha=0.5, width=0.75) +
  facet_wrap(~rhythmic_par, scales='free_x') + 
  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  guides(fill=FALSE) +
  labs(x='standard deviation', y='') +
  theme( panel.spacing = unit(1., "lines"),
         legend.position = "top")


# 4. COEFFICIENT OF VARIATION -> Normalization of variability to the fixed effect: Where are the clock genes?
# -----------------------------------------------------------------------------------------------------------
# Fig2B: Correlations of cv_phi vs cv_amp -> where are the DR/non_DR genes located?
df <- variation_A %>% mutate(variable=ifelse(variable=="A_T", "tissue", "subject")) %>% 
  dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_amp"="cv")) %>% 
  full_join(variation_phi %>% mutate(variable=ifelse(variable=="phi_T", "tissue", "subject")) %>% 
              dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_phi"="cv"))) %>%
  full_join(variation_magn %>% mutate(variable=ifelse(variable=="magn_T", "tissue", "subject")) %>% 
              dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_magn"="cv"))) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))
df$variable <- factor(df$variable, levels=c("tissue", "subject")) %>% forcats::fct_recode(layer="tissue")

fig2B <- ggplot(df %>% dplyr::mutate(clock_gene = ifelse(Symbol %in% clock_genes, TRUE, FALSE))) +
  geom_point(aes(x=cv_amp, y=cv_phi), alpha=0.75, color="grey", size=.8) + 
  geom_point(data=filter(df, Symbol %in% DR_genes),      
             aes(x=cv_amp, y=cv_phi, color=variable), alpha=0.5,size=.8) + 
  facet_wrap(~variable, scales="free", nrow = 2, ncol=1) +
  geom_point(data=filter(df, Symbol %in% clock_genes & Symbol %in% nonDR_genes),      
             aes(x=cv_amp, y=cv_phi), alpha=1, shape=21, color="black", size=3, fill="grey20") +
  geom_point(data=filter(df, Symbol %in% clock_genes & Symbol %in% DR_genes),      
             aes(x=cv_amp, y=cv_phi, fill=variable), alpha=1, shape=21, color="black", size=3) +
  geom_text_repel(data=filter(df, Symbol %in% clock_genes), 
                  aes(x=cv_amp, y=cv_phi, label=Symbol_it), 
                  force        = 0.5,
                  nudge_x      = 2.0,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.2, 
                  parse=TRUE, size=2.8, segment.color="grey50") +
  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"), guide="none") +
  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c")) +
  labs(x='amplitude coefficient\nof variation', y='phase coefficient of variation') +
  theme(legend.position = "right",
        strip.text = element_blank(),
        panel.spacing = unit(10.5, "lines"),
        aspect.ratio = .64) +
  scale_y_continuous(limits=c(0,0.07)) + scale_x_continuous(limits=c(0,2.8), breaks=c(0, 0.5, 1, 1.5))


# 5. DETERMINE FRACTIONS OF VARIANCE -> How much of the variability is associated to subject/layer in each gene?
# --------------------------------------------------------------------------------------------------------------
amp <- variation_full %>% 
  mutate(vble = paste(rhythmic_par, variable, sep="_")) %>% 
  dplyr::filter(rhythmic_par=="amplitude") %>% 
  dplyr::select(Symbol, vble, variance) %>% 
  spread(vble, variance) %>%
  tibble::column_to_rownames("Symbol") %>%
  mutate(total=rowSums(.)) %>% tibble::rownames_to_column("Symbol") %>%
  gather(vble, variance, -Symbol, -total)
amp_frac <- amp[,c(1,3,4,2)] %>% mutate(fraction = variance/total) %>% dplyr::select(Symbol, vble, fraction) %>%
  spread(vble, fraction) %>% tibble::column_to_rownames("Symbol")

phi <- variation_full %>% 
  mutate(vble = paste(rhythmic_par, variable, sep="_")) %>% 
  dplyr::filter(rhythmic_par=="phase") %>% 
  dplyr::select(Symbol, vble, variance) %>% 
  spread(vble, variance) %>%
  tibble::column_to_rownames("Symbol") %>%
  mutate(total=rowSums(.)) %>% tibble::rownames_to_column("Symbol") %>%
  gather(vble, variance, -Symbol, -total)
phi_frac <- phi[,c(1,3,4,2)] %>% mutate(fraction = variance/total) %>% dplyr::select(Symbol, vble, fraction) %>%
  spread(vble, fraction) %>% tibble::column_to_rownames("Symbol")

magn <- variation_full %>% 
  mutate(vble = paste(rhythmic_par, variable, sep="_")) %>% 
  dplyr::filter(rhythmic_par=="magnitude") %>% 
  dplyr::select(Symbol, vble, variance) %>% 
  spread(vble, variance) %>%
  tibble::column_to_rownames("Symbol") %>%
  mutate(total=rowSums(.)) %>% tibble::rownames_to_column("Symbol") %>%
  gather(vble, variance, -Symbol, -total)
magn_frac <- magn[,c(1,3,4,2)] %>% mutate(fraction = variance/total) %>% dplyr::select(Symbol, vble, fraction) %>%
  spread(vble, fraction) %>% tibble::column_to_rownames("Symbol")

df_fraction_variance <- amp_frac %>% tibble::rownames_to_column("Symbol") %>%
  full_join(phi_frac %>% tibble::rownames_to_column("Symbol") ) %>%
  full_join(magn_frac %>% tibble::rownames_to_column("Symbol") ) %>% 
  gather(key, value, -Symbol) %>% 
  tidyr::separate(key, c("rhythm_par","variable"), convert = TRUE, sep = "_") 
df_fraction_variance$rhythm_par = factor(df_fraction_variance$rhythm_par, levels=c('magnitude','amplitude','phase'))

# Fig2C: How much fraction of variance in magn/amp/phase is attributed to subject/layer in our whole set of rhythmic genes?
fig2C_1 <- ggplot(df_fraction_variance %>% filter(variable=="subject")) +
  geom_histogram(aes(value, y=..density.., fill=variable), 
                 alpha=0.6, position='identity', colour="white", bins=50) +
  geom_density(aes(value, fill=variable, group=variable), alpha=0.3) +
  facet_wrap(~rhythm_par, ncol=1, scales="free") + 
  coord_cartesian(ylim=c(0, 3.5)) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none") +
  labs(x='fraction of variance\nexplained by subject', y='density')  +
  theme(panel.spacing = unit(.4, "lines")) + scale_x_continuous(breaks=c(0, 0.5, 1.0))
#https://ggrepel.slowkow.com/articles/examples.html#align-labels-on-the-left-or-right-edge-1

# Label in which bin to the clock genes fall (in fig2C_1)
# bin information of magnitude panel
info_hist_magn <- layer_data(fig2C_1, 1) %>% filter(PANEL==1)
only_CGs_magn <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="magnitude") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_magn <- c()
for (v in only_CGs_magn$value){
  m <- (info_hist_magn$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_magn <- append(idxs_magn, m2)
}
only_CGs_magn$hist_xvalue <- info_hist_magn[idxs_magn, "x"]
only_CGs_magn$hist_yvalue <- info_hist_magn[idxs_magn, "y"]
only_CGs_magn$rhythm_par <- "magnitude"

# bin information of amplitude panel
info_hist_amp <- layer_data(fig2C_1, 1) %>% filter(PANEL==2)
only_CGs_amp <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="amplitude") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_amp <- c()
for (v in only_CGs_amp$value){
  m <- (info_hist_amp$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_amp <- append(idxs_amp, m2)
}
only_CGs_amp$hist_xvalue <- info_hist_amp[idxs_amp, "x"]
only_CGs_amp$hist_yvalue <- info_hist_amp[idxs_amp, "y"]
only_CGs_amp$rhythm_par <- "amplitude"

# bin information of phase panel
info_hist_phi <- layer_data(fig2C_1, 1) %>% filter(PANEL==3)
only_CGs_phi <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="phase") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_phi <- c()
for (v in only_CGs_phi$value){
  m <- (info_hist_phi$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_phi <- append(idxs_phi, m2)
}
only_CGs_phi$hist_xvalue <- info_hist_phi[idxs_phi, "x"]
only_CGs_phi$hist_yvalue <- info_hist_phi[idxs_phi, "y"]
only_CGs_phi$rhythm_par <- "phase"

only_CGs <- rbind(only_CGs_magn, only_CGs_amp, only_CGs_phi)
only_CGs$rhythm_par = factor(only_CGs$rhythm_par, levels=c('magnitude','amplitude','phase'))
only_CGs$hist_yvalue <- ifelse(only_CGs$hist_yvalue>3.5, 3.5, only_CGs$hist_yvalue)

# add information of clock genes in bins from fig2C_1
fig2C <- fig2C_1 +   
  geom_text_repel(data=only_CGs, aes(x = hist_xvalue, y = hist_yvalue, label=Symbol_it), size=2.8, parse = TRUE, 
                  max.overlaps = Inf, box.padding = 1., force = 2, force_pull   = 0, # do not pull toward data points
                  nudge_y = 4.5, direction = "x", angle  = 90, hjust = 0.5, segment.size = 0.2, max.iter = 1e4, max.time = 1, 
                  segment.color="grey50")  


##############################
##############################


# 7. LINEAR MIXED MODEL IN DERMIS VS. EPIDERMIS
# ----------------------------------------------
form <- ~ 1 + inphase + outphase + (1 + inphase + outphase|subject)

geneExpr.D <- yave$E %>% as.data.frame() %>% dplyr::select(contains("D")) %>%   #rhy genes in dermis
  filter(rownames(.) %in% filter(some_rhy, rhythmic_in_D)$ProbeName)
geneExpr.E <- yave$E %>% as.data.frame() %>% dplyr::select(contains("E")) %>%   #rhy genes in epidermis
  filter(rownames(.) %in% filter(some_rhy, rhythmic_in_E)$ProbeName)

fitList.D <- fitVarPartModel(geneExpr.D, form, info_exp %>% filter(tissue=="dermis"), showWarnings=TRUE, BPPARAM = bpparam) 
fitList.E <- fitVarPartModel(geneExpr.E, form, info_exp %>% filter(tissue=="epidermis"), showWarnings=TRUE, BPPARAM = bpparam) 

# Calculate variance in amplitude and phase from all rhythmic genes in D -> Error propagation
df_total.D = data.frame(ProbeName=NULL, Amp=NULL, var_A_subject=NULL,
                        phase=NULL, var_phi_subject=NULL, 
                        magn=NULL, var_magn_subject=NULL)

for (i in 1:length(fitList.D)){
  # fixed effects (coefficients)
  m_i <- fixef(fitList.D[[i]])[1] %>% as.numeric 
  a_i <- fixef(fitList.D[[i]])[2] %>% as.numeric 
  b_i <- fixef(fitList.D[[i]])[3] %>% as.numeric
  
  # amplitude and phase of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix sigma 
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList.D[[i]])))
  
  # define jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), nrow=1)
  Jac = matrix( c(0, a_i/A_i, b_i/A_i,
                  0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), nrow = 2, byrow = TRUE)
  
  # determine variance in amplitude and phase, separating tissue and subject contributions
  var_A_subject <- Jac_A %*% sigma_i %*% t(Jac_A)
  var_phi_subject <- Jac_phi %*% sigma_i %*% t(Jac_phi)
  var_S <- Jac %*% sigma_i %*% t(Jac) 
  
  df_i <- data.frame(ProbeName=names(fitList)[i], 
                     Amp=A_i,     var_A_subject=var_A_subject,
                     phase=phi_i, var_phi_subject=var_phi_subject,
                     magn=m_i,    var_magn_subject=sigma_i[1,1])
  df_total.D <- rbind(df_total.D, df_i)
}
# Save results of linear mixed model + variances of amplitude, phase 
if (!file.exists("visualize/data/variance_rhythmic_parameters_dermis.csv")){
  write.csv(df_total.D %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol)),
            "visualize/data/variance_rhythmic_parameters_dermis.csv")
}
df_total.D <- read.csv("visualize/data/variance_rhythmic_parameters_dermis.csv") %>% dplyr::select(-X)
hist(df_total.D$Amp, breaks=100)
df_total.D %<>% filter(Amp>.15) # filter genes with low amp_fit that result in high variability and "mask" variable genes
hist(df_total.D$Amp, breaks=100)
df_total.D$tissue <- "dermis"

# Calculate variance in amplitude and phase from all rhythmic genes in E -> Error propagation
df_total.E = data.frame(ProbeName=NULL, Amp=NULL, var_A_subject=NULL,
                        phase=NULL, var_phi_subject=NULL, 
                        magn=NULL, var_magn_subject=NULL)

for (i in 1:length(fitList.E)){
  # fixed effects (coefficients)
  m_i <- fixef(fitList.E[[i]])[1] %>% as.numeric 
  a_i <- fixef(fitList.E[[i]])[2] %>% as.numeric 
  b_i <- fixef(fitList.E[[i]])[3] %>% as.numeric
  
  # amplitude and phase of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix sigma 
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList.E[[i]])))
  
  # define jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), nrow=1)
  Jac = matrix( c(0, a_i/A_i, b_i/A_i,
                  0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), nrow = 2, byrow = TRUE)
  
  # determine variance in amplitude and phase, separating tissue and subject contributions
  var_A_subject <- Jac_A %*% sigma_i %*% t(Jac_A)
  var_phi_subject <- Jac_phi %*% sigma_i %*% t(Jac_phi)
  var_S <- Jac %*% sigma_i %*% t(Jac) 
  
  df_i <- data.frame(ProbeName=names(fitList)[i], 
                     Amp=A_i,     var_A_subject=var_A_subject,
                     phase=phi_i, var_phi_subject=var_phi_subject,
                     magn=m_i,    var_magn_subject=sigma_i[1,1])
  df_total.E <- rbind(df_total.E, df_i)
}
# Save results of linear mixed model + variances of amplitude, phase 
if (!file.exists("visualize/data/variance_rhythmic_parameters_epidermis.csv")){
  write.csv(df_total.E %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol)),
            "visualize/data/variance_rhythmic_parameters_epidermis.csv")
}
df_total.E <- read.csv("visualize/data/variance_rhythmic_parameters_epidermis.csv") %>% dplyr::select(-X)
hist(df_total.E$Amp, breaks=100)
df_total.E %<>% filter(Amp>.15) # filter genes with low amp_fit that result in high variability and "mask" variable genes
hist(df_total.E$Amp, breaks=100)
df_total.E$tissue <- "epidermis"

df_total <- rbind(df_total.D, df_total.E)


# Calculate sd and cv from variances, organize results for plots
variation_A   <- df_total %>% dplyr::select(ProbeName, Amp, var_A_subject, tissue) %>% 
  rename(c("variance"="var_A_subject")) %>%
  mutate(sd = sqrt(variance),
         cv = sd/Amp,
         rhythmic_par = "amplitude") %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_phi <- df_total %>% dplyr::select(ProbeName, phase, var_phi_subject, tissue) %>%
  rename(c("variance"="var_phi_subject")) %>%
  mutate(sd = sqrt(variance),
         cv = sd/24,
         rhythmic_par = "phase") %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_magn <- df_total %>% dplyr::select(ProbeName, magn, var_magn_subject, tissue) %>%
  rename(c("variance"="var_magn_subject")) %>%
  mutate(sd = sqrt(variance),
         cv = sd/magn,
         rhythmic_par = "magnitude") %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))

variation_complete <- rbind(variation_A %>% dplyr::rename(c("value_fit"="Amp")), 
                        variation_phi %>% dplyr::rename(c("value_fit"="phase")) %>%
                          mutate(value_fit = ifelse(value_fit < 0, value_fit + 24, value_fit))) %>% 
  rbind(variation_magn %>% dplyr::rename(c("value_fit"="magn"))) 
variation_complete$rhythmic_par <- factor(variation_complete$rhythmic_par, levels=c("magnitude", "amplitude", "phase"))

fig2D <- ggplot(variation_complete, aes(x=sd, color=tissue)) +
  stat_density(size=1.5, geom="line", position="identity") +
  facet_wrap(~rhythmic_par, scales="free") +
  labs(x='standard deviation', y='') +
  theme( panel.spacing = unit(1., "lines"),
         legend.position = "top")


#############################
#############################

# Arrange plots in a grid
# ------------------------
fig2_1 <- plot_grid(NULL, NULL, fig2A, nrow=3, rel_heights=c(2,.05,1), labels = c("A", "", "B"))
fig2_3 <- plot_grid(fig2B, NULL, fig2D, nrow=3, rel_heights=c(2,.05, 1), labels = c("D", "E"))
fig2 <- plot_grid(fig2_1, NULL, fig2C, NULL, fig2_3, ncol=5, rel_widths=c(1.,0.03,.66,.03,1), labels = c("", "", "C", "", ""))
fig2 %>% ggsave('figures/fig2.pdf', ., width = 11, height = 7.5)