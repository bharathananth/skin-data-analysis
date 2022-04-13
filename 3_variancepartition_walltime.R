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
suppressPackageStartupMessages(library(clusterProfiler))
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



# 2. EXECUTE VARIANCE PARTITION MODEL
# -----------------------------------
#form <- ~ (1|subject) + (1|tissue) + (1|time) + (1|time:tissue) + (1|time:subject) 
#        #time here is wall time (continuous vble [internal_time] can't be modeled as random)
form <- ~ 1 + inphase + outphase + (1 + inphase + outphase|subject) + (1 + inphase + outphase||tissue)

#varPart <- fitExtractVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE) 
#varPart <- varPart %>% as.data.frame() %>%  
#  mutate(ProbeName=rownames(.)) %>% inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>%
#  tibble::column_to_rownames("Symbol")
#vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID))
#    #sortCols: orders columns such that the first one is the one with largest (mean) variation, second is second largest, etc
#
#if (!file.exists("visualize/data/variancePartition_full.csv")){
#  write.csv(vp %>% tibble::rownames_to_column("Symbol"), "visualize/data/variancePartition_full.csv")
#}
# Save 50 genes with least variability across subjects -> supplementary Table 3
#if (!file.exists("figures/supp_table3.csv")){
#  write.csv(varPart %>% dplyr::arrange(subject) %>% head(50) %>% tibble::rownames_to_column("Symbol"), 
#            "figures/supp_table3.csv")
#}

fitList <- fitVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE) 


# CALCULATE VARIANCE IN AMPLITUDE AND PHASE FROM ALL RHYTHMIC GENES -> ERROR PROPAGATION
df_total = data.frame(ProbeName=NULL, Amp=NULL, var_A_layer=NULL, var_A_subject=NULL,
                      phase=NULL, var_phi_layer=NULL, var_phi_subject=NULL, 
                      magn=NULL, var_magn_layer=NULL, var_magn_subject=NULL)

for (i in 1:length(fitList)){
  
  # fixed effects (coefficients)
  m_i <- fixef(fitList[[i]])[1] %>% as.numeric 
  a_i <- fixef(fitList[[i]])[2] %>% as.numeric 
  b_i <- fixef(fitList[[i]])[3] %>% as.numeric
  
  # amplitude and phase of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix sigma (because there are >1 subject/tissue to do fits)
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList[[i]])))
  sigma_i_S <- sigma_i[1:3,1:3] #covariance_subject
  sigma_i_T <- sigma_i[4:6,4:6] #covariance_tissue
  
  # define jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), # dA/dm, dA/da, dA/db
                   nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), 
                     nrow=1)
  Jac = matrix( c(0, a_i/A_i, b_i/A_i,
                  0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), # dA/dm, dA/da, dA/db
                nrow = 2, byrow = TRUE) #!!!!!!!
  
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

# Save results of variance Partition
if (!file.exists("visualize/data/variance_rhythmic_parameters.csv")){
  write.csv(df_total %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol)),
            "visualize/data/variance_rhythmic_parameters.csv")
}
df_total <- read.csv("visualize/data/variance_rhythmic_parameters.csv") %>% dplyr::select(-X)
hist(df_total$Amp, breaks=100)
df_total %<>% filter(Amp>.15) # filter genes with low amp_fit that result in high variability and "mask" variable genes
hist(df_total$Amp, breaks=100)

# PLOT DISTRIBUTIONS OF RANDOM EFFECTS -> HOW VARIABLE ARE MAGNITUDES, AMPLITUDES AND PHASES ACROSS LAYERS/SUBJECTS?
#ranef_total_subj <- data.frame('(Intercept)'=NULL, 'inphase'=NULL, 'outphase'=NULL, 'gene'=NULL)
#ranef_total_lay  <- data.frame('(Intercept)'=NULL, 'inphase'=NULL, 'outphase'=NULL, 'gene'=NULL)
#for (g in c(1:length(fitList))){
#  
#  #ranef_g_subj <- ranef(fitList[[g]])$subject %>% as.data.frame
#  #ranef_g_lay <- ranef(fitList[[g]])$tissue%>% as.data.frame
#  
#  ranef_total_subj <- rbind(ranef_total_subj, ranef_g_subj)
#  ranef_total_lay <- rbind(ranef_total_lay, ranef_g_lay)
#}


ranef_total_subj <- data.frame('(Intercept)'=NULL, 'inphase'=NULL, 'outphase'=NULL, 'gene'=NULL)
for (g in c(1:length(fitList))){
  #new_gene <- ranef_g_subj[r,]
  new_gene <- ranef(fitList[[g]])$subject %>% as.data.frame
  new_gene$gene <- paste0('gene', g)
  ranef_total_subj <- rbind(ranef_total_subj, new_gene)
}
ranef_total_lay  <- data.frame('(Intercept)'=NULL, 'inphase'=NULL, 'outphase'=NULL, 'gene'=NULL)
for (g in c(1:length(fitList))){
  #new_gene <- ranef_g_lay[r,]
  new_gene <- ranef(fitList[[g]])$tissue %>% as.data.frame
  new_gene$gene <- paste0('gene', g)
  ranef_total_lay <- rbind(ranef_total_lay, new_gene)
}

ranef_total_subj$effect <- 'subject'
ranef_total_subj$factor <- rep(unique(info_exp$subject), times=length(fitList))
ranef_total_lay$effect <- 'layer'
ranef_total_lay$factor  <- rep(unique(info_exp$tissue), times=length(fitList))

ranef_total <- full_join(ranef_total_lay, ranef_total_subj) %>% dplyr::rename(c('magnitude'='(Intercept)'))
ranef_total$amplitude <- sqrt(ranef_total$inphase^2 + ranef_total$outphase^2) #log2
ranef_total$phase <- atan2(ranef_total$outphase, ranef_total$inphase)*12/pi
ranef_total$phase <- ifelse(ranef_total$phase + 8 < 0, ranef_total$phase + 8 + 24, ranef_total$phase + 8)
ranef_total_gath <- ranef_total %>% select(-amplitude, -phase) %>% gather(rhythmic_par, value, -gene, -effect, -factor)
#ranef_total_gath <- ranef_total %>% select(-inphase, -outphase) %>% gather(rhythmic_par, value, -gene, -effect, -factor)
ranef_total_gath$rhythmic_par = factor(ranef_total_gath$rhythmic_par, levels=c('magnitude','inphase','outphase'))

#p1 <- ggplot(ranef_total_gath %>% filter(rhythmic_par=="magnitude")) + 
#  geom_jitter(aes(y=effect, x=value, group=effect, color=effect, fill=effect), alpha=0.3, size=.5) +
#  geom_boxplot(aes(y=effect, x=value, fill=effect), alpha=0.1, width=0.33, fill="white", outlier.size = 1) +
#  #facet_wrap(~rhythmic_par, scales='free') + 
#  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  guides(color=FALSE, fill=FALSE) + 
#  labs(x='variability in m', y='') 
#p1 <- ggExtra::ggMarginal(p1, groupFill = TRUE, type="density", margins = "x", xparams = list(size = .3), size=1.5)
#
#p2 <- ggplot(ranef_total_gath %>% filter(rhythmic_par=="inphase")) + 
#  geom_jitter(aes(y=effect, x=value, group=effect, color=effect, fill=effect), alpha=0.1, size=.5) +
#  geom_boxplot(aes(y=effect, x=value, fill=effect), alpha=0.1, outlier.size = 1, width=0.33, fill="white") +
#  #facet_wrap(~rhythmic_par, scales='free') + 
#  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  guides(color=FALSE, fill=FALSE)  + labs(x='variability in a', y='')
#p2 <- ggExtra::ggMarginal(p2, groupFill = TRUE, type="density", margins = "x", xparams = list(size = .3), size=1.5)
#
#p3 <- ggplot(ranef_total_gath %>% filter(rhythmic_par=="outphase")) + 
#  geom_jitter(aes(y=effect, x=value, group=effect, color=effect, fill=effect), alpha=0.1, size=.5) +
#  geom_boxplot(aes(y=effect, x=value, fill=effect), alpha=0.1, width=0.33, fill="white", outlier.size = 1) +
#  #facet_wrap(~rhythmic_par, scales='free') + 
#  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
#  guides(color=FALSE, fill=FALSE) +
#  #scale_x_continuous(breaks=c(0,6,12,18,24)) + 
#  labs(x='variability in b', y='')
#p3 <- ggExtra::ggMarginal(p3, groupFill = TRUE, type="density", margins = "x", xparams = list(size = .3),  size=1.5)
#
#fig2A <- ggpubr::ggarrange(p1, NULL, p2, NULL, p3, nrow=1, ncol=5, widths = c(1, -0.5, 1, -0.5, 1))

ranef_total_gath$effect <- factor(ranef_total_gath$effect, levels=c("layer", "subject"))
fig2A <- ggplot(ranef_total_gath) + 
  geom_jitter(aes(y=effect, x=value, group=effect, color=effect, fill=effect), alpha=0.1, size=.5) +
  geom_boxplot(aes(y=effect, x=value, fill=effect), alpha=0.1, width=0.33, fill="white", outlier.size = 1) +
  facet_wrap(~rhythmic_par, scales='free') + 
  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  guides(color=FALSE, fill=FALSE) +
  #scale_x_continuous(breaks=c(0,6,12,18,24)) + 
  labs(x='variability', y='') +
  ggtitle('better distribution of sds?') +
  theme( panel.spacing = unit(2.5, "lines"))


# Save 50 genes with least magnitude variability (variance) across subjects -> supplementary Table 3
if (!file.exists("figures/supp_table3.csv")){
  write.csv(df_total %>% dplyr::arrange(var_magn_subject) %>% head(50), #%>% tibble::column_to_rownames("Symbol") %>% dplyr::select(-ProbeName),
            "figures/supp_table3.csv")
}

# Save 50 genes with largest magnitude variability (variance) across subjects -> supplementary Table 4
if (!file.exists("figures/supp_table4.csv")){
  write.csv(df_total %>% dplyr::arrange(desc(var_magn_subject)) %>% head(50),# %>% tibble::column_to_rownames("Symbol") %>% dplyr::select(-ProbeName),
            "figures/supp_table4.csv")
}



# Calculate sd and cv from variances, organize results for plots
variation_A   <- df_total %>% dplyr::select(ProbeName, Amp, var_A_subject, var_A_layer) %>% 
  gather(variable, variance, -ProbeName, -Amp) %>%
  mutate(variable = ifelse(variable=="var_A_subject", "A_S", "A_T"),
         sd = sqrt(variance),
         cv = sd/Amp) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_phi <- df_total %>% dplyr::select(ProbeName, phase, var_phi_subject, var_phi_layer) %>%
  gather(variable, variance, -ProbeName, -phase) %>%
  mutate(variable = ifelse(variable=="var_phi_subject", "phi_S", "phi_T"),
         sd = sqrt(variance),
         cv = sd/24) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_magn <- df_total %>% dplyr::select(ProbeName, magn, var_magn_subject, var_magn_layer) %>%
  gather(variable, variance, -ProbeName, -magn) %>%
  mutate(variable = ifelse(variable=="var_magn_subject", "magn_S", "magn_T"),
         sd = sqrt(variance),
         cv = sd/magn) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))

variation_full <- rbind(variation_A %>% dplyr::rename(c("value_fit"="Amp")), 
                        variation_phi %>% dplyr::rename(c("value_fit"="phase")) %>%
                          mutate(value_fit = value_fit + 8,
                                 value_fit = ifelse(value_fit < 0, value_fit + 24, value_fit))) %>% 
  rbind(variation_magn %>% dplyr::rename(c("value_fit"="magn"))) %>%
  tidyr::separate(variable, c("rhythmic_par","variable"), sep = "_", convert = TRUE) %>%
  mutate(rhythmic_par=ifelse(rhythmic_par == "A", "amplitude", 
                             ifelse(rhythmic_par=="phi", "phase", "magnitude")),
         variable=ifelse(variable == "S", "subject", "tissue")) 

variation_full$rhythmic_par <- factor(variation_full$rhythmic_par, levels=c("magnitude", "amplitude", "phase"))
variation_full$effect <- ifelse(variation_full$variable=="tissue", "layer", "subject")
variation_full$effect <- factor(variation_full$effect, levels=c("layer", "subject"))
fig2A_2 <- ggplot(variation_full %>% dplyr::select(rhythmic_par, sd, Symbol, effect)) +
  geom_jitter(aes(y=effect, x=sd, group=effect, color=effect, fill=effect), alpha=0.1, size=.5) +
  geom_boxplot(aes(y=effect, x=sd, fill=effect), alpha=0.1, width=0.33, fill="white", outlier.size = 1) +
  facet_wrap(~rhythmic_par, scales='free') + 
  scale_fill_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  scale_color_manual(values = c("layer" = "#d1495b", "subject" = "#00798c"))  +
  guides(color=FALSE, fill=FALSE) +
  #scale_x_continuous(breaks=c(0,6,12,18,24)) + 
  labs(x='standard deviation', y='') +
  theme( panel.spacing = unit(1.6, "lines"))

#plotVarPart(sortCols(df_total %>% dplyr::select(Symbol, var_A_layer, var_A_subject, var_phi_layer, 
#                                                var_phi_subject,  var_magn_layer, var_magn_subject) %>% 
#                       tibble::column_to_rownames("Symbol"))) + 
#  ggtitle('WHAT ABOUT RESIDUALS, x axis is wrong\nI dont understand this plot')

#df_forviolin <- variation_full %>% dplyr::mutate(vble = paste0(rhythmic_par, "_", variable)) %>% 
#  dplyr::select(Symbol, vble, cv)
#df_forviolin$vble = factor(df_forviolin$vble, levels=c('amplitude_subject','amplitude_tissue',
#                                                       'magnitude_tissue', 'magnitude_subject', 
#                                                       'phase_subject', 'phase_tissue'))
#ggplot(df_forviolin, aes(x=vble, y=cv, fill = factor(vble)) ) + 
#  geom_violin( scale="width" ) +
#  geom_boxplot(width=0.07, fill="grey", outlier.colour='black') + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ggtitle('Is this plot confusing? It somehow confuses me.. Can we compare different CVs?\nwhat about CV_Phase?')

# Fractions of variance 
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

fig2C_2 <- ggplot(df_fraction_variance) + geom_boxplot(aes(y=value, x=variable, fill=variable), alpha=0.6) + 
  facet_wrap(~rhythm_par) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) 
  

clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
ZZ_genes <- c(ZZ_genes_D, ZZ_genes_E) %>% unique

df_fraction_variance$rhythm_par = factor(df_fraction_variance$rhythm_par, levels=c('magnitude','amplitude','phase'))
fig2C_1 <- ggplot(df_fraction_variance %>% filter(variable=="subject")) +
  geom_histogram(aes(value, y=..density.., fill=variable), 
                 alpha=0.6, position='identity', colour="white", bins=50) +
  geom_density(aes(value, fill=variable, group=variable), alpha=0.3) + #, position="fill"
  #geom_text_repel(
  #  data=filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject")%>%
  #    mutate(Symbol_it = paste0("italic('", Symbol, "')")),
  #  aes(
  #    x = value,
  #    y = 2.5,
  #    label=Symbol_it
  #  ),
  #  size=2.8, parse = TRUE, max.overlaps = Inf, box.padding = 1.
  #) + 
  #stat_bin(bins=50) + ylim(c(0, 3)) +  
  #stat_bin(bins=50, geom="text", aes(label=..count..), vjust=-1.5) 
  #geom_vline(data=filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject"), 
  #           aes(xintercept=value), color="black", size=.1, linetype="dashed") +
  facet_wrap(~rhythm_par, scales="free") + 
  coord_cartesian(ylim=c(0, 3.5))+ # ylim(0, 3.5) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none") +
  labs(x='Fraction of variance explained by subject', y='Density')  +
  theme(panel.spacing = unit(2.8, "lines"))
#https://ggrepel.slowkow.com/articles/examples.html#align-labels-on-the-left-or-right-edge-1

# bin information of magnitude panel
info_hist_magn <- layer_data(fig2C, 1) %>% filter(PANEL==1)
only_CGs_magn <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="magnitude") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_magn <- c()
for (v in only_CGs_magn$value){
  m <- (info_hist_magn$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_magn <- append(idxs_magn, m2)
  print(m2)
}
only_CGs_magn$hist_xvalue <- info_hist_magn[idxs_magn, "x"]
only_CGs_magn$hist_yvalue <- info_hist_magn[idxs_magn, "y"]
only_CGs_magn$rhythm_par <- "magnitude"

# bin information of amplitude panel
info_hist_amp <- layer_data(fig2C, 1) %>% filter(PANEL==2)
only_CGs_amp <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="amplitude") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_amp <- c()
for (v in only_CGs_amp$value){
  m <- (info_hist_amp$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_amp <- append(idxs_amp, m2)
  print(m2)
}
only_CGs_amp$hist_xvalue <- info_hist_amp[idxs_amp, "x"]
only_CGs_amp$hist_yvalue <- info_hist_amp[idxs_amp, "y"]
only_CGs_amp$rhythm_par <- "amplitude"


# bin information of phase panel
info_hist_phi <- layer_data(fig2C, 1) %>% filter(PANEL==3)
only_CGs_phi <- filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject" & rhythm_par=="phase") %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% dplyr::select(Symbol, Symbol_it, rhythm_par, value) 
idxs_phi <- c()
for (v in only_CGs_phi$value){
  m <- (info_hist_phi$xmax - v) > 0
  m2 <- min(which(m == TRUE))
  idxs_phi <- append(idxs_phi, m2)
  print(m2)
}
only_CGs_phi$hist_xvalue <- info_hist_phi[idxs_phi, "x"]
only_CGs_phi$hist_yvalue <- info_hist_phi[idxs_phi, "y"]
only_CGs_phi$rhythm_par <- "phase"

only_CGs <- rbind(only_CGs_magn, only_CGs_amp, only_CGs_phi)

df_fraction_variance$rhythm_par = factor(df_fraction_variance$rhythm_par, levels=c('magnitude','amplitude','phase'))
only_CGs$rhythm_par = factor(only_CGs$rhythm_par, levels=c('magnitude','amplitude','phase'))
only_CGs$hist_yvalue <- ifelse(only_CGs$hist_yvalue>3.5, 3.5, only_CGs$hist_yvalue)
fig2C <- fig2C_1 +   
  geom_text_repel(
    data=only_CGs,# %>% filter(Symbol %in% c("ARNTL", "CRY1",  "DBP"  , "NPAS2")),
    aes(x = hist_xvalue, y = hist_yvalue, label=Symbol_it),
    size=2.8, parse = TRUE, max.overlaps = Inf, box.padding = 1.,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 1,
    direction    = "x",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1, segment.color="grey50")  #+
  geom_text_repel(
    data=only_CGs %>% filter(Symbol %in% c("NR1D1", "NR1D2", "PER1",  "PER2",  "PER3" )),
    aes(x = hist_xvalue, y = hist_yvalue, label=Symbol_it),
    size=2.8, parse = TRUE, max.overlaps = Inf, box.padding = 1.,
    #force_pull   = 0.05, # do not pull toward data points
    nudge_y      = 2,
    direction    = "x",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1, segment.color="grey")  

###########

fig2C_3 <- ggplot(df_fraction_variance ) +
  #geom_histogram(aes(value, y=..density..), 
  #               alpha=0.6, position='identity', colour="black", fill="white", bins=50) +
  stat_ecdf(aes(x = value, len=dim(df_total)[1], y = ..y..*len, color=variable), geom = "step") + 
  geom_vline(data=filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject"), 
             aes(xintercept=value), color="black", size=.1, linetype="dashed") +
  geom_text_repel(data=filter(df_fraction_variance, Symbol %in% clock_genes & variable=="subject") %>%
                    mutate(Symbol_it = paste0("italic('", Symbol, "')")),
                  aes(x=value, y=1000, label=Symbol_it), angle=90, size=2.5,
                  max.overlaps=Inf, box.padding=1., point.padding=.5, parse=TRUE, color="grey20") +
  facet_wrap(~rhythm_par) +
  scale_color_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
  labs(x='Fraction of variance explained by subject', y='Number of genes')  

#
## from the top 20 genes with highest amp_vbility across subjects, how does the variance in phase/magn across subj look like?
#ggpubr::ggarrange(
#  plotPercentBars(amp_frac %>% arrange(desc(amplitude_tissue)) %>% head(20) %>% arrange(rownames(.)) %>%
#                    dplyr::rename(c("subject"="amplitude_subject", "tissue"="amplitude_tissue"))
#                  ) + 
#    ylab('% of total variance') + 
#    scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
#    theme(axis.text.y = element_text(face="italic")) +
#    ggtitle('Amplitude'),
#  plotPercentBars(phi_frac %>% filter(rownames(.) %in% rownames(amp_frac %>% arrange(desc(amplitude_tissue)) %>% 
#                                                                  head(20))) %>% arrange(rownames(.)) %>%
#                    dplyr::rename(c("subject"="phase_subject", "tissue"="phase_tissue"))
#                  ) + 
#    ylab('% of total variance') + 
#    scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
#    theme(axis.text.y = element_text(face="italic")) +
#    ggtitle('Phase'),
#  plotPercentBars(magn_frac %>% filter(rownames(.) %in% rownames(amp_frac %>% arrange(desc(amplitude_tissue)) %>% 
#                                                                   head(20))) %>% arrange(rownames(.)) %>%
#                    dplyr::rename(c("subject"="magnitude_subject", "tissue"="magnitude_tissue"))
#                  ) + 
#    ylab('% of total variance') + 
#    scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
#    theme(axis.text.y = element_text(face="italic")) +
#    ggtitle('Magnitude'),
#  nrow=1, ncol=3, legend="right", common.legend=TRUE)

#ggpubr::ggarrange(plotPercentBars(phi_frac %>% arrange(desc(phase_subject)) %>% head(20)) +
#                    ylab('% of total variance'),
#                  plotPercentBars(phi_frac %>% arrange(desc(phase_tissue)) %>% head(20)) +
#                    ylab('% of total variance'),
#                  nrow=1, ncol=2, common.legend=TRUE, legend="right")


# Correlations of cv_amp, cv_phi and cv_magn
df <- variation_A %>% mutate(variable=ifelse(variable=="A_T", "tissue", "subject")) %>% 
  dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_amp"="cv")) %>% 
  full_join(variation_phi %>% mutate(variable=ifelse(variable=="phi_T", "tissue", "subject")) %>% 
              dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_phi"="cv"))) %>%
  full_join(variation_magn %>% mutate(variable=ifelse(variable=="magn_T", "tissue", "subject")) %>% 
              dplyr::select(Symbol, variable, cv) %>% dplyr::rename(c("cv_magn"="cv"))) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))

  
df$variable <- as.factor(df$variable, levels=c("tissue", "subject"))
pairs(df[,c(3:5)],
      col = alpha(c("#00798c", "#d1495b"), 0.4)[df$variable],   
      pch = 19,                                                
      cex = 0.6,       
      labels = c("cv_amp","cv_phi","cv_magn"),  
      gap = 0.3, 
      upper.panel = NULL)

# Same plot of cv_phi vs cv_amp -> where are the DR/non_DR genes located?
DR_genes    <- filter(some_rhy, diff_rhythmic)$Symbol
nonDR_genes <- filter(some_rhy, diff_rhythmic==FALSE)$Symbol

df$variable <- factor(df$variable, levels=c("tissue", "subject"))
fig2B <- ggplot(df %>% mutate(clock_gene = ifelse(Symbol %in% clock_genes, TRUE, FALSE))) +
  geom_point(aes(x=cv_amp, y=cv_phi), alpha=0.3, color="grey", size=.8) + 
  geom_point(data=filter(df, Symbol %in% DR_genes),      
             aes(x=cv_amp, y=cv_phi, color=variable), alpha=0.5,size=.8) + 
  facet_wrap(~variable, scales="free") +
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
                  parse=TRUE, size=2.8, segment.color="grey") +
                  #max.overlaps=Inf, box.padding=1.1, point.padding=.5, parse=TRUE, color="black") +
  scale_color_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none") +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
  #ggtitle(paste(length(which(df_total$ProbeName %in% filter(some_rhy, diff_rhythmic)$ProbeName)),
  #              '/', dim(df_total)[1], ' amp-filtered vP genes are DR')) +
  labs(x='amplitude coefficient of variation', y='phase coefficient of variation') +
  theme(legend.position = "right",
        strip.text = element_blank(),
        panel.spacing = unit(2.5, "lines")) +
  scale_y_continuous(limits=c(0,0.08)) + scale_x_continuous(limits=c(0,1.8))


fig2 <- plot_grid(fig2A_2, NULL, 
                  plot_grid(NULL, fig2B, NULL, ncol=3, rel_widths = c(.015,1,0.25)), 
                  NULL, fig2C, nrow=5, 
                  labels = c("A", "", "B", "", "C"), rel_heights = c(1,0.1,.9,0.1,1),align = "v")
#fig2 <- plot_grid(plot_grid(fig2A_2, fig2B, ncol=2, rel_widths = c(.6,.4), labels = c("A", "B")), 
#                  NULL, fig2C, nrow=3, 
#                  labels = c("", "", "C"), rel_heights = c(1,0.1,1),align = "v")
fig2 %>% ggsave('figures/fig2.pdf', ., width = 11, height = 9)


# Correlations of variances vs means
variation_full$rhythmic_par = factor(variation_full$rhythmic_par, levels=c('magnitude','amplitude','phase'))
variation_full$variable= factor(variation_full$variable, levels=c('tissue','subject'))
variation_full %<>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))

suppfig3_part1 <- ggplot(variation_full %>% dplyr::select(rhythmic_par, value_fit, variable)
#                   gather(effect, value, -ProbeName, -rhythmic_par, -variable, -cv, -sd, -Symbol, -Symbol_it)
                   ) +
  geom_histogram(aes(value_fit, fill=variable), color='white', bins=50, alpha=0.6) + 
  facet_wrap(~rhythmic_par, scales="free") +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none") #+
  ggtitle('distribution of fixed effects')

suppfig3_part2 <- ggplot(variation_full) + 
  facet_wrap(variable ~ rhythmic_par, scales="free") +
  facetted_pos_scales(
    y = list(scale_y_continuous(limits = c(0, 25)), scale_y_continuous(limits = c(0, 0.3)), scale_y_continuous(limits = c(0, 3)),
             scale_y_continuous(limits = c(0, 1.5)), scale_y_continuous(limits = c(0, 0.3)), scale_y_continuous(limits = c(0, 3)))
  ) +
  geom_point(aes(x= value_fit, y=variance), alpha=0.3, color="grey", size=0.6) + 
  geom_point(data=filter(variation_full, Symbol %in% DR_genes), aes(x= value_fit, y=variance, color=variable), 
             alpha=0.5, size=0.6)  +
  geom_point(data=filter(variation_full, Symbol %in% clock_genes & Symbol %in% nonDR_genes),      
             aes(x= value_fit, y=variance), alpha=1, shape=21, color="black", size=1, fill="grey20") +
  geom_point(data=filter(variation_full, Symbol %in% clock_genes & Symbol %in% DR_genes),      
             aes(x= value_fit, y=variance, fill=variable), alpha=1, shape=21, color="grey60", size=1) +
  geom_text_repel(data=filter(variation_full, Symbol %in% clock_genes), 
                  aes(x= value_fit, y=variance, label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.1, point.padding=.5, parse=TRUE, color="black", size=2.5) +
  scale_color_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none")+
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c"), guide="none") 

suppfig3_new <- plot_grid(suppfig3_part1, NULL, suppfig3_part2, labels=c("A","","B"), nrow=3, rel_heights=c(1,0.1,2))
suppfig3_new %>% ggsave('figures/suppfig3_new.pdf', ., width = 11, height = 9)
###############################################
#up to here 11.04.2022

# Time-telling genes
ZZ_genes_D <- c("KLF8", "PER3", "ABCC5", "RBM14", "EIF4A2", "HNRNPH3", "ZNF692", 
                "HNRNPDL", "APLN", "BECN1", "CMTM1", "OVGP1", "ARNTL", "ABCC5", "RBM14", "ZBTB16", "ADRA1B",
                "TTC32", "ZBED8", "DYNLL1", "NMNAT1", "AVIL", "TMEM104", "TMCO6", "IRF7") %>% unique
ZZ_genes_E <- c("ELL3", "CYP4B1", "ZNF143", "CYHR1", "RBM14", "ZNF692", "OVGP1", "POLA1", "SOX5", "FKBP5",
                "FOCAD", "ROR1", "TRIM35", "PLEKHA6", "ZBTB16") %>% unique

df_cv_magn_amp <- variation_full %>% 
  dplyr::select(rhythmic_par, variable, cv, Symbol, Symbol_it) %>% 
  filter(variable=="subject", rhythmic_par != "phase") %>% 
  spread(rhythmic_par, cv) 
df_cv_magn_phi <- variation_full %>% 
  dplyr::select(rhythmic_par, variable, cv, Symbol, Symbol_it) %>% 
  filter(variable=="subject", rhythmic_par != "amplitude") %>% 
  spread(rhythmic_par, cv) 

ggpubr::ggarrange(
  ggplot(df_cv_magn_amp) +
    geom_point(aes(x=magnitude, y=amplitude), color="#00798c", alpha=0.6) +
    geom_point(data=filter(df_cv_magn_amp, Symbol %in% ZZ_genes_D | Symbol %in% ZZ_genes_E),
               aes(x=magnitude, y=amplitude), color="red", alpha=1) +
   # geom_text_repel(data=filter(df_cv_magn_amp, Symbol %in% unique(c(ZZ_genes_D, ZZ_genes_E))),
   #                 aes(x=magnitude, y=amplitude, label=Symbol_it), 
   #                 color="black", max.overlaps=Inf, box.padding=1.1, point.padding=.5, parse=TRUE) +
    labs(x='magnitude CV', y='amplitude CV') + #BECN1 and FOCAD missing!!
    ggtitle('BECN1 and FOCAD missing!! (actually not rhy)'),
  ggplot(df_cv_magn_phi) +
    geom_point(aes(x=magnitude, y=phase), color="#00798c", alpha=0.6) +
    geom_point(data=filter(df_cv_magn_phi, Symbol %in% ZZ_genes_D | Symbol %in% ZZ_genes_E),
               aes(x=magnitude, y=phase), color="red", alpha=1) +
    #geom_text_repel(data=filter(df_cv_magn_phi, Symbol %in% unique(c(ZZ_genes_D, ZZ_genes_E))),
    #                aes(x=magnitude, y=phase, label=Symbol_it), 
    #                color="black", max.overlaps=Inf, box.padding=1.1, point.padding=.5, parse=TRUE) +
    labs(x='magnitude CV', y='phase CV') +
    ggtitle('BECN1 and FOCAD missing!! (actually not rhy)'),
  nrow=1, ncol=2)

ZZD_fraction_variance <- df_fraction_variance %>% filter(Symbol %in% ZZ_genes_D) %>% arrange(desc(value))
ZZE_fraction_variance <- df_fraction_variance %>% filter(Symbol %in% ZZ_genes_E) %>% arrange(desc(value))
ggplot(ZZD_fraction_variance) +
  geom_col(aes(x=value, y=reorder(Symbol, value), fill=variable)) + 
  facet_wrap(~rhythm_par, scales='free_y') +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
  theme(axis.text.y = element_text(face="italic"), aspect.ratio=1.7) +
  labs(x='Fraction of variance explained', y='')
ggplot(ZZE_fraction_variance) +
  geom_col(aes(x=value, y=Symbol, fill=variable)) + 
  facet_wrap(~rhythm_par) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c")) +
  theme(axis.text.y = element_text(face="italic"), aspect.ratio=1.7) +
  labs(x='Fraction of variance explained', y='')



##############################
##############################

## Figure 2B: How does each variable contribute to the variance 
## ------------------------------------------------------------
#fig2B <- plotVarPart(vp, label.angle=60) + 
#  theme_custom() + ylab("Percentage of\nvariance explained") + 
#  theme(aspect.ratio=0.5, legend.position = "none", ) + 
#  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", 
#                               "time" = "#edae49", "time:tissue" = "#F78828", 
#                               "time:subject" = "#39B600", "Residuals" = "grey80")) +
#  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subject\nmean\nvariation", 
#                               "time" = "Common\ncircadian\nvariation", "time:subject" = "Inter-subject\ncircadian\nvariation",
#                               "time:tissue" = "Inter-layer\ncircadian\nvariation", "Residuals" = "Residual\nvariation"))
#### The time variation doesn't seem to be very tissue- or subject-specific (since the variance of Time:Tissue
#### and Time:Subject is < than the Time variation). In other words, time can explan variation without tissue or subject


# Suppl. Figure 3A: Negative control -- time variance of 1000 arrhythmic genes expected ~0
# -----------------------------------------------------------------------------------
no_rhy <- dplyr::filter(results, adj_P_Val_D_or_E > 0.1)
geneExpr_arrhy <- yave$E %>% as.data.frame %>% filter(rownames(.) %in% no_rhy$ProbeName)

fitList <- fitVarPartModel(geneExpr_arrhy[1:1000,], form, info_exp, showWarnings=FALSE) 

# Calculate variance in amplitude and phase from all rhythmic genes -> Error propagation
df_total = data.frame(ProbeName=NULL, Amp=NULL, var_A_layer=NULL, var_A_subject=NULL,
                      phase=NULL, var_phi_layer=NULL, var_phi_subject=NULL, 
                      magn=NULL, var_magn_layer=NULL, var_magn_subject=NULL)

for (i in 1:length(fitList)){
  
  # fixed effects (coefficients)
  m_i <- fixef(fitList[[i]])[1] %>% as.numeric 
  a_i <- fixef(fitList[[i]])[2] %>% as.numeric 
  b_i <- fixef(fitList[[i]])[3] %>% as.numeric
  
  # amplitude and phase of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix sigma (because there are >1 subject/tissue to do fits)
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList[[i]])))
  sigma_i_S <- sigma_i[1:3,1:3] #covariance_subject
  sigma_i_T <- sigma_i[4:6,4:6] #covariance_tissue
  
  # define jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), # dA/dm, dA/da, dA/db
                   nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), 
                     nrow=1)
  Jac = matrix( c(0, a_i/A_i, b_i/A_i,
                  0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), # dA/dm, dA/da, dA/db
                nrow = 2, byrow = TRUE) #!!!!!!!
  
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
hist(df_total$Amp, breaks=100)
#df_total %<>% filter(Amp>.20) # filter genes with low amp_fit that result in high variability and "mask" variable genes
hist(df_total$Amp, breaks=100)


# Calculate sd and cv from variances, organize results for plots
variation_A   <- df_total %>% dplyr::select(ProbeName, Amp, var_A_subject, var_A_layer) %>% 
  gather(variable, variance, -ProbeName, -Amp) %>%
  mutate(variable = ifelse(variable=="var_A_subject", "A_S", "A_T"),
         sd = sqrt(variance),
         cv = sd/Amp) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_phi <- df_total %>% dplyr::select(ProbeName, phase, var_phi_subject, var_phi_layer) %>%
  gather(variable, variance, -ProbeName, -phase) %>%
  mutate(variable = ifelse(variable=="var_phi_subject", "phi_S", "phi_T"),
         sd = sqrt(variance),
         cv = sd/24) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_magn <- df_total %>% dplyr::select(ProbeName, magn, var_magn_subject, var_magn_layer) %>%
  gather(variable, variance, -ProbeName, -magn) %>%
  mutate(variable = ifelse(variable=="var_magn_subject", "magn_S", "magn_T"),
         sd = sqrt(variance),
         cv = sd/magn) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))

variation_full <- rbind(variation_A %>% dplyr::rename(c("value_fit"="Amp")), 
                        variation_phi %>% dplyr::rename(c("value_fit"="phase")) %>%
                          mutate(value_fit = value_fit + 8,
                                 value_fit = ifelse(value_fit < 0, value_fit + 24, value_fit))) %>% 
  rbind(variation_magn %>% dplyr::rename(c("value_fit"="magn"))) %>%
  tidyr::separate(variable, c("rhythmic_par","variable"), sep = "_", convert = TRUE) %>%
  mutate(rhythmic_par=ifelse(rhythmic_par == "A", "amplitude", 
                             ifelse(rhythmic_par=="phi", "phase", "magnitude")),
         variable=ifelse(variable == "S", "subject", "tissue")) 

plotVarPart(sortCols(df_total %>% dplyr::select(Symbol, var_A_layer, var_A_subject, var_phi_layer, 
                                                var_phi_subject,  var_magn_layer, var_magn_subject) %>% 
                       tibble::column_to_rownames("Symbol"))) + 
  ggtitle('WHAT ABOUT RESIDUALS, x axis is wrong\nI dont understand this plot')

df_forviolin <- variation_full %>% dplyr::mutate(vble = paste0(rhythmic_par, "_", variable)) %>% 
  dplyr::select(Symbol, vble, variance)
df_forviolin$vble = factor(df_forviolin$vble, levels=c('amplitude_subject','amplitude_tissue',
                                                       'magnitude_tissue', 'magnitude_subject', 
                                                       'phase_subject', 'phase_tissue'))
ggplot(df_forviolin, aes(x=vble, y=variance, fill = factor(vble)) ) + 
  geom_violin( scale="width" ) +
  geom_boxplot(width=0.07, fill="grey", outlier.colour='black') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Is this plot confusing? It somehow confuses me.. Can we compare different CVs?\nwhat about CV_Phase?')



varPart <- fitExtractVarPartModel(geneExpr_arrhy[1:1000,], form, info_exp) 
varPart <- varPart %>% as.data.frame() %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% tibble::column_to_rownames("Symbol") #%>%
vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID)) 

suppfig3A <- plotVarPart(vp) + #ggtitle('Variance partition on 1000 non-rhythmic genes') + 
  theme_custom() + theme(aspect.ratio=0.5, legend.position = "none") +
  ylab("Percentage of variance explained") #+ 
  #scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", 
  #                             "time" = "#edae49", "time:tissue" = "#F78828", 
  #                             "time:subject" = "#39B600", "Residuals" = "grey80")) +
  #scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subj.\nmean\nvariation", 
  #                            "time" = "Common\ncircadian\nvariation", "time:subject" = "Inter-subj.\ncircadian\nvariation",
  #                            "time:tissue" = "Inter-layer\ncircadian\nvariation", "Residuals" = "Residual\nvariation"))

# Suppl. Figure 3B: positive control -- high time variance of clock genes expected 
# --------------------------------------------------------------------------------
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
geneExpr_clockgenes <- geneExpr %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% 
  filter(Symbol %in% clock_genes) %>% tibble::column_to_rownames("ProbeName") %>%
  dplyr::select(-Symbol, -EntrezID)

varPart <- fitExtractVarPartModel(geneExpr_clockgenes, form, info_exp) 
varPart <- varPart %>% as.data.frame() %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% tibble::column_to_rownames("Symbol")
vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID)) 

suppfig3B <- plotVarPart(vp) + #ggtitle('Variance partition on clock genes') + 
  theme_custom() + theme(aspect.ratio=0.5, legend.position = "none") +
  ylab("Percentage of variance explained") + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", 
                               "time" = "#edae49", "time:tissue" = "#F78828", 
                               "time:subject" = "#39B600", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subj.\nmean\nvariation", 
                              "time" = "Common\ncircadian\nvariation", "time:subject" = "Inter-subj.\ncircadian\nvariation",
                              "time:tissue" = "Inter-layer\ncircadian\nvariation", "Residuals" = "Residual\nvariation"))

#---------------------------------
#---------------------------------

# 3. VARIANCE PARTITION IN DERMIS VS. EPIDERMIS
# ---------------------------------------------
form <- ~ (1|time) + (1|subject) #+ (1|(time:subject)) -> not possible because of only 1 value per 'group'

geneExpr.D <- yave$E %>% as.data.frame() %>% dplyr::select(contains("D")) %>%   #rhy genes in dermis
  filter(rownames(.) %in% filter(some_rhy, A_D > amp_cutoff & adj_P_Val < fdr_cutoff)$ProbeName)
geneExpr.E <- yave$E %>% as.data.frame() %>% dplyr::select(contains("E")) %>%   #rhy genes in epidermis
  filter(rownames(.) %in% filter(some_rhy, A_E > amp_cutoff & adj_P_Val < fdr_cutoff)$ProbeName)

varPart.D <- fitExtractVarPartModel(geneExpr.D, form, info_exp %>% filter(tissue=="dermis"), showWarnings=TRUE) %>%
  as.data.frame() %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% tibble::column_to_rownames("Symbol") 
varPart.E <- fitExtractVarPartModel(geneExpr.E, form, info_exp %>% filter(tissue=="epidermis")) %>%
  as.data.frame() %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% tibble::column_to_rownames("Symbol") 

vp.D <- sortCols(varPart.D %>% dplyr::select(-ProbeName, -EntrezID)) 
vp.E <- sortCols(varPart.E %>% dplyr::select(-ProbeName, -EntrezID)) 


# Figures 2C and 2D: variance partition done in dermis vs epidermis separately
# ----------------------------------------------------------------------------
fig2D <- plotVarPart(vp.D) + 
  theme_custom() + theme(aspect.ratio=0.6, legend.position = "none") +
  ylab("Percentage of\nvariance explained\n(rhythmic genes in dermis)") +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#1B9E77", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subject\nmean dermis\nvariation", 
                              "time" = "Dermis\ncircadian\nvariation", "Residuals" = "Residual\ndermis\nvariation"))

fig2E <- plotVarPart(vp.E) + 
  theme_custom() + theme(aspect.ratio=0.6, legend.position = "none") +
  ylab("Percentage of\nvariance explained\n(rhythmic genes in epidermis)") +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#D95F02", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subject\nmean epid.\nvariation", 
                              "time" = "Epidermis\ncircadian\nvariation", "Residuals" = "Residual\nepidermis\nvariation"))


if (!file.exists("visualize/data/variancePartition_dermis.csv")){
  write.csv(vp.D %>% tibble::rownames_to_column("Symbol"), "visualize/data/variancePartition_dermis.csv")
  write.csv(vp.E %>% tibble::rownames_to_column("Symbol"), "visualize/data/variancePartition_epidermis.csv")}

  
# Figures 2E, F: Genes with highest time-variance in each tissue
# --------------------------------------------------------------
number_most_vargenes <- 20
most_timevar.D   <- varPart.D %>% dplyr::arrange(desc(time)) %>% head(number_most_vargenes)
most_timevar.E   <- varPart.E %>% dplyr::arrange(desc(time)) %>% head(number_most_vargenes) 

# order data frames so that common genes appear at end of plot
common_genes <- inner_join(most_timevar.D %>% tibble::rownames_to_column() %>% dplyr::select(rowname, ProbeName, EntrezID), 
                           most_timevar.E %>% tibble::rownames_to_column() %>% dplyr::select(rowname, ProbeName, EntrezID)) 
most_timevar.D <- rbind(anti_join(most_timevar.D  %>% tibble::rownames_to_column(), common_genes), 
                        inner_join(most_timevar.D %>% tibble::rownames_to_column(), common_genes) %>% arrange(rowname)) %>%
  tibble::column_to_rownames("rowname")
most_timevar.E <- rbind(anti_join(most_timevar.E  %>% tibble::rownames_to_column(), common_genes),  
                        inner_join(most_timevar.E %>% tibble::rownames_to_column(), common_genes) %>% arrange(rowname))  %>%
  tibble::column_to_rownames("rowname")

fig2F <- plotPercentBars(most_timevar.D %>% dplyr::select(-ProbeName, -EntrezID)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="bottom", aspect.ratio=2.3,
        rect = element_rect(fill="transparent")) + 
  ylab('Percentage of variance explained\n(rhythmic genes in dermis)') + 
  geom_vline(xintercept=dim(common_genes)[1]+0.5) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#1B9E77", "Residuals" = "grey80"),
                    labels = c("tissue" = "Inter-layer\nmean dermis\nvariation", 
                               "subject" = "Inter-subject\nmean dermis\nvariation", 
                               "time" = "Dermis\ncircadian\nvariation", "Residuals" = "Residual\ndermis\nvariation")) 

fig2G <- plotPercentBars(most_timevar.E %>% dplyr::select(-ProbeName, -EntrezID)) + 
  theme_custom() +
  theme(axis.text.y = element_text(face="italic"), legend.position="bottom", aspect.ratio=2.3, 
        rect = element_rect(fill="transparent")) + 
  ylab('Percentage of variance explained\n(rhythmic genes in epidermis)') + 
  geom_vline(xintercept=dim(common_genes)[1]+0.5) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#D95F02", "Residuals" = "grey80"),
                    labels = c("tissue" = "Inter-layer\nmean epid.\nvariation", 
                               "subject" = "Inter-subject\nmean epid.\nvariation", 
                              "time" = "Epidermis\ncircadian\nvariation", "Residuals" = "Residual\nepidermis\nvariation"))

fig2FG <- plot_grid(NULL, fig2F, NULL, fig2G, nrow=1, ncol=4, labels=c("F", "", "G", ""), rel_widths=c(-0.2,1,-0.2,1))


##############################
##############################


# Arrange plots in grid
# ---------------------
fig2_1 <- plot_grid(NULL, fig2A, NULL, fig2B, nrow=1, ncol=4, labels=c("A", "", "", "B"), rel_widths=c(0.095,0.66,0.1,1))
fig2_2 <- plot_grid(fig2C, labels=c('C'))
fig2_3 <- plot_grid(fig2D, fig2E, nrow=2, ncol=1, labels=c("D", "E"), rel_heights=c(1,1), align='v')
fig2_4 <- plot_grid(fig2_3, fig2FG, ncol=2, nrow=1, rel_widths=c(0.66,1))
fig2   <- plot_grid(fig2_1, fig2_2, fig2_4, nrow=3, ncol=1, rel_heights=c(0.7,1.,1.), align='v')
fig2 %>% ggsave('figures/fig2.pdf', ., width = 11, height = 11.5)


sfig3_1 <- plot_grid(NULL, suppfig3A, NULL, suppfig3B, nrow=1, ncol=4, labels=c("A", "", "B", ""), rel_widths=c(0.1,1,0.1,1))
sfig3_2 <- plot_grid(suppfig3C, nrow=1, ncol=1, labels=c("C"))
sfig3_part1 <- plot_grid(sfig3_1, NULL, sfig3_2, nrow=3, rel_heights=c(0.6,0.1,0.85))
sfig3_part1 %>% ggsave('figures/suppfig3_1.pdf', ., width = 11, height = 9)

sfig3_3 <- plot_grid(NULL, suppfig3D, nrow=1, ncol=2, labels=c("D", ""), rel_widths=c(0.08,1))
sfig3_4 <- plot_grid(NULL, suppfig3E, nrow=1, ncol=2, labels=c("E", ""), rel_widths=c(0.02,1))
sfig3_part2 <- plot_grid(sfig3_3, NULL, sfig3_4, nrow=3, rel_heights=c(0.85,0.1,1.1))
sfig3_part2 %>% ggsave('figures/suppfig3_2.pdf', ., width = 11, height = 12)
