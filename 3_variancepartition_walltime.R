suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(hgug4112a.db))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr) )
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(TissueEnrich))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidytext))
suppressPackageStartupMessages(library(lme4))
#suppressPackageStartupMessages(library(clusterProfiler))

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

results <- readRDS("visualize/data/results_populationrhy_walltime.rds") %>% 
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
info_exp <- info_exp %>% full_join(experiment) %>% mutate(tissue = ifelse(tissue=="D", "dermis", "epidermis"),
                                                          time = as.character(time))

# Figure 2A (introductory panel)
# ------------------------------
geneExpr_gath <- geneExpr %>% tibble::rownames_to_column("gene") %>% gather(key, value, -gene) %>% 
  tidyr::separate(key, c("tissuetime","subject"), sep = "_", convert=TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1)  %>% 
  full_join(experiment %>% dplyr::select(tissue, time, subject)) %>% mutate(time=as.character(time)) 

df1 <- geneExpr_gath %>% gather(key, category, -gene, -value) 
df2 <- df1 %>% group_by(category) %>% 
  dplyr::summarise(magnitude = mean(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE)) %>% as.data.frame() %>%
  inner_join(df1 %>% dplyr::select(category, key)) %>% unique()
df3 <- df2 %>% group_by(key)%>% 
  dplyr::summarise(mean_magn = mean(magnitude, na.rm=TRUE), sd_magn = sd(magnitude, na.rm=TRUE)) %>% as.data.frame()
df3 %<>% mutate(key=paste0(key,"_avg")) %>% rename(c("magnitude"="mean_magn", "sd"="sd_magn")) 

df <- gtools::smartbind(df2, df3) #%>%
temp = data.frame(key=unique(df$key))
temp$group = sub("_.*", "", temp$key)

df %<>% full_join(temp) %>% arrange(group) 
#df %<>%
#  group_split(group) %>% 
#  purrr::map_dfr(~ add_row(.x, .after = Inf)) %>% as.data.frame()
#df$key = replace_na(df$key, "empty")

fig2A <- ggplot() + #geom_vline(xintercept=dplyr::summarise(df2, mean=mean(magnitude))$mean, linetype="dashed", size=0.5, color="gray80") +
  geom_point(data = df %>% filter(key %in% c("time", "tissue", "subject")), 
             aes(x=magnitude, y=key, group=key, color=group),
             shape=21, size=2, position = position_dodge(width = 1.9), show.legend = F) +
  geom_errorbar(data= df %>% filter(key %in% c("time_avg", "tissue_avg", "subject_avg")),
                aes(y=key, xmin = magnitude-sd, xmax = magnitude + sd, color=group), width = 0.3, size=1,
                position = position_dodge(width = 1.9), show.legend = F) +
  geom_point(data= df %>% filter(key %in% c("time_avg", "tissue_avg", "subject_avg")),
             aes(x=magnitude, y=key, color=group), size=4.5, position = position_dodge(width = 1.9)) +
  scale_color_manual(labels = c("tissue" = "Inter-layer\nmean variation", "subject" = "Inter-subject\nmean variation", 
                                "time" = "Common\ncircadian variation"),
                     values = c("tissue" = "#d1495b", "subject" = "#00798c", "time" = "#edae49" )) + #facet_grid(rows=vars(group), scales="free") + 
  theme(legend.position="right", 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        strip.text = element_blank(),
        aspect.ratio=1,
        panel.spacing.y = unit(-15, "lines"),
        legend.spacing.y = unit(1.0, 'cm')) + 
  scale_y_discrete(expand=expansion(c(.6,.6))) +xlab('\ngene expression levels')

#--------------------------------
#--------------------------------


# 2. EXECUTE VARIANCE PARTITION MODEL (WITH INTERACTION TERMS)
# ------------------------------------------------------------
info_exp %<>% mutate(MSF_sc_dec = hms(MSF_sc),
                     MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + second(MSF_sc_dec) / 360),2),
                     inphase = cos(2*pi*as.numeric(time)/24),
                     outphase = sin(2*pi*as.numeric(time)/24))
#form <- ~ (1|subject) + (1|tissue) + (1|time) + (1|time:tissue) + (1|time:subject) 
#        #time here is wall time (continuous vble [internal_time] can't be modeled as random)
form <- ~ (1 + inphase + outphase|subject) + (1 + inphase + outphase|tissue)
form <- ~ 1 + inphase + outphase + (1 + inphase + outphase|subject) + (1 + inphase + outphase|tissue)

varPart <- fitExtractVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE) 
varPart <- varPart %>% as.data.frame() %>%  
  mutate(ProbeName=rownames(.)) %>% inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>%
  tibble::column_to_rownames("Symbol")
vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID))
    #sortCols: orders columns such that the first one is the one with largest (mean) variation, second is second largest, etc

if (!file.exists("visualize/data/variancePartition_full.csv")){
  write.csv(vp %>% tibble::rownames_to_column("Symbol"), "visualize/data/variancePartition_full.csv")
}
# Save 50 genes with least variability across subjects -> supplementary Table 3
if (!file.exists("figures/supp_table3.csv")){
  write.csv(varPart %>% dplyr::arrange(subject) %>% head(50) %>% tibble::rownames_to_column("Symbol"), 
            "figures/supp_table3.csv")
}
##############################
###NEW
fitList <- fitVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE) 

#for gene i=1
m_i <- fixef(fitList[[1]])[1]; a_i <- fixef(fitList[[1]])[2]; b_i <- fixef(fitList[[1]])[3]
A_i <- sqrt(a_i^2 + b_i^2) %>% as.numeric

Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), # dA/dm, dA/da, dA/db
                 nrow = 1)
sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList[[1]])))
sigma_i_S <- sigma_i[1:3,1:3]
sigma_i_T <- sigma_i[4:6,4:6]

var_A_S <- Jac_A %*% sigma_i_S %*% t(Jac_A); sd_A_S <- sqrt(var_A_S); cv_A_S <- sd_A_S/A_i 
var_A_T <- Jac_A %*% sigma_i_T %*% t(Jac_A); sd_A_T <- sqrt(var_A_T); cv_A_T <- sd_A_T/A_i 

# calculate variance in amplitude from all rhythmic genes
df_total = data.frame(ProbeName=NULL, Amp=NULL, var_A_T=NULL, var_A_S=NULL,
                      phase=NULL, var_phi_T=NULL, var_phi_S=NULL)
#df_total = data.frame(var_A_T=NULL, sd_A_T=NULL, cv_A_T=NULL, var_A_S=NULL, sd_A_S=NULL, cv_A_S=NULL)
for (i in 1:length(fitList)){
  # fixed effects (coefficients)
  m_i <- fixef(fitList[[i]])[1] %>% as.numeric; 
  a_i <- fixef(fitList[[i]])[2] %>% as.numeric; 
  b_i <- fixef(fitList[[i]])[3] %>% as.numeric
  
  # amplitude and phase of the fit for gene i across subjects and tissues
  A_i   <- sqrt(a_i^2 + b_i^2) 
  phi_i <- atan2(b_i, a_i)*12/pi
  
  # with estimation of coefficients comes a covariance matrix sigma (because there are >1 subject/tissue to do fits)
  sigma_i <-  as.matrix(Matrix::bdiag(VarCorr(fitList[[1]])))
  sigma_i_S <- sigma_i[1:3,1:3] #covariance_subject
  sigma_i_T <- sigma_i[4:6,4:6] #covariance_tissue
  
  # define jacobian
  Jac_A <- matrix( c(0, a_i/A_i, b_i/A_i), # dA/dm, dA/da, dA/db
                   nrow = 1)
  Jac_phi <- matrix( c(0, (1/(1+(b_i/a_i)^2)) * (-b_i/(a_i^2)), (1/(1+(b_i/a_i)^2)) * (1/a_i)), 
                     nrow=1)
  
  # determine variance, sd and cv in amplitude and phase, separating tissue and subject contributions
  var_A_S <- Jac_A %*% sigma_i_S %*% t(Jac_A)#; sd_A_S <- sqrt(var_A_S); cv_A_S <- sd_A_S/A_i 
  var_A_T <- Jac_A %*% sigma_i_T %*% t(Jac_A)#; sd_A_T <- sqrt(var_A_T); cv_A_T <- sd_A_T/A_i 
  var_phi_S <- Jac_phi %*% sigma_i_S %*% t(Jac_phi)
  var_phi_T <- Jac_phi %*% sigma_i_T %*% t(Jac_phi)
  
  #df_i <- data.frame(var_A_T=var_A_T, sd_A_T=sd_A_T, cv_A_T=cv_A_T, var_A_S=var_A_S, sd_A_S=sd_A_S, cv_A_S=cv_A_S)
  df_i <- data.frame(ProbeName=names(fitList)[i], 
                     Amp=A_i,     var_A_T=var_A_T,     var_A_S=var_A_S,
                     phase=phi_i, var_phi_T=var_phi_T, var_phi_S=var_phi_S)
  df_total <- rbind(df_total, df_i)
}

##############################
##############################

# Figure 2B: How does each variable contribute to the variance 
# ------------------------------------------------------------
fig2B <- plotVarPart(vp, label.angle=60) + 
  theme_custom() + ylab("Percentage of\nvariance explained") + 
  theme(aspect.ratio=0.5, legend.position = "none", ) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", 
                               "time" = "#edae49", "time:tissue" = "#F78828", 
                               "time:subject" = "#39B600", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subject\nmean\nvariation", 
                               "time" = "Common\ncircadian\nvariation", "time:subject" = "Inter-subject\ncircadian\nvariation",
                               "time:tissue" = "Inter-layer\ncircadian\nvariation", "Residuals" = "Residual\nvariation"))
### The time variation doesn't seem to be very tissue- or subject-specific (since the variance of Time:Tissue
### and Time:Subject is < than the Time variation). In other words, time can explan variation without tissue or subject


# Learning about the most varying genes in each variable - Figure 2B: Which are the genes that vary most in each category?
# ------------------------------------------------------------------------------------------------------------------------
number_most_vargenes <- 10

most_tissuevar   <- varPart %>% dplyr::arrange(desc(tissue)) %>% mutate(var="Tissue")
most_subjvar     <- varPart %>% dplyr::arrange(desc(subject)) %>% mutate(var="Subject")
most_timevar     <- varPart %>% dplyr::arrange(desc(time)) %>% mutate(var="Time")
most_resvar      <- varPart %>% dplyr::arrange(desc(Residuals)) %>% mutate(var="Residuals")
most_timesubjvar <- varPart[order(-varPart[,"time:subject"]),] %>% mutate(var="Time:Subject")
most_timetissvar <- varPart[order(-varPart[,"time:tissue"]),] %>% mutate(var="Time:Tissue")
most_var <- rbind(most_tissuevar[1:number_most_vargenes,], most_subjvar[1:number_most_vargenes,], 
                  most_timevar[1:number_most_vargenes,], most_resvar[1:number_most_vargenes,],
                  most_timesubjvar[1:number_most_vargenes,], most_timetissvar[1:number_most_vargenes,])

fig2C_1 <- plotPercentBars(most_var %>% filter(var=="Tissue" | var=="Subject") %>% 
                             dplyr::select(-ProbeName, -EntrezID, -var)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="right", aspect.ratio=1.7) + 
  ylab(' ') + 
  geom_vline(xintercept=c(10.5, 20.5, 30.5)) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", "time" = "#edae49", "Residuals" = "grey80",
                               "time:tissue" = "#F78828", "time:subject" = "#39B600"),
                    labels=c("tissue" = "Inter-layer mean variation", "subject" = "Inter-subject mean variation", 
                             "time" = "Common circadian variation", "time:subject" = "Inter-subject circadian variation",
                             "time:tissue" = "Inter-layer circadian variation", "Residuals" = "Residual variation")) 

fig2C_2 <- plotPercentBars(most_var %>% filter(var=="Time" | var=="Residuals") %>% 
                             dplyr::select(-ProbeName, -EntrezID, -var)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="right", aspect.ratio=1.7) + 
  ylab('Percentage of variance explained') + geom_vline(xintercept=c(10.5, 20.5, 30.5)) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", "time" = "#edae49", "Residuals" = "grey80",
                               "time:tissue" = "#F78828", "time:subject" = "#39B600"),
                    labels=c("tissue" = "Inter-layer mean variation", "subject" = "Inter-subject mean variation", 
                             "time" = "Common circadian variation", "time:subject" = "Inter-subject circadian variation",
                             "time:tissue" = "Inter-layer circadian variation", "Residuals" = "Residual variation"))

fig2C_3 <- plotPercentBars(most_var %>% filter(var=="Time:Subject" | var=="Time:Tissue") %>% 
                             dplyr::select(-ProbeName, -EntrezID, -var)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="right", aspect.ratio=1.7, 
        rect = element_rect(fill="transparent")) + 
  ylab(' ') + 
  geom_vline(xintercept=c(10.5, 20.5, 30.5)) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", "time" = "#edae49", "Residuals" = "grey80",
                               "time:tissue" = "#F78828", "time:subject" = "#39B600"),
                    labels=c("tissue" = "Inter-layer mean variation", "subject" = "Inter-subject mean variation", 
                             "time" = "Common circadian variation", "time:subject" = "Inter-subject circadian variation",
                             "time:tissue" = "Inter-layer circadian variation", "Residuals" = "Residual variation"))

fig2C <- ggpubr::ggarrange(fig2C_1, NULL, fig2C_2, NULL, fig2C_3, nrow=1, ncol=5, 
                           common.legend=TRUE, legend="right", widths=c(1.09,0.1,1,0.1,0.98)) 



# Learning about the most varying genes in each variable - Suppfig3C-E: GO analysis
# ---------------------------------------------------------------------------------
number_most_vargenes <- 200
most_var <- rbind(most_tissuevar[1:number_most_vargenes,], most_subjvar[1:number_most_vargenes,], 
                  most_resvar [1:number_most_vargenes,], most_timesubjvar[1:number_most_vargenes,],
                  most_timetissvar[1:number_most_vargenes,])

go <- NULL
for (i in unique(most_var$var)){
  go_i <- enrichGO(most_var %>% filter(var==i) %$% EntrezID, 
                   universe = varPart$EntrezID, ont="ALL",
                   'org.Hs.eg.db', pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>% as.data.frame() %>%
    tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
    tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/") %>%
    mutate(P.DE=as.numeric(pvalue))
  go_i <- go_i %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N) %>% as.data.frame() %>% mutate(var=i, P.DE=round(P.DE,3))
  go <- rbind(go, go_i)
}

suppfig3C <- ggplot(go %>% filter(var=="Tissue") %>% mutate(P.DE=ifelse(P.DE == 0, 0.00001, P.DE)), 
                    aes(x=-log10(P.DE), y=reorder_within(Description, -log10(P.DE), var), color=ONTOLOGY, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(ONTOLOGY), switch="y") + 
  labs(color="Ontology") + scale_y_reordered() +
  theme_bw() + ggtitle(paste0('GO analysis on top ', number_most_vargenes, ' most layer-variable genes')) +
  theme(#aspect.ratio=2.3,
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size=10, face='bold')) + 
  scale_colour_manual(values = c("BP" = "lightpink2", "CC" = "dodgerblue3", "MF" = "palevioletred4")) 

suppfig3D <- ggplot(go %>% filter(var=="Time:Tissue") %>% mutate(P.DE=ifelse(P.DE == 0, 0.00001, P.DE)), 
                    aes(x=-log10(P.DE), y=reorder_within(Description, -log10(P.DE), var), color=ONTOLOGY, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(ONTOLOGY), switch="y") + 
  labs(color="Ontology") + scale_y_reordered() +
  theme_bw() + ggtitle(paste0('GO analysis on top ', number_most_vargenes, ' most time:layer-variable genes')) +
  theme(#aspect.ratio=2.3,
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size=10, face='bold')) + 
  scale_colour_manual(values = c("BP" = "lightpink2", "CC" = "dodgerblue3", "MF" = "palevioletred4")) 

suppfig3E <- ggplot(go %>% filter(var=="Residuals") %>% mutate(P.DE=ifelse(P.DE == 0, 0.00001, P.DE)), 
                    aes(x=-log10(P.DE), y=reorder_within(Description, -log10(P.DE), var), color=ONTOLOGY, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(ONTOLOGY), switch="y") + 
  labs(color="Ontology") + scale_y_reordered() +
  theme_bw() + ggtitle(paste0('GO analysis on top ', number_most_vargenes, ' most residuals-variable genes')) +
  theme(#aspect.ratio=2.3,
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size=10, face='bold')) + 
  scale_colour_manual(values = c("BP" = "lightpink2", "CC" = "dodgerblue3", "MF" = "palevioletred4")) 



# Suppl. Figure 3A: Negative control -- time variance of arrhythmic genes expected ~0
# -----------------------------------------------------------------------------------
no_rhy <- dplyr::filter(results, adj_P_Val > 0.1)
geneExpr_arrhy <- yave$E %>% as.data.frame %>% filter(rownames(.) %in% no_rhy$ProbeName)

varPart <- fitExtractVarPartModel(geneExpr_arrhy[1:1000,], form, info_exp) 
varPart <- varPart %>% as.data.frame() %>% mutate(ProbeName=rownames(.)) %>%
  inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>% tibble::column_to_rownames("Symbol") #%>%
vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID)) 

suppfig3A <- plotVarPart(vp) + #ggtitle('Variance partition on 1000 non-rhythmic genes') + 
  theme_custom() + theme(aspect.ratio=0.5, legend.position = "none") +
  ylab("Percentage of variance explained") + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", 
                               "time" = "#edae49", "time:tissue" = "#F78828", 
                               "time:subject" = "#39B600", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subj.\nmean\nvariation", 
                              "time" = "Common\ncircadian\nvariation", "time:subject" = "Inter-subj.\ncircadian\nvariation",
                              "time:tissue" = "Inter-layer\ncircadian\nvariation", "Residuals" = "Residual\nvariation"))

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
