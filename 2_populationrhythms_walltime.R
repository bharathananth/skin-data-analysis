suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(hgug4112a.db))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr)) 
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(ggthemes))

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


#################################
#################################


# 1. SET CUTOFFS
# --------------
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R

#--------------------------------
#--------------------------------

# 2. READ FILES
# -------------
info_subjects <- read.csv("resources/info_subjects_short.csv") %>% dplyr::select(-X) # read info of subjects
experiment <- readRDS("visualize/data/experiment.rds") %>% full_join(info_subjects) # read sample details from column names
yave <- readRDS("visualize/data/rawdata.rds") # read y gene expression data (without outlier removal)
#yave_D <- readRDS("visualize/data/rawdata_dermis.rds") # read y gene expression data (without outlier removal)
#yave_E <- readRDS("visualize/data/rawdata_epidermis.rds") # read y gene expression data (without outlier removal)

# Remove outliers in yave
ind <- which(colnames(yave) == PCA_outliers)    
yave <- yave[, -ind] ; yave_E <- yave[,-(ind-77)]
experiment <- experiment[-ind,]
nrow(experiment) == ncol(yave)   #input-check
wts <- NULL #since we are removing PCA-outliers, weights from lmFit (performed later) are set to NULL 

#--------------------------------
#--------------------------------

# 3. DO LINEAR MODEL FITS AND SAVE RESULTS OF FITS
# ------------------------------------------------

# Prepare sample details for design matrix
tissue  <- factor(experiment$tissue)
time    <- experiment$time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design matrix
design <- model.matrix(~ 0 + tissue + tissue:inphase + tissue:outphase) %>% #H0: rhythms are different across tissues
  as.data.frame() %>% rename(c("tissueD_inphase"="tissueD:inphase", "tissueE_inphase"="tissueE:inphase",
                               "tissueD_outphase"="tissueD:outphase", "tissueE_outphase"="tissueE:outphase")) %>% 
  as.matrix()

# duplicate Correlations
dupcor <- duplicateCorrelation(yave, design, block=subject)

# fits
#wts <- limma::arrayWeights(yave, model.matrix(~tissue + subject))
fit <- limma::lmFit(yave, design, weights = wts, block=subject, correlation=dupcor$consensus)
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# test null hypothesis
rhy_indices <- which(grepl("phase",colnames(design)))
results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
  set_colnames(gsub("\\.","_", colnames(.))) 
results$adj_P_Val_D_or_E <- results$adj_P_Val
rhy_D_or_E <- results[(results$adj_P_Val_D_or_E < fdr_cutoff) & # fdr < cutoff
                        (pmax(sqrt(results$tissueD_inphase^2 + results$tissueD_outphase^2), #amplitude_D or amplitude_E > cutoff
                              sqrt(results$tissueE_inphase^2 + results$tissueE_outphase^2)) > amp_cutoff),]

# analysis of differentially rhythmic genes
diff_rhy_contrast <- limma::makeContrasts(tissueD_inphase - tissueE_inphase, tissueD_outphase - tissueE_outphase, levels = design)
diff_rhy_fit <- limma::contrasts.fit(fit, diff_rhy_contrast)
diff_rhy_fit <- limma::eBayes(diff_rhy_fit, robust = TRUE, trend = TRUE)
diff_rhy_results <- limma::topTable(diff_rhy_fit, number = Inf, sort.by = "none")
diff_rhy_results <- diff_rhy_results[rhy_D_or_E$ProbeName, ]

# organize results
rhy_D_or_E$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$P.Value, method = "BH")
rhy_D_or_E$diff_rhythmic <- rhy_D_or_E$adj_p_val_DR < fdr_cutoff

results <- full_join(results, rhy_D_or_E) %>% select(-adj_P_Val)

# save results of fits
if (!file.exists("visualize/data/results_populationrhy_walltime.rds")){
  saveRDS(results, file = "visualize/data/results_populationrhy_walltime.rds")
}  

#--------------------------------
#--------------------------------

# 4. READ RESULTS AND PLAY WITH THEM: 
# ----------------------------------
# We are comparing the results from the analysis done with internal time vs. the analysis done with wall time
# That's why we load the results from fig1 here as well
results_internaltime <- readRDS("visualize/data/results_populationrhy_internaltime.rds") %>%
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #note log2 values!
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi,
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)
results_internaltime$rhythmic_in_D <- ifelse(results_internaltime$diff_rhythmic==FALSE, TRUE, 
                                             ifelse(results_internaltime$diff_rhythmic==TRUE & results_internaltime$A_D > amp_cutoff, 
                                                    TRUE, FALSE)) 
results_internaltime$rhythmic_in_E <- ifelse(results_internaltime$diff_rhythmic==FALSE, TRUE, 
                                             ifelse(results_internaltime$diff_rhythmic==TRUE & results_internaltime$A_E > amp_cutoff, 
                                                    TRUE, FALSE))

results_walltime <- readRDS("visualize/data/results_populationrhy_walltime.rds") %>%
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #note log2 values!
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi,
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)
  # note phase = 0 means that the respective gene peaks at 8AM (time of first sampling)
results_walltime$rhythmic_in_D <- ifelse(results_walltime$diff_rhythmic==FALSE, TRUE, 
                                         ifelse(results_walltime$diff_rhythmic==TRUE & results_walltime$A_D > amp_cutoff, TRUE, FALSE)) 
results_walltime$rhythmic_in_E <- ifelse(results_walltime$diff_rhythmic==FALSE, TRUE, 
                                         ifelse(results_walltime$diff_rhythmic==TRUE & results_walltime$A_E > amp_cutoff, TRUE, FALSE))
results <- results_walltime


rhy_results_walltime <- results %>% filter(rhythmic_in_D | rhythmic_in_E) %>%
  dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, rhythmic_in_D, rhythmic_in_E) %>% 
  gather(tissue, rhythmic, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName,) %>% 
  mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis")) %>%
  full_join(results %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, A_D, A_E) %>% 
              gather(tissue, amp_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))) %>%
  full_join(results %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, phaseD, phaseE) %>% 
              gather(tissue, phase_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))) 
rhy_results_walltime <- rhy_results_walltime %>% anti_join( rhy_results_walltime %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) )

rhy_results_internaltime <- results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>%
  dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, rhythmic_in_D, rhythmic_in_E) %>% 
  gather(tissue, rhythmic, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName,) %>% 
  mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis")) %>%
  full_join(results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, A_D, A_E) %>% 
              gather(tissue, amp_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))) %>%
  full_join(results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, phaseD, phaseE) %>% 
              gather(tissue, phase_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))) 
rhy_results_internaltime <- rhy_results_internaltime %>% anti_join( rhy_results_internaltime %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) )


# Supplementary figure 1A: comparison of rhythmic genes with wall vs internal time
yrhy_internaltime <- list(dermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                          filter(rhy_results_internaltime, tissue=="dermis")$Symbol], 
                          epidermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                             filter(rhy_results_internaltime, tissue=="epidermis")$Symbol])
yrhy_walltime <- list(dermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                      filter(rhy_results_walltime, tissue=="dermis")$Symbol], 
                      epidermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                         filter(rhy_results_walltime, tissue=="epidermis")$Symbol])
yrhy_dermis <- list(internal = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                       filter(rhy_results_internaltime, tissue=="dermis")$Symbol], 
                    wall = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                       filter(rhy_results_walltime, tissue=="dermis")$Symbol])
yrhy_epidermis <- list(internal = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                      filter(rhy_results_internaltime, tissue=="epidermis")$Symbol], 
                       wall = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                  filter(rhy_results_walltime, tissue=="epidermis")$Symbol])


suppfig1A_1 <- ggvenn(yrhy_internaltime, fill_color = c("#1B9E77", "#D95F02"), text_size = 3,
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Internal\ntime"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig1A_2 <- ggvenn(yrhy_walltime, fill_color = c("#1B9E77", "#D95F02"), text_size = 3,
                      stroke_size = 0.25, set_name_size = 4, stroke_linetype = "dashed") +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Wall\ntime"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig1A_3 <- ggvenn(yrhy_dermis, fill_color = c("#1B9E77","#1B9E77"), text_size = 3, #stroke_linetype = c("solid", "dashed"), 
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Dermis"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig1A_4 <- ggvenn(yrhy_epidermis, fill_color = c("#D95F02","#D95F02"), text_size=3, #stroke_linetype = c("solid", "dashed"), 
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Epidermis"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold')

suppfig1A <- plot_grid(suppfig1A_1, NULL, suppfig1A_2, NULL, suppfig1A_3, NULL, suppfig1A_4, 
                       rel_heights = c(1, -0.25, 1, -0.25, 1, -0.25, 1),
                       ncol=1, nrow=7, labels=c("A", ""))


#########

# Supplementary Figure 1B: correlation of amplitudes estimated from wall vs internal time analyses
rhy_results <- rbind(rhy_results_internaltime %>% mutate(analysis="internal_time"),
                     rhy_results_walltime %>% mutate(analysis="wall_time"))
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP",
                 "NFIL3", "TEF", "BHLHE40", "BHLHE41")#, "HLF"

toplot <- rhy_results %>% select(Symbol, tissue, amp_value, analysis) %>% spread(analysis, amp_value) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')"),
         FCamp_internal_time = 2^(2*internal_time),
         FCamp_wall_time = 2^(2*wall_time))

highamp_cutoff <- 0.75
suppfig1B <- ggplot(toplot, aes(x=2^(2*internal_time), y=2^(2*wall_time), color=tissue)) +
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=FCamp_wall_time, y=FCamp_internal_time), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3, 
                  segment.color="black", color="black", parse=TRUE) + 
  facet_wrap(~tissue, scales="free") + guides(color=FALSE) +
  scale_y_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  scale_x_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  xlab('Amplitude fold change, internal time') + #xlab(bquote(~log[2]*' (fold amplitude) internal time')) + 
  xlab('Amplitude fold change, wall time') + #ylab(bquote(~log[2]*' (fold amplitude) wall time')) + 
  theme_custom() 

amp_dermis <- filter(toplot, tissue=="dermis") %>% na.omit()
print(paste0(which(amp_dermis$internal_time > amp_dermis$wall_time) %>% length(), "/", dim(amp_dermis)[1], 
             " genes (rhythmic in dermis) have a higher amplitude when estimated from internal time"))

amp_epidermis <- filter(toplot, tissue=="epidermis")
print(paste0(which(amp_epidermis$internal_time > amp_epidermis$wall_time) %>% length(), "/", dim(amp_epidermis)[1], 
             " genes (rhythmic in epidermis) have a higher amplitude when estimated from internal time"))


#########

# Supplementary Figure 1C: correlation of phases estimated from wall vs internal time analyses
toplot <- rhy_results %>% select(Symbol, tissue, phase_value, analysis) %>% 
  #mutate(phase_clock1 = phase_value + 8) %>% 
  #mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>% select(-phase_value) %>%
  mutate(phase_value = phase_value%%24) %>% 
  spread(analysis, phase_value) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) 

suppfig1C <- ggplot(toplot, aes(x=internal_time, y=wall_time, color=tissue)) +
  geom_abline(slope=1, intercept=4, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=internal_time, y=wall_time), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE) + guides(color=FALSE) +
  facet_wrap(~tissue, scales="free") + 
  scale_y_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  #xlab(bquote('time after'*~MSF[sc]*' (h), internal time')) +  ylab(bquote('time after'*~MSF[sc]*' (h), wall time'))  + 
  xlab('Phase (h), internal time') +  ylab('Phase (h), wall time')  + 
  theme_custom() 


#########

# Supplementary Figure 1D: correlation of p values from wall vs internal time analyses
toplot <- rhy_results %>% select(Symbol, tissue, adj_P_Val_D_or_E, analysis) %>% spread(analysis, adj_P_Val_D_or_E) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')"))

suppfig1D <- ggplot(toplot, aes(x=-log10(internal_time), y=-log10(wall_time), color=tissue)) +
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=-log10(internal_time), y=-log10(wall_time)), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE) + guides(color=FALSE) +
  facet_wrap(~tissue, scales="free") + guides(color=FALSE) +
  #scale_y_continuous(limits=c(0.0,0.05), trans='log10') +
  #scale_x_continuous(limits=c(0.26,1.2), trans='log10') +
  xlab(bquote(~-log[10]*' adjusted'~italic(' p')~'value (internal time)')) + 
  ylab(bquote(~-log[10]*' adjusted'~italic(' p')~'value (wall time)')) + theme_custom() 

################
################

sfig1_1 <- plot_grid(suppfig1A, nrow=1, labels = c("", ""))
sfig1_2 <- plot_grid(suppfig1B, NULL, suppfig1C, NULL, suppfig1D, nrow=5, 
                     labels = c("B", "", "C", "", "D"), rel_heights = c(1,0.1,1,0.1,1))
sfig1   <- plot_grid(sfig1_1, sfig1_2, ncol=2)

sfig1 %>% ggsave('figures/suppfig1.pdf', ., width = 11, height = 8.5)
