# go to directory where the renv is located and set it as working directory
# note that 0_preana.R should be run before this file (to pre-process microarray gene expression data)

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
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggtext))

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
info_subjects <- read.csv("data/info_subjects_short.csv") %>% dplyr::select(-X) # read info of subjects
experiment <- readRDS("data/experiment.rds") %>% full_join(info_subjects) # read sample details from column names
yave <- readRDS("data/rawdata.rds") # read y gene expression data (without outlier removal)

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
layer  <- factor(experiment$layer)
time    <- experiment$time
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design matrix
design <- model.matrix(~ 0 + layer + layer:inphase + layer:outphase) %>% #H0: rhythms are different across layers
  as.data.frame() %>% dplyr::rename(c("layerD_inphase"="layerD:inphase", "layerE_inphase"="layerE:inphase",
                                      "layerD_outphase"="layerD:outphase", "layerE_outphase"="layerE:outphase")) %>% 
  as.matrix()

# duplicate Correlations
dupcor <- duplicateCorrelation(yave, design, block=subject)

# fits
fit <- limma::lmFit(yave, design, weights = wts, block=subject, correlation=dupcor$consensus)
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# test null hypothesis
rhy_indices <- which(grepl("phase",colnames(design)))
results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
  set_colnames(gsub("\\.","_", colnames(.))) 
results$adj_P_Val_D_or_E <- results$adj_P_Val
rhy_D_or_E <- results[(results$adj_P_Val_D_or_E < fdr_cutoff) & # fdr < cutoff
                        (pmax(sqrt(results$layerD_inphase^2 + results$layerD_outphase^2), 
                              #minimum amplitude requirement: amplitude_D or amplitude_E > cutoff
                              sqrt(results$layerE_inphase^2 + results$layerE_outphase^2)) > amp_cutoff),]

# analysis of differentially rhythmic genes
diff_rhy_contrast <- limma::makeContrasts(layerD_inphase - layerE_inphase, layerD_outphase - layerE_outphase, 
                                          levels = design)
diff_rhy_fit <- limma::contrasts.fit(fit, diff_rhy_contrast)
diff_rhy_fit <- limma::eBayes(diff_rhy_fit, robust = TRUE, trend = TRUE)
diff_rhy_results <- limma::topTable(diff_rhy_fit, number = Inf, sort.by = "none")
diff_rhy_results <- diff_rhy_results[rhy_D_or_E$ProbeName, ]

# organize results
rhy_D_or_E$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$P.Value, method = "BH")
rhy_D_or_E$diff_rhythmic <- rhy_D_or_E$adj_p_val_DR < fdr_cutoff

results <- full_join(results, rhy_D_or_E) %>% dplyr::select(-adj_P_Val)

# save results of fits
if (!file.exists("data/results_populationrhy_walltime.rds")){
  saveRDS(results, file = "data/results_populationrhy_walltime.rds")
}  

#--------------------------------
#--------------------------------

# 4. READ RESULTS: 
# ----------------
# Compare the results from the analysis done with internal time vs. the analysis done with wall time: load results from fig1 
results_internaltime <- readRDS("data/results_populationrhy_internaltime.rds") %>%
  dplyr::mutate(A_D = sqrt(layerD_inphase^2 + layerD_outphase^2), #note log2 values!
                A_E = sqrt(layerE_inphase^2 + layerE_outphase^2),
                phaseD = atan2(layerD_outphase, layerD_inphase)*12/pi,
                phaseE = atan2(layerE_outphase, layerE_inphase)*12/pi)
results_internaltime$rhythmic_in_D <- ifelse(results_internaltime$diff_rhythmic==FALSE, TRUE, 
                                             ifelse(results_internaltime$diff_rhythmic==TRUE & 
                                                      results_internaltime$A_D > amp_cutoff, 
                                                    TRUE, FALSE)) 
results_internaltime$rhythmic_in_E <- ifelse(results_internaltime$diff_rhythmic==FALSE, TRUE, 
                                             ifelse(results_internaltime$diff_rhythmic==TRUE & 
                                                      results_internaltime$A_E > amp_cutoff, 
                                                    TRUE, FALSE))

results_walltime <- readRDS("data/results_populationrhy_walltime.rds") %>%
  dplyr::mutate(A_D = sqrt(layerD_inphase^2 + layerD_outphase^2), #note log2 values!
                A_E = sqrt(layerE_inphase^2 + layerE_outphase^2),
                phaseD = atan2(layerD_outphase, layerD_inphase)*12/pi,
                phaseE = atan2(layerE_outphase, layerE_inphase)*12/pi)
results_walltime$rhythmic_in_D <- ifelse(results_walltime$diff_rhythmic==FALSE, TRUE, 
                                         ifelse(results_walltime$diff_rhythmic==TRUE & results_walltime$A_D > amp_cutoff, TRUE, FALSE)) 
results_walltime$rhythmic_in_E <- ifelse(results_walltime$diff_rhythmic==FALSE, TRUE, 
                                         ifelse(results_walltime$diff_rhythmic==TRUE & results_walltime$A_E > amp_cutoff, TRUE, FALSE))
results <- results_walltime


rhy_results_walltime <- results %>% filter(rhythmic_in_D | rhythmic_in_E) %>%
  dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, rhythmic_in_D, rhythmic_in_E) %>% 
  gather(layer, rhythmic, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName,) %>% 
  dplyr::mutate(layer = ifelse(layer=="rhythmic_in_D", "dermis", "epidermis")) %>%
  full_join(results %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, A_D, A_E) %>% 
              gather(layer, amp_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(layer = ifelse(layer=="A_D", "dermis", "epidermis"))) %>%
  full_join(results %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, phaseD, phaseE) %>% 
              gather(layer, phase_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(layer = ifelse(layer=="phaseD", "dermis", "epidermis"))) 
rhy_results_walltime <- rhy_results_walltime %>% 
  anti_join( rhy_results_walltime %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) )

rhy_results_internaltime <- results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>%
  dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, rhythmic_in_D, rhythmic_in_E) %>% 
  gather(layer, rhythmic, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName,) %>% 
  dplyr::mutate(layer = ifelse(layer=="rhythmic_in_D", "dermis", "epidermis")) %>%
  full_join(results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, A_D, A_E) %>% 
              gather(layer, amp_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(layer = ifelse(layer=="A_D", "dermis", "epidermis"))) %>%
  full_join(results_internaltime %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, phaseD, phaseE) %>% 
              gather(layer, phase_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(layer = ifelse(layer=="phaseD", "dermis", "epidermis"))) 
rhy_results_internaltime <- rhy_results_internaltime %>% 
  anti_join( rhy_results_internaltime %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) )


# Supplementary Figure 3A: comparison of rhythmic genes with wall vs internal time
yrhy_internaltime <- list(dermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                          filter(rhy_results_internaltime, layer=="dermis")$Symbol], 
                          epidermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                             filter(rhy_results_internaltime, layer=="epidermis")$Symbol])
yrhy_walltime <- list(dermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                      filter(rhy_results_walltime, layer=="dermis")$Symbol], 
                      epidermis = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                         filter(rhy_results_walltime, layer=="epidermis")$Symbol])
yrhy_dermis <- list(internal = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                       filter(rhy_results_internaltime, layer=="dermis")$Symbol], 
                    wall = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                       filter(rhy_results_walltime, layer=="dermis")$Symbol])
yrhy_epidermis <- list(internal = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                      filter(rhy_results_internaltime, layer=="epidermis")$Symbol], 
                       wall = yave$genes$ProbeName[yave$genes$Symbol %in% 
                                                  filter(rhy_results_walltime, layer=="epidermis")$Symbol])


suppfig3A_1 <- ggvenn(yrhy_internaltime, fill_color = c("#1B9E77", "#D95F02"), text_size = 3,
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Internal\ntime"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig3A_2 <- ggvenn(yrhy_walltime, fill_color = c("#1B9E77", "#D95F02"), text_size = 3,
                      stroke_size = 0.25, set_name_size = 4, stroke_linetype = "dashed") +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Wall\ntime"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig3A_3 <- ggvenn(yrhy_dermis, fill_color = c("#1B9E77","#1B9E77"), text_size = 3, 
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Dermis"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold') 
suppfig3A_4 <- ggvenn(yrhy_epidermis, fill_color = c("#D95F02","#D95F02"), text_size=3, 
                      stroke_size = 0.25, set_name_size = 4) +
  ggtext::geom_richtext(data = tibble(x=-2, y=0, s="Epidermis"), aes(x, y, label = s), 
                        angle=90, fill = NA, label.color = NA, size=5, fontface='bold')

suppfig3A <- plot_grid(suppfig3A_1, NULL, suppfig3A_2, NULL, suppfig3A_3, NULL, suppfig3A_4, 
                       rel_heights = c(1, -0.25, 1, -0.25, 1, -0.25, 1),
                       ncol=1, nrow=7, labels=c("A", ""))


#########

# Supplementary Figure 3B: correlation of amplitudes estimated from wall vs internal time analyses
rhy_results <- rbind(rhy_results_internaltime %>% dplyr::mutate(analysis="internal_time"),
                     rhy_results_walltime %>% dplyr::mutate(analysis="wall_time"))
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP",
                 "NFIL3", "TEF", "BHLHE40", "BHLHE41")#, "HLF"

toplot <- rhy_results %>% dplyr::select(Symbol, layer, amp_value, analysis) %>% spread(analysis, amp_value) %>% 
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')"),
                FCamp_internal_time = 2^(2*internal_time),
                FCamp_wall_time = 2^(2*wall_time))

highamp_cutoff <- 0.75
suppfig3B <- ggplot(toplot, aes(x=2^(2*internal_time), y=2^(2*wall_time), color=layer)) +
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=FCamp_wall_time, y=FCamp_internal_time), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3, 
                  segment.color="black", color="black", parse=TRUE) + 
  facet_wrap(~layer, scales="free") + guides(color=FALSE) +
  scale_y_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  scale_x_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  xlab('Amplitude fold change, internal time') + 
  ylab('Amplitude fold change, wall time') + 
  theme_custom() 


#########

# Supplementary Figure 3C: correlation of phases estimated from wall vs internal time analyses
toplot <- rhy_results %>% dplyr::select(Symbol, layer, phase_value, analysis) %>% 
  dplyr::mutate(phase_value = phase_value%%24) %>% 
  spread(analysis, phase_value) %>% 
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')")) 

suppfig3C <- ggplot(toplot, aes(x=internal_time, y=wall_time, color=layer)) +
  geom_abline(slope=1, intercept=4, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=internal_time, y=wall_time), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE) + guides(color=FALSE) +
  facet_wrap(~layer, scales="free") + 
  scale_y_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  xlab(bquote('Internal time phase [peak time after'*~MSF[sc]*' (h)]')) + ylab("Wall time phase [peak time (h)]")  + 
  theme_custom() 


#########

# Supplementary Figure 3D: correlation of p values from wall vs internal time analyses
toplot <- rhy_results %>% dplyr::select(Symbol, layer, adj_P_Val_D_or_E, analysis) %>% spread(analysis, adj_P_Val_D_or_E) %>% 
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')"))

suppfig3D <- ggplot(toplot, aes(x=-log10(internal_time), y=-log10(wall_time), color=layer)) +
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(alpha=0.3) + 
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=-log10(internal_time), y=-log10(wall_time)), color="black") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1., point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE) + guides(color=FALSE) +
  facet_wrap(~layer, scales="free") + guides(color=FALSE) +
  xlab(bquote(~-log[10]*' adj.'~italic(' p')~'value (internal time)')) + 
  ylab(bquote(~-log[10]*' adj.'~italic(' p')~'value (wall time)')) + theme_custom() 

################
################


# ARRANGE PLOTS IN GRID
# ---------------------
sfig3_1 <- plot_grid(suppfig3A, nrow=1, labels = c("", ""))
sfig3_2 <- plot_grid(suppfig3B, NULL, suppfig3C, NULL, suppfig3D, nrow=5, 
                     labels = c("B", "", "C", "", "D"), rel_heights = c(1,0.1,1,0.1,1))
sfig3   <- plot_grid(sfig3_1, sfig3_2, ncol=2)

sfig3 %>% ggsave('figures/suppfig3.pdf', ., width = 11, height = 8.5)


##########
##########
