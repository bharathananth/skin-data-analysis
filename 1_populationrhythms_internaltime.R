"
Note that preana.R should be run before this file, in order to normalize etc the microarray gene expression data
"

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
suppressPackageStartupMessages(library(tidytext))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(mgsub))

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


##################################
##################################


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
experiment <- readRDS("visualize/data/experiment.rds") %>% # read sample details from column names
  full_join(info_subjects) %>% # we're going to correct wall time (sampling time) to internal time
  dplyr::mutate(MSF_sc_dec = hms(MSF_sc)) %>% 
  dplyr::mutate(
    MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + second(MSF_sc_dec) / 360),2),
    diff_to_refsubj = MSF_sc_dec - median(MSF_sc_dec),
    internal_time = time - diff_to_refsubj) 
saveRDS(experiment %>% select(-MSF_sc_dec), "visualize/data/experiment.rds")
yave <- readRDS("visualize/data/rawdata.rds") # read y gene expression data (without outlier removal)
yave_D <- readRDS("visualize/data/rawdata_dermis.rds") # read gene expression data in dermis  "
yave_E <- readRDS("visualize/data/rawdata_epidermis.rds") # read gene expression data in epidermis " 

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
time    <- experiment$internal_time
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
if (!file.exists("visualize/data/results_populationrhy_internaltime.rds")){
  saveRDS(results, file = "visualize/data/results_populationrhy_internaltime.rds")
}  

#--------------------------------
#--------------------------------

# 4. READ RESULTS AND PLAY WITH THEM
# ----------------------------------
results <- readRDS("visualize/data/results_populationrhy_internaltime.rds") %>%
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp**2 = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)
# phase = 0 means that the respective gene peaks at 8AM (time of first sampling)

results$rhythmic_in_D <- ifelse(results$diff_rhythmic==FALSE, TRUE, # if DR=FALSE, the gene is rhythmic in both layers 
                                ifelse(results$diff_rhythmic==TRUE & results$A_D > amp_cutoff, TRUE, FALSE)) 
                                #if DR=TRUE, a gene is rhythmic in a layer when its amplitude is higher than cutoff
results$rhythmic_in_E <- ifelse(results$diff_rhythmic==FALSE, TRUE, 
                                ifelse(results$diff_rhythmic==TRUE & results$A_E > amp_cutoff, TRUE, FALSE))

#--------------------------------
#--------------------------------

# 5. FIGURES
# ----------

# Fig1B: Demographics of data
lims <- as.POSIXct(strptime(c(paste0(Sys.Date(), " 02:15:00"), paste0(Sys.Date(), " 07:30:00")), 
                            format = "%Y-%m-%d %H:%M", tz="UTC"))
fig1B <- ggplot(experiment) +
  geom_point(aes(as.POSIXct(x = MSF_sc, format = "%H:%M", tz="UTC"), age, shape=sex), size=2) + 
  scale_y_continuous(breaks=seq(18,32,2)) +
  scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M", limits=lims) +
  theme(aspect.ratio=1) + 
  xlab("Mid sleep time") + ylab("Age (years)") + theme_custom()


#####

# Fig1C: How does number of rhythmic genes in dermis and epidermis change with FDR
#results_dermis <- results %>% 
#  dplyr::select(Symbol, adj_P_Val_D_or_E, A_D, ProbeName, rhythmic_in_D, rhythmic_in_E) %>%
#  gather(tissue, amp_value, -adj_P_Val_D_or_E, -Symbol, -ProbeName, -rhythmic_in_D, -rhythmic_in_E) %>%
#  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
#results_epidermis <- results %>% 
#  dplyr::select(Symbol, adj_P_Val_D_or_E, A_E, ProbeName, rhythmic_in_D, rhythmic_in_E) %>%
#  gather(tissue, amp_value, -adj_P_Val_D_or_E, -Symbol, -ProbeName, -rhythmic_in_D, -rhythmic_in_E) %>%
#  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
#results_org <- rbind(results_epidermis, results_dermis)
#no_dermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="dermis"))[1] # #of genes passing Amp threshold in dermis
#no_epidermis <- dim(results_org %>% filter(amp_value > amp_cutoff & tissue=="epidermis"))[1]
#no_common <- dim(results_org %>% filter(amp_value > amp_cutoff) %>% group_by(Symbol) %>% filter(n()>1) %>% as.data.frame())[1]
#
#fig1C <- ggplot() + 
#  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
#  stat_ecdf(data = results_org %>% mutate(len = ifelse(tissue=="dermis", no_dermis, no_epidermis)) %>% 
#              filter(amp_value > amp_cutoff),
#            aes(x = ifelse(amp_value > amp_cutoff, adj_P_Val_D_or_E), len=len, y = ..y..*len, color=tissue), geom = "step") +  
#  stat_ecdf(data = results_org %>% filter(amp_value > amp_cutoff) %>% group_by(Symbol) %>% filter(n()>1) %>% 
#              as.data.frame() %>% mutate(len = no_common/2, tissue="both"),
#            aes(x = ifelse(amp_value > amp_cutoff, adj_P_Val_D_or_E), len=len, y = ..y..*len, color=tissue), 
#            geom = "step") + # #8491B4B2  
#  scale_x_continuous(breaks = seq(0.00, 0.10, by=0.05)) + 
#  coord_cartesian(xlim = c(-0.0005, 0.103), ylim = c(-0.02,1500), expand=FALSE)  +
#  xlab("False discovery rate\n") + ylab("Number of rhythmic genes") + theme_custom() +
#  scale_colour_manual(values = c("dermis" = "#1B9E77", "epidermis" = "#D95F02", "both" = "#7570B3"))  

results_dermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val_D_or_E, rhythmic_in_D, ProbeName) %>%
  gather(tissue, rhythmic, -adj_P_Val_D_or_E, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis"))
results_epidermis <- results %>% 
  dplyr::select(Symbol, rhythmic_in_E, adj_P_Val_D_or_E, ProbeName) %>%
  gather(tissue, rhythmic, -adj_P_Val_D_or_E, -Symbol, -ProbeName) %>%
  mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis"))
results_org <- rbind(results_epidermis, results_dermis)
no_dermis <- dim(results_org %>% filter(rhythmic == TRUE & tissue =="dermis"))[1] # #of genes rhythmic in dermis
no_epidermis <- dim(results_org %>% filter(rhythmic == TRUE & tissue =="epidermis"))[1]
no_common <- dim(results_org %>% filter(rhythmic) %>% group_by(Symbol) %>% filter(n()>1) %>% as.data.frame())[1] / 2

fig1C <- ggplot() + 
  geom_vline(aes(xintercept=fdr_cutoff), color="grey", linetype="dashed") +
  stat_ecdf(data = results_org %>% mutate(len = ifelse(tissue=="dermis", no_dermis, no_epidermis)) %>% 
              filter(rhythmic),
            aes(x = ifelse(rhythmic, adj_P_Val_D_or_E), len=len, y = ..y..*len, color=tissue), geom = "step") +  
  stat_ecdf(data = results_org %>% filter(rhythmic) %>% group_by(Symbol) %>% filter(n()>1) %>% 
              as.data.frame() %>% mutate(len = no_common, tissue="both"),
            aes(x = ifelse(rhythmic, adj_P_Val_D_or_E), len=len, y = ..y..*len, color=tissue), 
            geom = "step") + # #8491B4B2  
  scale_x_continuous(breaks = seq(0.00, 0.10, by=0.05)) + 
  coord_cartesian(xlim = c(-0.0005, 0.103), ylim = c(-0.02,1500), expand=FALSE)  +
  xlab("False discovery rate\n") + ylab("Number of rhythmic genes") + theme_custom() +
  scale_colour_manual(values = c("dermis" = "#1B9E77", "epidermis" = "#D95F02", "both" = "#7570B3"))  

# Inset of Fig1C: show number of rhythmic genes in each category for our FDR 
#results_amp <- results %>% 
#  dplyr::select(Symbol, adj_P_Val_D_or_E, A_D, A_E, AveExpr, ProbeName) %>%
#  gather(tissue, amp_value, -adj_P_Val_D_or_E, -Symbol, -AveExpr, -ProbeName) %>%
#  filter(amp_value > amp_cutoff) %>%
#  mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))
#results_phase <- results %>% 
#  filter(pmax(A_E, A_D) > amp_cutoff) %>%
#  dplyr::select(Symbol, adj_P_Val_D_or_E, phaseD, phaseE, ProbeName) %>%
#  gather(tissue, phase_value, -adj_P_Val_D_or_E, -Symbol, -ProbeName) %>%
#  mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))
#results_passAmpcutoff <- inner_join(results_amp, results_phase) 
#  # full_join if we want all clock genes, independently of amp > cutoff
#  # results_passAmpcutoff contains genes with amp>cutoff, but no filtering with fdr has been done yet
#rhy_results <- results_passAmpcutoff %>% filter(adj_P_Val_D_or_E < fdr_cutoff)
#rhy_results_ext <- results %>% filter(rhythmic_in_D == TRUE | rhythmic_in_E == TRUE) #another way of visualizing rhythmic genes
#
#fig1C_inset <- ggplot(data=data.frame(var=c("dermis", "epidermis", "both"),
#                                      n=c(dim(rhy_results %>% filter(tissue=="dermis"))[1],
#                                          dim(rhy_results %>% filter(tissue=="epidermis"))[1],
#                                          dim(results %>% filter(A_D>amp_cutoff & A_E>amp_cutoff & adj_P_Val_D_or_E<fdr_cutoff))[1])),
#                      aes(x=var, y=n, fill=var)) + 
#  geom_bar(stat="identity", width=0.65) + 
#  scale_fill_manual(values = c("dermis" = "#1B9E77", "epidermis" = "#D95F02", "both" = "#7570B3")) + 
#  scale_x_discrete(limits=c("dermis", "epidermis", "both")) +
#  labs(x="", y='') +
#  theme_custom() + theme(legend.position="none", aspect.ratio=1.75,
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         #rect = element_rect(fill = "transparent"),
#                         axis.title.y = element_text(size=10),
#                         axis.text.y = element_text(size=8)) + 
#  scale_y_continuous(limits = c(0, 1300), breaks = c(0, 500, 1000)) 


rhy_results <- results %>% filter(rhythmic_in_D | rhythmic_in_E) 

fig1C_inset <- ggplot(data=data.frame(var=c("dermis", "epidermis", "both"),
                                      n=c(dim(rhy_results %>% filter(rhythmic_in_D))[1],
                                          dim(rhy_results %>% filter(rhythmic_in_E))[1],
                                          dim(rhy_results %>% filter(rhythmic_in_D & rhythmic_in_E))[1])),
                      aes(x=var, y=n, fill=var)) + 
  geom_bar(stat="identity", width=0.65) + 
  scale_fill_manual(values = c("dermis" = "#1B9E77", "epidermis" = "#D95F02", "both" = "#7570B3")) + 
  scale_x_discrete(limits=c("dermis", "epidermis", "both")) +
  labs(x="", y='') +
  theme_custom() + theme(legend.position="none", aspect.ratio=1.75,
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         rect = element_rect(fill = "transparent"),
                         axis.title.y = element_text(size=10),
                         axis.text.y = element_text(size=8)) + 
  scale_y_continuous(limits = c(0, 1500), breaks = c(0, 500, 1000, 1500)) 

fig1C <- ggdraw() + draw_plot(fig1C) + draw_plot(fig1C_inset, x = 0.5722, y = .11, width = .4, height = .4)


#####

# Supplementary figure 2A: Heatmap of circadian rhythmic genes: z scores, acropase-ordered
rhy_results %<>% #arrange dataframe to be used later on
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
rhy_results <- rhy_results %>% anti_join( rhy_results %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) ) #remove non_rhy genes

toplot <- yave$E %>% transform(ProbeName = yave$genes$ProbeName, 
                               Symbol = yave$genes$Symbol) %>% as_tibble() %>%
  tidyr::gather(junk, value, -ProbeName, -Symbol) %>%
  tidyr::separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
  inner_join(rhy_results %>% mutate(tissue=ifelse(tissue=="dermis", "D", "E")) ) %>%
  group_by(ProbeName, Symbol, tissue, subject) %>%
  mutate(z.score=(value-mean(value)) / sd(value)) %>% as.data.frame() %>%
  dplyr::group_by(time, Symbol, tissue) %>%
  summarise(z.score=mean(z.score))

toplot %<>% as.data.frame() %>% mutate(tissue=ifelse(tissue=="D", "dermis", "epidermis")) %>%
  full_join(rhy_results %>% select(Symbol, phase_value, tissue)) %>% 
  rename(c("phase"="phase_value")) %>% arrange(phase) 
toplot$Symbol_ord <- factor(toplot$Symbol, levels = rev(unique(toplot$Symbol)), ordered=TRUE)

suppfig2A <- ggplot(toplot %>% arrange(phase), 
                aes(x=time, y=Symbol_ord)) +
  geom_tile(aes(fill=z.score)) + 
  labs(x="time (h)", y="", fill=expression(italic('z')~'score')) + 
  guides(color = FALSE) +
  facet_wrap(~tissue, scales="free_y", ncol=2, nrow=1) + theme_custom() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        legend.position = "right", legend.title = element_text(color="black")) + 
  scale_fill_distiller(palette = "RdBu", limits=c(-1.75,1.75), breaks = c(1.5, 0, -1.5)) + 
  scale_x_continuous(breaks = c(8, 12, 16, 20, 24, 28, 32), labels = c(8, 12, 16, 20, 24, 28, 32)) 


#####

#Supplementary figure 2B: What happens with the core clock genes in dermis and epidermis?
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")
lims <- as.POSIXct(strptime(c("1970-01-01 00:00:00","1970-01-01 23:59:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
exps <- as.POSIXct(strptime(c("1969-12-31 23:50:00","1970-01-02 00:10:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
suppfig2B <- ggplot(rhy_results %>% filter(Symbol %in% clock_genes) %>%
                    mutate(phase_clock1 = phase_value + 8) %>% 
                    mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
                    mutate(phase_clock = format(.POSIXct(3600*phase_clock1, "UTC"), "%Y-%m-%d %H:%M:%S")) %>%
                    mutate(time_of_day=hms::hms(second(phase_clock),minute(phase_clock),hour(phase_clock))) %>%
                    mutate(Symbol_it = paste0("italic('", Symbol, "')")),                  
                aes(x=phase_clock1, y=amp_value)) +
  geom_point(aes(color=tissue)) +
  coord_polar(start = 0, direction=1) + 
  scale_x_continuous(limits=c(0,24), breaks = c(0, 6, 12, 18), labels = c(0, 6, 12, 18)) +
  scale_y_continuous(breaks = seq(0.2,0.8,0.2)) +
  facet_wrap(~tissue, ncol=2, nrow=1) +
  geom_segment(aes(x=phase_clock1, y=0, 
                   xend=phase_clock1, yend=amp_value, color=tissue)) +
  geom_text_repel(aes(label=Symbol_it), max.overlaps=Inf, box.padding=1, size=3.5, point.padding=.5,
                  segment.color="grey70", color="grey50", parse=TRUE) +
  xlab("") + ylab(bquote(~log[2]*' fold amplitude')) + guides(color = FALSE) + 
  theme_custom() + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())


#####

# Figure 1D: Scatterplot of FC Amp vs acrophase and histograms: all rhythmic genes in D (left) & all rhythmic genes in D (right)
toplot <- rhy_results %>% filter(tissue=="dermis") %>% mutate(FC_amp = 2^(2*amp_value)) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"),
         clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))
fig1D_1 <- ggplot(toplot, aes(x=phase_clock1, y=FC_amp, color=tissue, shape=clock_gene)) +
  geom_point(alpha=0.3) +
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=phase_clock1, y=FC_amp), color="black", shape=8) +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.1, point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + scale_y_continuous(trans='log2') + 
  ylab("Amplitude fold change") + xlab("Phase (h)") + theme_custom() + theme(aspect.ratio = 0.7) +
  guides(color = FALSE) + 
  ggtitle("\ndermis") 
fig1D_1 <- ggExtra::ggMarginal(fig1D_1, type = 'histogram', fill="#1B9E77", color="white")

toplot <- rhy_results %>% filter(tissue=="epidermis") %>% mutate(FC_amp = 2^(2*amp_value)) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
  mutate(Symbol_it = paste0("italic('", Symbol, "')"),
         clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))
fig1D_2 <- ggplot(toplot, aes(x=phase_clock1, y=FC_amp, color=tissue, shape=clock_gene)) +
  geom_point(alpha=0.3, color="#D95F02") +
  geom_point(data = filter(toplot, Symbol %in% clock_genes), aes(x=phase_clock1, y=FC_amp), color="black", shape=8) +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.1, point.padding=.5, size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + scale_y_continuous(trans='log2')+
  ylab("Amplitude fold change") + xlab("Phase (h)") + theme_custom() + theme(aspect.ratio = 0.7) +
  guides(color = FALSE) + 
  ggtitle("\nepidermis") 
fig1D_2 <- ggExtra::ggMarginal(fig1D_2, type = 'histogram', color="white", fill="#D95F02")

fig1D <- plot_grid(fig1D_1, fig1D_2, ncol=2, nrow=1)


#####

# Figure 1E: Correlation of amplitudes of genes rhythmic in BOTH tissues
both_rhy <- rhy_results[which(rhy_results$Symbol %in% rhy_results[duplicated(rhy_results$Symbol),"Symbol"]),] %>% 
  arrange(desc(amp_value)) 
  #we check for which symbols are duplicated in the rhy_results dataframe -> that means they are rhythmic in both layers
both_rhy_amp <- both_rhy %>% dplyr::select(Symbol, tissue, amp_value) %>% spread(tissue, amp_value) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% arrange(desc(dermis)) %>%
  mutate(clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))

# from the genes that are rhythmic in both layers, how many are DIFFERENTIALLY rhythmic?
both_rhy_DR <- results %>% filter(diff_rhythmic==TRUE) %>% filter(Symbol %in% both_rhy$Symbol)

fig1E <- ggplot(both_rhy_amp, aes(x=dermis, y=epidermis, shape=clock_gene)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="#7570B3", alpha=0.3) +
  geom_point(data = filter(both_rhy_amp, Symbol %in% clock_genes), 
             aes(x=dermis, y=epidermis), color="black", shape=8) + # #DC0000B2
  geom_point(data=both_rhy_amp %>% filter(Symbol %in% both_rhy_DR$Symbol), 
             aes(x=dermis, y=epidermis, shape=clock_gene), color="coral1") +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = filter(both_rhy_amp, Symbol %in% clock_genes),  box.padding=1., max.overlaps=10, 
                  size=3., aes(x=dermis, y=epidermis, label=Symbol_it), color="black", parse=TRUE, point.padding = .5) +
  coord_fixed() + theme_bw() + 
  xlab(bquote(~log[2]*' fold amp. dermis')) + ylab(bquote(~log[2]*' fold amp. epidermis')) +
  #xlab("Amplitude dermis") + ylab("Amplitude epidermis") +
  theme_custom() + 
  scale_x_continuous(limits=c(0.26,1.2), breaks = seq(0.2, 1.2, by=0.2), trans='log2') +
  scale_y_continuous(limits=c(0.26,1.2), breaks = seq(0.2, 1.2, by=0.2), trans='log2') 

print(paste0(which(both_rhy_amp$epidermis > both_rhy_amp$dermis) %>% length(), "/", dim(both_rhy_amp)[1], 
             " genes (rhythmic in both tissues) with higher amplitude in epidermis than dermis"))


#####

# Figure 1F: Correlation of phases of genes rhythmic in BOTH tissues
both_rhy_phase <- both_rhy %>% dplyr::select(Symbol, tissue, phase_value) %>% 
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>% 
  select(-phase_value) %>% spread(tissue, phase_clock1) %>% 
  mutate(Symbol_it = paste0("italic('", Symbol, "')"),
         clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))

fig1F <- ggplot(both_rhy_phase, aes(x=dermis, y=epidermis, shape=clock_gene)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="#7570B3", alpha=0.3) +
  geom_point(data = filter(both_rhy_phase, Symbol %in% clock_genes), 
             aes(x=dermis, y=epidermis), shape=8, color="black") + # #DC0000B2
  geom_point(data=both_rhy_phase %>% filter(Symbol %in% both_rhy_DR$Symbol), 
             aes(x=dermis, y=epidermis, shape=clock_gene), color="coral1") +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = filter(both_rhy_phase, Symbol %in% clock_genes), box.padding=1.2, max.overlaps=Inf, size=3.,
                  aes(x=dermis, y=epidermis, label=Symbol_it), color="black", parse=TRUE, point.padding = .5,
                  xlim = c(0, 24)) +
  coord_fixed() + theme_bw() + xlab("Phase dermis (h)") + ylab("Phase epidermis (h)") +
  theme_custom() +
  scale_x_continuous(breaks=c(0,6,12,18,24), labels=c(0,6,12,18,24)) +
  scale_y_continuous(breaks=c(0,6,12,18,24), labels=c(0,6,12,18,24))


#####

# Figure 1G: GO analysis -> top 20 GOBP terms 
# with clusterProfiler (but results confirmed with limma, msigdbr and other online tools)
gD <- enrichGO(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID,
                 universe = yave_D$genes$EntrezID, ont="BP",
                 'org.Hs.eg.db', pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")

gE <- enrichGO(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID,
                 universe = yave_E$genes$EntrezID, ont="BP",
                 'org.Hs.eg.db', pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")

if (!file.exists("visualize/data/enrichment/rhygenes_dermis_inttime.txt")){
  write.table(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$ENSEMBL,
              "visualize/data/enrichment/rhygenes_dermis_inttime.txt", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$ENSEMBL,
              "visualize/data/enrichment/rhygenes_epidermis_inttime.txt", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(yave_D$genes$ENSEMBL, "visualize/data/enrichment/exprgenes_dermis_inttime.txt", 
              sep=',', row.names = FALSE, col.names=FALSE, quote=FALSE)
  write.table(yave_E$genes$ENSEMBL, "visualize/data/enrichment/exprgenes_epidermis_inttime.txt", 
              sep=',', row.names = FALSE, col.names=FALSE, quote=FALSE)}

g <- gE %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="epidermis") %>% as.data.frame() %>%
  rbind(gD %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="dermis") %>% as.data.frame()) %>%
  mutate(qvalue = scales::scientific(qvalue, digits = 3),
         pvalue = scales::scientific(pvalue, digits = 3),
         p.adjust = scales::scientific(p.adjust, digits = 3),
         P.DE=as.numeric(pvalue)) %>% 
  select(-ID, -junk, -Count) %>% rename(c("Term"="Description"))


fig1G <- ggplot(g, aes(x=-log10(P.DE), y=reorder_within(Term, -log10(P.DE), tissue), color=tissue, size=hits)) + 
  geom_point() +  
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(2.0,5)) + 
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  guides(color = FALSE) + 
  theme_custom() + scale_y_reordered() +
  theme(aspect.ratio=3.2, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 


#####

# Supplementary Figure 2F: Heatmaps of Spearman correlation between clock genes in dermis and epidermis
timeseries_cg <- yave[which(yave$genes$Symbol %in% clock_genes),]$E %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  rename(c("ProbeName"="rowname")) %>% gather(tissuetime_subj, value, -ProbeName) %>%
  inner_join(yave$genes %>% as.data.frame() %>% select(ProbeName, Symbol)) %>% select(-ProbeName) %>%
  tidyr::separate(tissuetime_subj, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1)

timeseries_cg_raw <- timeseries_cg
timeseries_cg %<>%
  group_by(tissue, subject, Symbol) %>%
  mutate(value=value-mean(value)) %>% as.data.frame()

data_corrmat <- timeseries_cg %>% 
  spread(Symbol, value)
data_corrmat_raw <- timeseries_cg_raw %>% 
  spread(Symbol, value)

corrmat_D <- data_corrmat %>% filter(tissue=="D") %>% mutate(subjtime = paste(subject, time, sep="_")) %>%
  select(-tissue, -time, -subject) %>% tibble::column_to_rownames('subjtime') %>% cor(., method="spearman") %>% as.matrix()
corrmat_E <- data_corrmat %>% filter(tissue=="E") %>% mutate(subjtime = paste(subject, time, sep="_")) %>%
  select(-tissue, -time, -subject) %>% tibble::column_to_rownames('subjtime') %>% cor(., method="spearman") %>% as.matrix() 

corrmat <- corrmat_D %>% as.data.frame() %>% mutate(tissue="\ndermis") %>% tibble::rownames_to_column() %>% 
  gather(key, value, -tissue,-rowname) %>%
  full_join(corrmat_E %>% as.data.frame() %>% mutate(tissue="\nepidermis") %>% tibble::rownames_to_column() %>% 
              gather(key, value, -tissue,-rowname))
corrmat$key <- factor(corrmat$key, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "CRY1", 
                                            "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP"))
corrmat$rowname<- factor(corrmat$rowname, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "CRY1", 
                                                   "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP"))

my.lines <- data.frame(x=c(.5,6.5), y=c(6.5,.5), 
                       xend=c(13.5,6.5), yend=c(6.5,13.5))
suppfig2F <- ggplot(data = corrmat, aes(x=key, y=rowname, fill=value)) + facet_wrap(~tissue) +
  geom_tile() + theme_bw() + xlab("") + ylab("") + labs(fill=expression("Spearman's"~rho)) +
  scale_fill_distiller(palette = "RdBu", limits=c(-1.05,1.05), breaks = c(1, 0.5, 0, -0.5, -1)) + theme_custom() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face="italic", hjust=1), axis.text.y = element_text(face="italic"),
        legend.title = element_text(vjust = 0.8), legend.position = "right") + 
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F) 


#####

# Supplementary figure 2C: KEGG analysis -> top 20 pathways + with qval<0.05
# with clusterProfiler (also done with limma and msigdbr but less overlap than GOBP enrichment)
kD <- enrichKEGG(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID,
                   universe = yave_D$genes$EntrezID,
                   organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")

kE <- enrichKEGG(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID,
                   universe = yave_E$genes$EntrezID,
                   organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1.0, qvalueCutoff = 1.0, minGSSize = 5) %>%
  as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/")
  
k <- kE %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="epidermis") %>% as.data.frame() %>%
  rbind(kD %>% top_n(20, wt=-pvalue) %>% mutate(hits=DE*100/N, tissue="dermis") %>% as.data.frame()) %>%
  mutate(qvalue = scales::scientific(qvalue, digits = 3),
         pvalue = scales::scientific(pvalue, digits = 3),
         p.adjust = scales::scientific(p.adjust, digits = 3),
         P.DE=as.numeric(pvalue)) %>% 
  select(-ID, -junk, -Count) %>% rename(c("Pathway"="Description"))


suppfig2C <- ggplot(k, aes(x=-log10(P.DE), y=reorder_within(Pathway, -log10(P.DE), tissue), color=tissue, size=hits)) + 
  geom_point() +  
  facet_wrap(~tissue, scales='free_y') + expand_limits(x=c(0.75,4.0)) + 
  labs(x=bquote(~-log[10]*italic(' p')~'value'), y="KEGG pathway", size="Percentage of hits\nfrom each term") + 
  guides(color = FALSE) +
  theme_custom() + scale_y_reordered() +
  theme(aspect.ratio=2.3, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 


#####

# Supplementary figure 2D, E: PSEA
## in a terminal, execute the .jar PSEA file: > java -jar PSEA-master/PSEA1.1_VectorGraphics.jar 
## (go to directory where PSEA is)
## PSEA parameter choice: (0, 24, 5, 10000, 0.05) -> gene sets downloaded from https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
fdr_cutoff_PSEA <- 0.1 #less significant cutoff for PSEA analysis

# Repeat DR analysis with a less stringent fdr cutoff
rhy_D_or_E <- results[(results$adj_P_Val_D_or_E < fdr_cutoff_PSEA) & # fdr < cutoff
                        (pmax(results$A_D, results$A_E) > amp_cutoff),]
diff_rhy_contrast <- limma::makeContrasts(tissueD_inphase - tissueE_inphase, tissueD_outphase - tissueE_outphase, levels = design)
diff_rhy_fit <- limma::contrasts.fit(fit, diff_rhy_contrast)
diff_rhy_fit <- limma::eBayes(diff_rhy_fit, robust = TRUE, trend = TRUE)
diff_rhy_results <- limma::topTable(diff_rhy_fit, number = Inf, sort.by = "none")
diff_rhy_results <- diff_rhy_results[rhy_D_or_E$ProbeName, ]

rhy_D_or_E$adj_p_val_DR <- stats::p.adjust(diff_rhy_results$P.Value, method = "BH")
rhy_D_or_E$diff_rhythmic <- rhy_D_or_E$adj_p_val_DR < fdr_cutoff

# which genes are rhythmic in D, E or both?
rhy_D_or_E$rhythmic_in_D <- ifelse(rhy_D_or_E$diff_rhythmic==FALSE, TRUE, 
                                   ifelse(rhy_D_or_E$diff_rhythmic==TRUE & rhy_D_or_E$A_D > amp_cutoff, TRUE, FALSE))
rhy_D_or_E$rhythmic_in_E <- ifelse(rhy_D_or_E$diff_rhythmic==FALSE, TRUE, 
                                   ifelse(rhy_D_or_E$diff_rhythmic==TRUE & rhy_D_or_E$A_E > amp_cutoff, TRUE, FALSE))

# save for PSEA analysis
if (!file.exists("visualize/data/phases_fig1_D.csv")){ 
 
  rhy_D_or_E %>% filter(rhythmic_in_D) %>%
    mutate(phase_clock1 = phaseD + 8) %>% 
    mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
    select(Symbol, phase_clock1) %>% mutate(phase_clock1=round(phase_clock1, 0) %>% as.numeric()) %>%
    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_D.txt", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
}

if (!file.exists("visualize/data/phases_fig1_E.csv")){ 
  rhy_D_or_E %>% filter(rhythmic_in_E) %>%
    mutate(phase_clock1 = phaseE + 8) %>% 
    mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>%
    select(Symbol, phase_clock1) %>% mutate(phase_clock1=round(phase_clock1, 0) %>% as.numeric()) %>%
    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_E.txt", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
}

# Supplementary figure 2D: PSEA with KEGG
q_cutoff <- 0.25# 0.25 for KEGG, 0.05 for GO:BP

PSEA_d_K <- read.csv("visualize/data/PSEA/dermis_C2all/results.txt", sep='\t') %>% 
  filter(Set.N >= 5)
PSEA_d_K <- PSEA_d_K[PSEA_d_K[,4] <  q_cutoff, ] %>% mutate(tissue = "dermis")
PSEA_e_K <- read.csv("visualize/data/PSEA/epidermis_C2all/results.txt", sep='\t')%>% 
  filter(Set.N >= 5)
PSEA_e_K <- PSEA_e_K[PSEA_e_K[,4] <  q_cutoff, ] %>% mutate(tissue = "epidermis")

PSEA_K <- full_join(PSEA_d_K %>% select(Set.ID, Set.N, Vector.average.value, tissue), 
                    PSEA_e_K %>% select(Set.ID, Set.N, Vector.average.value, tissue)) %>%
  mutate(term = gsub("_", " ", Set.ID)) %>% mutate(term=gsub("KEGG", "", term)) %>% mutate(term = term %>% tolower()) %>%
  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))

m_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
m_kegg <- m_c2 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame() %>% 
  inner_join(PSEA_K %>% rename(c("gs_name"="Set.ID"))) %>%
  inner_join(rhy_results %>% rename(c("gene_symbol"="Symbol"))) %>% #see which rhythmic genes belong to which category
  select(-gs_name, - ProbeName, -amp_value, -AveExpr) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>% select(-phase_value)

suppfig2D <- ggplot(m_kegg) + 
  geom_point(aes(x=phase_clock1, y=reorder(term, Vector.average.value), color=tissue), alpha=0.5, size=1, shape=4) + 
  geom_point(aes(x=Vector.average.value, y=reorder(term, Vector.average.value), 
                 fill=tissue, size=Set.N), color="black", shape=21) + 
  facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  xlab("Phase (h)") + ylab("KEGG pathway") + labs(size='# of genes') + guides(color=FALSE, fill=FALSE) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     plot.title = element_text(hjust = 0.5, size=10, face='bold'),
                     legend.position="bottom",
                     legend.spacing.x = unit(1.0, "mm")) +
  expand_limits(y=0) + 
  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  scale_size(breaks=c(5,10,15)) + scale_fill_brewer(palette="Dark2")

# Supplementary figure 2E: PSEA with GOBP
q_cutoff <- 0.05# 0.25 for KEGG, 0.05 for GO:BP

PSEA_d_G <- read.csv("visualize/data/PSEA/dermis_C5GOBP/results.txt", sep='\t') %>% 
  filter(Set.N >= 5)
PSEA_d_G <- PSEA_d_G[PSEA_d_G[,4] <  q_cutoff, ] %>% mutate(tissue = "dermis")
PSEA_e_G <- read.csv("visualize/data/PSEA/epidermis_C5GOBP/results.txt", sep='\t')%>% 
  filter(Set.N >= 5)
PSEA_e_G <- PSEA_e_G[PSEA_e_G[,4] <  q_cutoff, ] %>% mutate(tissue = "epidermis")

PSEA_G <- full_join(PSEA_d_G %>% select(Set.ID, Set.N, Vector.average.value, tissue), 
                    PSEA_e_G %>% select(Set.ID, Set.N, Vector.average.value, tissue)) %>%
  mutate(term = gsub("_", " ", Set.ID)) %>% mutate(term=gsub("GOBP", "", term)) %>% mutate(term = term %>% tolower()) %>%
  mutate(term = mgsub::mgsub(term, c("rna ", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i "), 
                             c("RNA ", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ")))

m_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
m_gobp <- m_c5 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame() %>% 
  inner_join(PSEA_G %>% rename(c("gs_name"="Set.ID"))) %>%
  inner_join(rhy_results %>% rename(c("gene_symbol"="Symbol")))  %>% select(-gs_name, - ProbeName, -amp_value, -AveExpr) %>%
  mutate(phase_clock1 = phase_value + 8) %>% 
  mutate(phase_clock1 = ifelse(phase_clock1 < 0, phase_clock1 + 24, phase_clock1)) %>% select(-phase_value)

suppfig2E <- ggplot(m_gobp) + 
  geom_point(aes(x=phase_clock1, y=reorder(term, Vector.average.value), color=tissue), alpha=0.5, size=1, shape=4) + 
  geom_point(aes(x=Vector.average.value, y=reorder(term, Vector.average.value), fill=tissue, size=Set.N), 
             color="black", shape=21) + 
  facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  xlab("Phase (h)") + ylab("GO:BP term") + labs(size='# of genes') + guides(color=FALSE, fill=FALSE) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     plot.title = element_text(hjust = 0.5, size=10, face='bold'),
                     legend.position="bottom",
                     legend.spacing.x = unit(1.0, "mm")) +
  expand_limits(y=0) + scale_size(breaks=c(25,75,125)) +
  scale_x_continuous(limits=c(0,24), breaks = seq(0, 24, by=6)) +
  scale_fill_brewer(palette="Dark2")


#################################
#################################


# ARRANGE PLOTS IN GRID
# ---------------------
fig1_1 <- plot_grid(NULL, fig1B, fig1C, ncol=3, nrow=1, labels=c("A", "B", "C"), rel_widths = c(2,0.65,1))
fig1_2 <- plot_grid(fig1D, NULL, 
                    fig1E + ggtitle("\n"), NULL, fig1F + ggtitle("\n"), 
                    ncol=5, labels=c("D", "", "E", "", "F"), rel_widths = c(2.7, 0.02, 0.9, 0.02, 0.86))
fig1_3 <- plot_grid(fig1G, labels="G")
fig1 <- plot_grid(fig1_1, NULL, fig1_2, NULL, fig1_3, align='v', nrow=5, 
                  rel_heights = c(1.5, 0.05,1.6, 0.02, 1.8))

fig1 %>% ggsave('figures/fig1.pdf', ., width = 11, height = 11)

###

sfig2_1 <- plot_grid(suppfig2A, NULL, suppfig2B, ncol=3, nrow=1, labels=c("A", "", "B"), rel_widths = c(1,0.1,0.9))
sfig2_2 <- plot_grid(suppfig2C, labels="C")
sfig2_3 <- plot_grid(NULL, NULL, labels=c("D", "E"), rel_widths = c(1,1.15))
suppfig2D <- plot_grid(suppfig2D, NULL, labels=c("", ""), ncol=1, rel_heights = c(1,1))
sfig2_4 <- plot_grid(suppfig2D, suppfig2E, labels=c("", ""), rel_widths = c(1,1.15))
sfig2 <- plot_grid(sfig2_1, NULL, sfig2_2, sfig2_3, sfig2_4, align='v', nrow=5, 
                   rel_heights = c(1.5, 0.1, 1.8, 0.1, 2.))#, 0.1, 1.5))

sfig2 %>% ggsave('figures/suppfig2.pdf', ., width = 11, height = 12.8)
