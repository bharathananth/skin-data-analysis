# go to directory of skin-data-analysis-renv and set it as working directory
# note that 0_preana.R should be run before this file (to pre-process microarray gene expression data)
renv::activate('../skin-data-analysis-renv/') 

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
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(tidytext))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(mgsub))
suppressPackageStartupMessages(library(ReactomePA))

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
experiment <- readRDS("results/experiment.rds") %>% # read sample details from column names
  full_join(info_subjects) %>% # we're going to correct wall time (sampling time) to internal time
  dplyr::mutate(MSF_sc_dec = lubridate::hms(MSF_sc)) %>% 
  dplyr::mutate(
    MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + lubridate::second(MSF_sc_dec) / 360),2),
    diff_to_refsubj = MSF_sc_dec - median(MSF_sc_dec),
    internal_time_ref = time - diff_to_refsubj) %>%
  dplyr::mutate(internal_time = time - MSF_sc_dec)
yave <- readRDS("results/rawdata.rds") # read y gene expression data (without outlier removal)

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
  as.data.frame() %>% dplyr::rename(c("tissueD_inphase"="tissueD:inphase", "tissueE_inphase"="tissueE:inphase",
                                      "tissueD_outphase"="tissueD:outphase", "tissueE_outphase"="tissueE:outphase")) %>% 
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
                        # minimum amplitude requirement: amplitude_D or amplitude_E > cutoff:
                        (pmax(sqrt(results$tissueD_inphase^2 + results$tissueD_outphase^2), 
                              sqrt(results$tissueE_inphase^2 + results$tissueE_outphase^2)) > amp_cutoff),]

# analysis of differentially rhythmic genes
diff_rhy_contrast <- limma::makeContrasts(tissueD_inphase - tissueE_inphase, tissueD_outphase - tissueE_outphase, 
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
if (!file.exists("results/results_populationrhy_internaltime.rds")){
  saveRDS(results, file = "results/results_populationrhy_internaltime.rds")
}  

#--------------------------------
#--------------------------------

# 4. READ RESULTS:
# ----------------
results <- readRDS("results/results_populationrhy_internaltime.rds") %>%
  dplyr::mutate(A_D = sqrt(tissueD_inphase^2 + tissueD_outphase^2), #Amp**2 = cos**2 + sin**2 (note the log2 values)
                A_E = sqrt(tissueE_inphase^2 + tissueE_outphase^2),
                phaseD = atan2(tissueD_outphase, tissueD_inphase)*12/pi, #atan2 takes two arguments (y,x), atan takes the angle
                phaseE = atan2(tissueE_outphase, tissueE_inphase)*12/pi)

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
  xlab(bquote(''*~MSF[sc])) + ylab("Age (years)") + theme_custom()


#####

# Fig1C: How does number of rhythmic genes in dermis and epidermis change with FDR
results_dermis <- results %>% 
  dplyr::select(Symbol, adj_P_Val_D_or_E, rhythmic_in_D, ProbeName) %>%
  tidyr::gather(tissue, rhythmic, -adj_P_Val_D_or_E, -Symbol, -ProbeName) %>%
  dplyr::mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis"))
results_epidermis <- results %>% 
  dplyr::select(Symbol, rhythmic_in_E, adj_P_Val_D_or_E, ProbeName) %>%
  tidyr::gather(tissue, rhythmic, -adj_P_Val_D_or_E, -Symbol, -ProbeName) %>%
  dplyr::mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis"))
results_org <- rbind(results_epidermis, results_dermis)
no_dermis <- dim(results_org %>% filter(rhythmic == TRUE & tissue =="dermis"))[1] # #of genes rhythmic in dermis
no_epidermis <- dim(results_org %>% filter(rhythmic == TRUE & tissue =="epidermis"))[1]
no_common <- dim(results_org %>% filter(rhythmic) %>% group_by(Symbol) %>% filter(n()>1) %>% as.data.frame())[1] / 2
no_common_DR <- dim(results %>% filter(rhythmic_in_D, rhythmic_in_E, diff_rhythmic))[1]

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
rhy_results <- results %>% dplyr::filter(rhythmic_in_D | rhythmic_in_E) 

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
  tidyr::gather(tissue, rhythmic, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName,) %>% 
  dplyr::mutate(tissue = ifelse(tissue=="rhythmic_in_D", "dermis", "epidermis")) %>%
  full_join(results %>% dplyr::filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, A_D, A_E) %>% 
              tidyr::gather(tissue, amp_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(tissue = ifelse(tissue=="A_D", "dermis", "epidermis"))) %>%
  full_join(results %>% filter(rhythmic_in_D | rhythmic_in_E) %>% 
              dplyr::select(Symbol, adj_P_Val_D_or_E, diff_rhythmic, AveExpr, ProbeName, phaseD, phaseE) %>% 
              tidyr::gather(tissue, phase_value, -Symbol, -diff_rhythmic, -adj_P_Val_D_or_E, -AveExpr,-ProbeName) %>% 
              dplyr::mutate(tissue = ifelse(tissue=="phaseD", "dermis", "epidermis"))) 
rhy_results <- rhy_results %>% anti_join( rhy_results %>% filter(rhythmic == FALSE & diff_rhythmic==TRUE) ) #remove non_rhy genes

toplot <- yave$E %>% transform(ProbeName = yave$genes$ProbeName, 
                               Symbol = yave$genes$Symbol) %>% as_tibble() %>%
  tidyr::gather(junk, value, -ProbeName, -Symbol) %>%
  tidyr::separate(junk, c("tissuetime","subject"), sep = "_", convert = TRUE) %>%
  tidyr::separate(tissuetime, c("tissue","time"), convert = TRUE, sep = 1) %>%
  full_join(experiment %>% dplyr::select(subject, time, internal_time)) %>%
  inner_join(rhy_results %>% mutate(tissue=ifelse(tissue=="dermis", "D", "E")) ) %>%
  group_by(ProbeName, Symbol, tissue, subject) %>%
  dplyr::mutate(z.score=(value-mean(value)) / sd(value)) %>% as.data.frame() %>%
  dplyr::group_by(time, Symbol, tissue) %>%
  summarise(z.score=mean(z.score))

toplot %<>% as.data.frame() %>% dplyr::mutate(tissue=ifelse(tissue=="D", "dermis", "epidermis")) %>%
  full_join(rhy_results %>% dplyr::select(Symbol, phase_value, tissue)) %>% 
  dplyr::rename(c("phase"="phase_value")) %>% arrange(phase) 
toplot$Symbol_ord <- factor(toplot$Symbol, levels = rev(unique(toplot$Symbol)), ordered=TRUE)

suppfig2A <- ggplot(toplot %>% arrange(phase) %>% as.data.frame(), 
                aes(x=time-8, y=Symbol_ord)) +
  geom_tile(aes(fill=z.score)) + 
  labs(x=bquote('time after'*~MSF[sc]*' (h)'), y="", fill=expression(italic('z')~'score')) + 
  guides(color = "none") +
  facet_wrap(~tissue, scales="free_y", ncol=2, nrow=1) + theme_custom() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        legend.position = "right", legend.title = element_text(color="black")) + 
  scale_fill_distiller(palette = "RdBu", limits=c(-1.75,1.75), breaks = c(1.5, 0, -1.5)) + 
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24), labels = c(0, 4, 8, 12, 16, 20, 24)) 


#####

#Supplementary figure 2B: What happens with the core clock genes in dermis and epidermis?
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP",
                 "NFIL3", "TEF", "BHLHE40", "BHLHE41")#, "HLF")
  #FIG2 from https://www.annualreviews.org/doi/10.1146/annurev-physiol-073109-130051?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed
lims <- as.POSIXct(strptime(c("1970-01-01 00:00:00","1970-01-01 23:59:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
exps <- as.POSIXct(strptime(c("1969-12-31 23:50:00","1970-01-02 00:10:00"), format = "%Y-%m-%d %H:%M", tz="UTC"))
suppfig2B <- ggplot(rhy_results %>% dplyr::filter(Symbol %in% clock_genes) %>%
                      dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')")),                  
                aes(x=phase_value%%24, y=2^(2*amp_value))) +
  geom_point(aes(color=tissue)) +
  coord_polar(start = 0, direction=1) + 
  scale_x_continuous(limits=c(0,24), breaks = c(0, 6, 12, 18), labels = c(0, 6, 12, 18)) +
  facet_wrap(~tissue, ncol=2, nrow=1) +
  geom_segment(aes(x=phase_value%%24, y=0, 
                   xend=phase_value%%24, yend=2^(2*amp_value), color=tissue)) +
  geom_text_repel(aes(label=Symbol_it), max.overlaps=Inf, box.padding=1, size=3.5, point.padding=.5,
                  segment.color="grey70", color="grey50", parse=TRUE) +
  xlab("") + ylab("Amplitude fold change") +  guides(color = "none") + 
  theme_custom() + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())


#####

# Figure 1D: Scatterplot of FC Amp vs acrophase and histograms: all rhythmic genes in D (left) & E (right)
toplot <- rhy_results %>% dplyr::filter(tissue=="dermis") %>% dplyr::mutate(FC_amp = 2^(2*amp_value)) %>%
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')"),
                clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))
fig1D_1 <- ggplot(toplot, aes(x=phase_value%%24, y=FC_amp, color=tissue, shape=clock_gene)) +
  geom_hline(yintercept=2^(2*amp_cutoff), linetype="dashed", "grey30") +
  geom_point(alpha=0.3) +
  geom_point(data = dplyr::filter(toplot, Symbol %in% clock_genes), 
             aes(x=phase_value%%24, y=FC_amp), color="black", shape=8) +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = dplyr::filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.1, point.padding=.5,  size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + 
  scale_y_continuous(limits=c(1,5), breaks=seq(1:4), trans='log2') +
  ylab("Amplitude fold change") + xlab(bquote('time after'*~MSF[sc]*' (h)')) + 
  theme_custom() + theme(aspect.ratio = 0.7) +
  guides(color = "none") + 
  ggtitle("\ndermis") +
  annotate(geom='text', label=paste0(' ', no_dermis, ' genes'), x=-Inf, y=Inf, hjust=0, vjust=1, color="#1B9E77")
fig1D_1 <- ggExtra::ggMarginal(fig1D_1, type = 'histogram', fill="#1B9E77", color="white")

toplot <- rhy_results %>% dplyr::filter(tissue=="epidermis") %>% dplyr::mutate(FC_amp = 2^(2*amp_value)) %>%
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')"),
                clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))
fig1D_2 <- ggplot(toplot, aes(x=phase_value%%24, y=FC_amp, color=tissue, shape=clock_gene)) +
  geom_hline(yintercept=2^(2*amp_cutoff), linetype="dashed", "grey30") +
  geom_point(alpha=0.3, color="#D95F02") +
  geom_point(data = dplyr::filter(toplot, Symbol %in% clock_genes), 
             aes(x=phase_value%%24, y=FC_amp), color="black", shape=8) +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = dplyr::filter(toplot, Symbol %in% clock_genes), aes(label=Symbol_it), 
                  max.overlaps=Inf, box.padding=1.1, point.padding=.5, size=3,
                  segment.color="black", color="black", parse=TRUE)  +
  scale_x_continuous(limits=c(0,24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24) ) + 
  scale_y_continuous(limits=c(1,5), breaks=seq(1:4), trans='log2')+
  ylab("Amplitude fold change") +xlab(bquote('time after'*~MSF[sc]*' (h)')) + 
  theme_custom() + theme(aspect.ratio = 0.7) +
  guides(color = "none") + 
  ggtitle("\nepidermis") +
  annotate(geom='text', label=paste0(' ', no_epidermis, ' genes'), x=-Inf, y=Inf, hjust=0, vjust=1, color="#D95F02")
fig1D_2 <- ggExtra::ggMarginal(fig1D_2, type = 'histogram', color="white", fill="#D95F02")

fig1D <- plot_grid(fig1D_1, fig1D_2, ncol=2, nrow=1)


#####

# Figure 1E: Correlation of amplitudes of genes rhythmic in BOTH tissues
both_rhy <- rhy_results[which(rhy_results$Symbol %in% rhy_results[duplicated(rhy_results$Symbol),"Symbol"]),] %>% 
  arrange(desc(amp_value)) 
  #check for which symbols are duplicated in the rhy_results dataframe -> that means they are rhythmic in both layers
both_rhy_amp <- both_rhy %>% dplyr::select(Symbol, tissue, amp_value) %>% tidyr::spread(tissue, amp_value) %>% 
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')")) %>% arrange(desc(dermis)) %>%
  dplyr::mutate(clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))
both_rhy_amp$ampFC_dermis <- 2^(2*both_rhy_amp$dermis)
both_rhy_amp$ampFC_epidermis <- 2^(2*both_rhy_amp$epidermis)

# from the genes that are rhythmic in both layers, how many are DIFFERENTIALLY rhythmic?
both_rhy_DR <- results %>% dplyr::filter(diff_rhythmic==TRUE) %>% dplyr::filter(Symbol %in% both_rhy$Symbol)

fig1E <- ggplot(both_rhy_amp, aes(x=ampFC_dermis, y=ampFC_epidermis, shape=clock_gene)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="#7570B3", alpha=0.3) +
  geom_point(data = dplyr::filter(both_rhy_amp, Symbol %in% clock_genes), 
             aes(x=ampFC_dermis, y=ampFC_epidermis), color="black", shape=8) + # #DC0000B2
  geom_point(data=both_rhy_amp %>% dplyr::filter(Symbol %in% both_rhy_DR$Symbol), 
             aes(x=ampFC_dermis, y=ampFC_epidermis, shape=clock_gene), color="coral1") +
  scale_shape_manual(values=c(16,8), guide=FALSE) +
  geom_text_repel(data = dplyr::filter(both_rhy_amp, Symbol %in% clock_genes), box.padding=1., 
                  size=3., aes(x=ampFC_dermis, y=ampFC_epidermis, label=Symbol_it), color="black", 
                  parse=TRUE, point.padding = .5, max.overlaps=Inf) +
  coord_fixed() + theme_bw() + 
  xlab('Amplitude fold change, dermis') + ylab('Amplitude fold change, epidermis') +
  theme_custom() + 
  scale_x_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  scale_y_continuous(limits=c(1, 5), breaks=seq(1:4), trans='log2') +
  annotate(geom='text', label=paste0(' ', no_common, ' genes'), x=1.0, y=5, hjust=0, vjust=1, color="#7570B3") +
  annotate(geom='text', label=paste0(' ', no_common_DR, ' genes'), x=1.00, y=4.2, hjust=0, vjust=1, color="coral1")


#####

# Figure 1F: Correlation of phases of genes rhythmic in BOTH tissues
both_rhy_phase <- both_rhy %>% dplyr::select(Symbol, tissue, phase_value) %>% 
  dplyr::mutate(phase_value = phase_value%%24) %>%
  tidyr::spread(tissue, phase_value) %>%
  dplyr::mutate(Symbol_it = paste0("italic('", Symbol, "')"),
         clock_gene=ifelse(Symbol %in% clock_genes, TRUE, FALSE))

fig1F <- ggplot(both_rhy_phase, aes(x=dermis, y=epidermis, shape=clock_gene)) + 
  geom_abline(slope=1, intercept=0, lty='dashed', color="gray") + 
  geom_point(color="#7570B3", alpha=0.3) +
  geom_point(data = dplyr::filter(both_rhy_phase, Symbol %in% clock_genes), 
             aes(x=dermis, y=epidermis), shape=8, color="black") + 
  geom_point(data=both_rhy_phase %>% dplyr::filter(Symbol %in% both_rhy_DR$Symbol), 
             aes(x=dermis, y=epidermis, shape=clock_gene), color="coral1") +
  scale_shape_manual(values=c(16,8), guide="none") +
  geom_text_repel(data = dplyr::filter(both_rhy_phase, Symbol %in% clock_genes), box.padding=1.2, max.overlaps=Inf, size=3.,
                  aes(x=dermis, y=epidermis, label=Symbol_it), color="black", parse=TRUE, point.padding = .5,
                  xlim = c(0, 24)) +
  coord_fixed() + theme_bw() + xlab("Phase dermis (h)") + ylab("Phase epidermis (h)") +
  theme_custom() +
  scale_x_continuous(breaks=c(0,6,12,18,24), labels=c(0,6,12,18,24)) +
  scale_y_continuous(breaks=c(0,6,12,18,24), labels=c(0,6,12,18,24))


#####

# Figure 1G: Reactome pathway analysis  
# with ReactomePA (but results confirmed with limma, msigdbr, clusterProfiler and other online tools)
rE <- enrichPathway(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="epidermis")$Symbol,]$genes$EntrezID,
                    universe = yave$genes$EntrezID, pvalueCutoff = 0.05, qvalueCutoff = 1.0, minGSSize = 20)

rD <- enrichPathway(yave[yave$genes$Symbol %in% filter(rhy_results, tissue=="dermis")$Symbol,]$genes$EntrezID,
                    universe = yave$genes$EntrezID, pvalueCutoff = 0.05, qvalueCutoff = 1.0, minGSSize = 20)

df_rD <- setReadable(rD, 'org.Hs.eg.db', 'ENTREZID') 
fig1G <- cnetplot(df_rD, color_gene = "#1B9E77", cex_label_gene = .7,layout="kk", foldChange = NULL)
fig1G <- fig1G + theme(legend.position = "none")


rDE <- enrichPathway(gene = filter(results, !is.na(adj_p_val_DR) & diff_rhythmic)$EntrezID, 
                     universe = filter(results, !is.na(adj_p_val_DR))$EntrezID, minGSSize = 20, pvalueCutoff = 0.05, qvalueCutoff = 1.0)
df_rDE <- setReadable(rDE, 'org.Hs.eg.db', 'ENTREZID') %>% as.data.frame() %>%
  tidyr::separate(GeneRatio, c("DE","junk"), convert = TRUE, sep = "/") %>% 
  tidyr::separate(BgRatio, c("N","junk"), convert = TRUE, sep = "/") %>%
  dplyr::mutate(hits=DE*100/N) %>%
  dplyr::mutate(Description = Description %>% tolower())
fig1H <- ggplot(df_rDE, aes(x=-log10(p.adjust), y=reorder(Description, -log10(p.adjust)), size=DE)) + 
  geom_point( color='coral1') +  
  expand_limits(x=c(0,2.0)) + 
  labs(x=bquote(~-log[10]*' adj.'~italic('p')~'value'), y="Pathway", size="Percentage\nof hits") + 
  guides(color = "none") +
  theme_custom() + scale_y_reordered() +
  theme(aspect.ratio=.35, legend.position = "right", legend.title = element_text(color="black"),
        panel.grid.major = element_line(), panel.grid.minor = element_line()) 


#####

#Supplementary figure 2C: Phase Set Enrichment Analysis
## in a terminal, execute the .jar PSEA file: > java -jar PSEA-master/PSEA1.1_VectorGraphics.jar 
## (go to directory where PSEA is)
## PSEA parameter choice: (0, 24, 5, 10000, 0.05) -> gene sets downloaded from https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
rhy_D_or_E <- results[(results$adj_P_Val_D_or_E < fdr_cutoff) & # fdr < cutoff
                        (pmax(results$A_D, results$A_E) > amp_cutoff),]
rhy_D_or_E$rhythmic_in_D <- ifelse(rhy_D_or_E$diff_rhythmic==FALSE, TRUE, 
                                   ifelse(rhy_D_or_E$diff_rhythmic==TRUE & rhy_D_or_E$A_D > amp_cutoff, TRUE, FALSE))
rhy_D_or_E$rhythmic_in_E <- ifelse(rhy_D_or_E$diff_rhythmic==FALSE, TRUE, 
                                   ifelse(rhy_D_or_E$diff_rhythmic==TRUE & rhy_D_or_E$A_E > amp_cutoff, TRUE, FALSE))

# save for PSEA analysis
dir.create("visualize/data/PSEA",showWarnings = FALSE)

if (!file.exists("visualize/data/PSEA/phasesforPSEA_fig1_D.txt")){ 
  rhy_D_or_E %>% dplyr::filter(rhythmic_in_D) %>%
    dplyr::mutate(phase_value = phaseD%%24) %>% 
    dplyr::select(Symbol, phase_value) %>% dplyr::mutate(phase_value=round(phase_value, 0) %>% as.numeric()) %>%
    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_D.txt", 
                row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
}
if (!file.exists("visualize/data/PSEA/phasesforPSEA_fig1_E.txt")){ 
  rhy_D_or_E %>% dplyr::filter(rhythmic_in_E) %>%
    dplyr::mutate(phase_value = phaseE%%24) %>% 
    dplyr::select(Symbol, phase_value) %>% dplyr::mutate(phase_value=round(phase_value, 0) %>% as.numeric()) %>%
    write.table(file = "visualize/data/PSEA/phasesforPSEA_fig1_E.txt", 
                row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
}


# PSEA with Reactome pathways
q_cutoff <- 0.05 

PSEA_d <- read.csv("visualize/data/PSEA/dermis_C2all/results.txt", sep='\t') %>% dplyr::filter(Set.N >= 5)
PSEA_d <- PSEA_d[PSEA_d[,6] <  q_cutoff, ] %>% mutate(tissue = "dermis")

PSEA_e <- read.csv("visualize/data/PSEA/epidermis_C2all/results.txt", sep='\t')%>% dplyr::filter(Set.N >= 5)
PSEA_e <- PSEA_e[PSEA_e[,6] <  q_cutoff, ] %>% mutate(tissue = "epidermis")

PSEA <- full_join(PSEA_d %>% select(Set.ID, Set.N, Vector.average.value, tissue), 
                  PSEA_e %>% select(Set.ID, Set.N, Vector.average.value, tissue)) %>%
  dplyr::mutate(term = gsub("_", " ", Set.ID)) %>% mutate(term=gsub("REACTOME", "", term)) %>% 
  dplyr::mutate(term = term %>% tolower()) %>%
  dplyr::mutate(term = mgsub::mgsub(term, c("rna ", "rna", "dna ", "tgf ", "nod ", "jak stat ", "ecm ", " i ", "ub ", "ii ",
                                            "tnfr2 ", "nf kb ", "dag ", "ip3 ", "ca2", "igf ", "igfbps", "netrin ", "upr" ),
                                    c("RNA ", "RNA", "DNA ", "TGF ", "NOD ", "JAK STAT ","ECM ", " I ", "Ub ", "II ",
                                      "TNFR2 ", "NF KB ", "DAG ", "IP3 ", "Ca2+", "IGF ", "IGFBPs", "Netrin ", "UPR")))

m_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") #https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
df_PSEA <- m_c2 %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame() %>% 
  inner_join(PSEA %>% rename(c("gs_name"="Set.ID"))) %>%
  inner_join(rhy_results %>% rename(c("gene_symbol"="Symbol"))) %>% #see which rhythmic genes belong to which category
  select(-gs_name, - ProbeName, -amp_value, -AveExpr) %>%
  mutate(phase_clock1 = phase_value%%24) %>% select(-phase_value)

suppfig2C <- ggplot(df_PSEA) + 
  geom_point(aes(x=phase_clock1, y=reorder(term, Vector.average.value), color=tissue), alpha=0.5, size=1, shape=4) + 
  geom_point(aes(x=Vector.average.value, y=reorder(term, Vector.average.value), 
                 fill=tissue, size=Set.N), color="black", shape=21) + 
  facet_grid(scales="free_y", space="free",rows=vars(tissue), switch="y") +
  xlab("Phase (h)") + ylab("Reactome pathway") + labs(size='# of genes') + guides(color=FALSE, fill=FALSE) +
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
  #scale_size(breaks=c(5,10,15)) + 
  scale_fill_brewer(palette="Dark2")

#################################
#################################


# ARRANGE PLOTS IN GRID
# ---------------------
fig1_1 <- plot_grid(NULL, fig1B, fig1C, ncol=3, nrow=1, labels=c("A", "B", "C"), rel_widths = c(2,0.65,1))
fig1_2 <- plot_grid(fig1D, NULL, 
                    fig1E + ggtitle("\n"), NULL, fig1F + ggtitle("\n"), 
                    ncol=5, labels=c("D", "", "E", "", "F"), rel_widths = c(2.7, 0.02, 0.9, 0.02, 0.86))
fig1H_2 <- plot_grid(fig1H, NULL, nrow=2, ncol=1, rel_heights = c(1,0.8))
fig1_3 <- plot_grid(NULL, fig1G, NULL, fig1H_2, nrow=1, labels=c("G","", "H", ""), rel_widths = c(.05,.8,.05,1))
fig1 <- plot_grid(fig1_1, NULL, fig1_2, NULL, fig1_3, align='v', nrow=5, 
                  rel_heights = c(1.5, 0.05,1.6, 0.05, 1.1))

fig1 %>% ggsave('figures/fig1.pdf', ., width = 11, height = 9.5)

###

suppfig2_1 <- plot_grid(suppfig2A, NULL, suppfig2B, ncol=3, nrow=1, labels=c("A", "", "B"), rel_widths = c(1,0.2,1.))
suppfig2_2 <- plot_grid(suppfig2C, ncol=2, nrow=1, labels=c("C", ""), rel_widths = c(1,0.2))
suppfig2 <- plot_grid(suppfig2_1, NULL, suppfig2_2, ncol=1, nrow=3, rel_heights = c(.6,.1,1))
suppfig2 %>% ggsave('figures/suppfig2.pdf', ., width = 11, height = 8)

##########
##########

renv::deactivate()
