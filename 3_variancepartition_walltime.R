suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(hgug4112a.db))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr) )
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(TissueEnrich))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidytext))

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
some_rhy <- dplyr::filter(results, pmax(A_D, A_E) > amp_cutoff & adj_P_Val < fdr_cutoff) #some_rhy as in fig1 analysis
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

#--------------------------------
#--------------------------------


# 2. EXECUTE VARIANCE PARTITION MODEL (WITH INTERACTION TERMS)
# ------------------------------------------------------------
info_exp %<>% mutate(MSF_sc_dec = hms(MSF_sc),
                     MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + second(MSF_sc_dec) / 360),2))
form <- ~ (1|subject) + (1|tissue) + (1|time) + (1|time:tissue) + (1|time:subject) #+ (1|MSF_sc_dec)         
        #time here is wall time (continuous vble [internal_time] can't be modeled as random)

varPart <- fitExtractVarPartModel(geneExpr, form, info_exp, showWarnings=FALSE) 
varPart <- varPart %>% as.data.frame() %>%  
  mutate(ProbeName=rownames(.)) %>% inner_join(genes %>% dplyr::select(ProbeName, Symbol, EntrezID)) %>%
  tibble::column_to_rownames("Symbol")
vp <- sortCols(varPart %>% dplyr::select(-ProbeName, -EntrezID))
    #sortCols: orders columns such that the first one is the one with largest (mean) variation, second is second largest, etc

if (!file.exists("visualize/data/variancePartition_full.csv")){
  write.csv(vp %>% tibble::rownames_to_column("Symbol"), "visualize/data/variancePartition_full.csv")}


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

fig2A_1 <- plotPercentBars(most_var %>% filter(var=="Tissue" | var=="Subject") %>% 
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

fig2A_2 <- plotPercentBars(most_var %>% filter(var=="Time" | var=="Residuals") %>% 
                             dplyr::select(-ProbeName, -EntrezID, -var)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="right", aspect.ratio=1.7) + 
  ylab('Percentage of variance explained') + geom_vline(xintercept=c(10.5, 20.5, 30.5)) + 
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "#00798c", "time" = "#edae49", "Residuals" = "grey80",
                               "time:tissue" = "#F78828", "time:subject" = "#39B600"),
                    labels=c("tissue" = "Inter-layer mean variation", "subject" = "Inter-subject mean variation", 
                             "time" = "Common circadian variation", "time:subject" = "Inter-subject circadian variation",
                             "time:tissue" = "Inter-layer circadian variation", "Residuals" = "Residual variation"))

fig2A_3 <- plotPercentBars(most_var %>% filter(var=="Time:Subject" | var=="Time:Tissue") %>% 
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

fig2A <- ggpubr::ggarrange(fig2A_1, NULL, fig2A_2, NULL, fig2A_3, nrow=1, ncol=5, 
                           common.legend=TRUE, legend="right", widths=c(1.09,0.1,1,0.1,0.98)) 



# Learning about the most varying genes in each variable - Suppfig3A-D: GO analysis
# ---------------------------------------------------------------------------------
number_most_vargenes <- 200
most_var <- rbind(most_tissuevar[1:number_most_vargenes,], most_subjvar[1:number_most_vargenes,], 
                  most_resvar [1:number_most_vargenes,], most_timesubjvar[1:number_most_vargenes,],
                  most_timetissvar[1:number_most_vargenes,])

go <- NULL
for (i in unique(most_var$var)){
  go_i <- goana(most_var %>% filter(var==i) %$% EntrezID, 
                universe = varPart$EntrezID) 
  go_i <- go_i %>% top_n(20, wt=-P.DE) %>% mutate(hits=DE*100/N) %>% as.data.frame() %>% mutate(var=i, P.DE=round(P.DE,3))
  go <- rbind(go, go_i)
}

suppfig3C <- ggplot(go %>% filter(var=="Tissue") %>% mutate(P.DE=ifelse(P.DE == 0, 0.00001, P.DE)), 
                    aes(x=-log10(P.DE), y=reorder_within(Term, -log10(P.DE), var), color=Ont, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' q')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(Ont), switch="y") + 
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
                    aes(x=-log10(P.DE), y=reorder_within(Term, -log10(P.DE), var), color=Ont, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' q')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(Ont), switch="y") + 
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
                    aes(x=-log10(P.DE), y=reorder_within(Term, -log10(P.DE), var), color=Ont, size=hits)) + 
  geom_point() + expand_limits(x=c(1,6)) +
  labs(x=bquote(~-log[10]*italic(' q')~'value'), y="GO:BP term", size="Percentage of hits\nfrom each term") + 
  facet_grid(scales="free_y", space="free",rows=vars(Ont), switch="y") + 
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



# Suppl. Figure 3E: Negative control -- time variance of arrhythmic genes expected ~0
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

# Suppl. Figure 3F: positive control -- high time variance of clock genes expected 
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
fig2C <- plotVarPart(vp.D) + 
  theme_custom() + theme(aspect.ratio=0.6, legend.position = "none") +
  ylab("Percentage of\nvariance explained\n(rhythmic genes in dermis)") +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#1B9E77", "Residuals" = "grey80")) +
  scale_x_discrete(labels = c("tissue" = "Inter-layer\nmean\nvariation", "subject" = "Inter-subject\nmean dermis\nvariation", 
                              "time" = "Dermis\ncircadian\nvariation", "Residuals" = "Residual\ndermis\nvariation"))

fig2D <- plotVarPart(vp.E) + 
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

fig2E <- plotPercentBars(most_timevar.D %>% dplyr::select(-ProbeName, -EntrezID)) + 
  theme_custom() + 
  theme(axis.text.y = element_text(face="italic"), legend.position="bottom", aspect.ratio=2.3,
        rect = element_rect(fill="transparent")) + 
  ylab('Percentage of variance explained\n(rhythmic genes in dermis)') + 
  geom_vline(xintercept=dim(common_genes)[1]+0.5) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#1B9E77", "Residuals" = "grey80"),
                    labels = c("tissue" = "Inter-layer\nmean dermis\nvariation", 
                               "subject" = "Inter-subject\nmean dermis\nvariation", 
                               "time" = "Dermis\ncircadian\nvariation", "Residuals" = "Residual\ndermis\nvariation")) 

fig2F <- plotPercentBars(most_timevar.E %>% dplyr::select(-ProbeName, -EntrezID)) + 
  theme_custom() +
  theme(axis.text.y = element_text(face="italic"), legend.position="bottom", aspect.ratio=2.3, 
        rect = element_rect(fill="transparent")) + 
  ylab('Percentage of variance explained\n(rhythmic genes in epidermis)') + 
  geom_vline(xintercept=dim(common_genes)[1]+0.5) +
  scale_fill_manual(values = c("tissue" = "#d1495b", "subject" = "gray48", "time" = "#D95F02", "Residuals" = "grey80"),
                    labels = c("tissue" = "Inter-layer\nmean epid.\nvariation", 
                               "subject" = "Inter-subject\nmean epid.\nvariation", 
                              "time" = "Epidermis\ncircadian\nvariation", "Residuals" = "Residual\nepidermis\nvariation"))

fig2EF <- plot_grid(NULL, fig2E, NULL, fig2F, nrow=1, ncol=4, labels=c("E", "", "F", ""), rel_widths=c(-0.2,1,-0.2,1))


##############################
##############################


# Arrange plots in grid
# ---------------------
fig2_1 <- plot_grid(fig2A, labels=c('A'))
fig2_2 <- plot_grid(NULL, fig2B, NULL, nrow=1, ncol=3, labels=c("B", "", ""), rel_widths=c(0.095,1,1))
fig2_3 <- plot_grid(fig2C, fig2D, nrow=2, ncol=1, labels=c("C", "D"), rel_heights=c(1,1), align='v')
fig2_4 <- plot_grid(fig2_3, fig2EF, ncol=2, nrow=1, rel_widths=c(0.66,1))
fig2   <- plot_grid(fig2_1, fig2_2, fig2_4, nrow=3, ncol=1, rel_heights=c(1.,0.7,1.), align='v')
fig2 %>% ggsave('figures/fig2.pdf', ., width = 11, height = 11.5)


sfig3_1 <- plot_grid(NULL, suppfig3A, NULL, suppfig3B, nrow=1, ncol=4, labels=c("A", "", "B", ""), rel_widths=c(0.1,1,0.1,1))
sfig3_2 <- plot_grid(suppfig3C, nrow=1, ncol=1, labels=c("C"))
sfig3_3 <- plot_grid(NULL, suppfig3D, nrow=1, ncol=2, labels=c("D", ""), rel_widths=c(0.08,1))
sfig3_4 <- plot_grid(NULL, suppfig3E, nrow=1, ncol=2, labels=c("E", ""), rel_widths=c(0.02,1))
sfig3 <- plot_grid(sfig3_1, NULL, sfig3_2, NULL, sfig3_3, NULL, sfig3_4, nrow=7, rel_heights=c(0.6,0.1,0.85,0.1,0.85,0.1,1.1))
sfig3 %>% ggsave('figures/suppfig3.pdf', ., width = 11, height = 22)
