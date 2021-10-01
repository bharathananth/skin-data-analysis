suppressPackageStartupMessages(library(hms))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(zeitzeiger)) #https://github.com/hugheylab/zeitzeiger
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))

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


# 1. READ FILES
# -------------
info_subjects <- read.csv("resources/info_subjects_short.csv") %>% dplyr::select(-X) # read info of subjects
experiment <- readRDS("visualize/data/experiment.rds") %>% full_join(info_subjects)  # read sample details from column names
yave <- readRDS("visualize/data/rawdata.rds") # read y gene expression data (without outlier removal)


# Remove outliers in yave
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R

ind <- which(colnames(yave) == PCA_outliers)    
yave <- yave[, -ind] 
experiment <- experiment[-ind,]
nrow(experiment) == ncol(yave)   #input-check
wts <- NULL #since we are removing PCA-outliers, weights from lmFit (performed later) are set to NULL 

tissue  <- factor(experiment$tissue)
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 
time    <- experiment$time
internal_time <- experiment$internal_time


#--------------------------
#--------------------------


# 2. ZEITZEIGER (https://zeitzeiger.hugheylab.org/articles/introduction.html#load-the-necessary-packages-1)
# -------------

# Dataframe preparation for Zeitzeiger: cols = genes, rows = all observations (across all subjs), separately for each tissue
xD <- yave$E %>% as.data.frame() %>% select(contains("D")) %>% t()
rownames_ord <- paste0(sapply(strsplit(rownames(xD), split="_"), "[", 2), '_', sapply(strsplit(rownames(xD), split="_"), "[", 1)) 
rownames(xD) <- rownames_ord
xD <- xD[order(rownames(xD)), ] %>% as.matrix() #ZeitZeiger takes input in this form

xE <- yave$E %>% as.data.frame() %>% select(contains("E")) %>% t()
rownames_ord <- paste0(sapply(strsplit(rownames(xE), split="_"), "[", 2), '_', sapply(strsplit(rownames(xE), split="_"), "[", 1)) 
rownames(xE) <- rownames_ord
xE <- xE[order(rownames(xE)), ]%>% as.matrix()


# Crossvalidation internal time
# -----------------------------
registerDoParallel(cores=2)
sumabsv <- c(1, 1.5, 2, 3)
nSpc <- 1:6
nFolds <- length(unique(subject))
nObs_D <- dim(xD)[1]; nObs_E <- dim(xE)[1]
nFeatures_D <- dim(xD)[2]; nFeatures_E <- dim(xE)[2]

foldid_D <- rep(1:nFolds, each=length(unique(time)))
foldid_E <- rep(1:nFolds, each=length(unique(time)))[-49] #E32_P109 is outlier so it is removed 

# Zeitzeiger takes time between 0 and 1
time_D <- (experiment %>% arrange(subject) %>% filter(tissue=="D") %$% internal_time)/24 #rearrange time
time_D <- ifelse(time_D>1, time_D-1, time_D)
time_E <- (experiment %>% arrange(subject) %>% filter(tissue=="E") %$% internal_time)/24 
time_E <- ifelse(time_E>1, time_E-1, time_E)

fitResultList_D <- zeitzeigerFitCv(xD, time_D, foldid_D)
fitResultList_E <- zeitzeigerFitCv(xE, time_E, foldid_E)

spcResultList_D <- list()
spcResultList_E <- list()
for (ii in 1:length(sumabsv)) {
  spcResultList_D[[ii]] <- zeitzeigerSpcCv(fitResultList_D, sumabsv=sumabsv[ii])
  spcResultList_E[[ii]] <- zeitzeigerSpcCv(fitResultList_E, sumabsv=sumabsv[ii])}

predResultList_D <- list()
predResultList_E <- list()
for (ii in 1:length(sumabsv)) {
  predResultList_D[[ii]] <- zeitzeigerPredictCv(xD, time_D, foldid_D, spcResultList_D[[ii]], nSpc=nSpc)
  predResultList_E[[ii]] <- zeitzeigerPredictCv(xE, time_E, foldid_E, spcResultList_E[[ii]], nSpc=nSpc)}


# Plot the error for each set of parameter values
# -----------------------------------------------
# Before plotting, we need to reorganize the output, making a data.frame with the information for each prediction
timePredList_D <- lapply(predResultList_D, function(a) a$timePred)
timePredList_E <- lapply(predResultList_E, function(a) a$timePred)

cvResult_D <- data.frame(do.call(rbind, timePredList_D),
                         timeObs=rep(time_D, length(sumabsv)),
                         sumabsv=rep(sumabsv, each = length(time_D)),
                         obs=rep(1:nObs_D, length(sumabsv)),
                         stringsAsFactors=FALSE)
cvResult_E <- data.frame(do.call(rbind, timePredList_E),
                         timeObs=rep(time_E, length(sumabsv)),
                         sumabsv=rep(sumabsv, each = length(time_E)),
                         obs=rep(1:nObs_E, length(sumabsv)),
                         stringsAsFactors=FALSE)

cvResultGath_D           <- gather(cvResult_D, key=nSpc, value=timePred, -obs, -timeObs, -sumabsv)
cvResultGath_D$nSpc      <- as.integer(sapply(as.character(cvResultGath_D$nSpc), function(a) substr(a, 2, nchar(a))))
cvResultGath_D$sumabsv   <- factor(cvResultGath_D$sumabsv)
cvResultGath_D$timeError <- getCircDiff(cvResultGath_D$timePred, cvResultGath_D$timeObs)

cvResultGath_E           <- gather(cvResult_E, key=nSpc, value=timePred, -obs, -timeObs, -sumabsv)
cvResultGath_E$nSpc      <- as.integer(sapply(as.character(cvResultGath_E$nSpc), function(a) substr(a, 2, nchar(a))))
cvResultGath_E$sumabsv   <- factor(cvResultGath_E$sumabsv)
cvResultGath_E$timeError <- getCircDiff(cvResultGath_E$timePred, cvResultGath_E$timeObs)

# Now calculate the median absolute error for each set of parameter values
cvResultGathGroup_D = cvResultGath_D %>% group_by(sumabsv, nSpc) %>% 
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()
cvResultGathGroup_D[,"tissue"] <- "dermis"
cvResultGathGroup_E = cvResultGath_E %>% group_by(sumabsv, nSpc) %>%
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()
cvResultGathGroup_E[,"tissue"] <- "epidermis"
cvResultGathGroup <- rbind(cvResultGathGroup_D, cvResultGathGroup_E)

suppfig4A <- ggplot(cvResultGathGroup) + facet_wrap(~tissue, scales='free') + theme_custom() +
  geom_point(aes(x=nSpc, y=medae, shape=sumabsv, color=sumabsv), size=3) +
  labs(x='Number of SPCs', y='Median absolute error') + 
  scale_y_continuous(limits=c(0.03, 0.16)) + scale_x_continuous(breaks=c(1,2,3,4,5,6)) + expand_limits(x=c(1,6)) +
  theme(legend.title = element_text(),
        legend.position = "right",
        panel.grid.major = element_line(),
        axis.line=element_line()) + scale_color_brewer(palette="Set1")


# Train a model on the full dataset
# ---------------------------------
sumabsv_D <- 2; sumabsv_E <- 2

fitResultFinal_D <- zeitzeigerFit(xD, time_D)
fitResultFinal_E <- zeitzeigerFit(xE, time_E)
spcResultFinal_D <- zeitzeigerSpc(fitResultFinal_D$xFitMean, fitResultFinal_D$xFitResid, 
                                  sumabsv=sumabsv_D) #sumabsv=2 + 2SPCs gives low MAE in D (doesn't improve much after)
spcResultFinal_E <- zeitzeigerSpc(fitResultFinal_E$xFitMean, fitResultFinal_E$xFitResid, 
                                  sumabsv=sumabsv_E) #sumabsv=2 + 2SPCs gives least MAE in E (doesn't improve much after)

dfVar_D <- data.frame(spc = 1:length(spcResultFinal_D$d), propVar = spcResultFinal_D$d^2 / sum(spcResultFinal_D$d^2))
dfVar_E <- data.frame(spc = 1:length(spcResultFinal_E$d), propVar = spcResultFinal_E$d^2 / sum(spcResultFinal_E$d^2))
dfVar_D[,"tissue"] <- "dermis"; dfVar_E[,"tissue"] <- "epidermis"
dfVar <- rbind(dfVar_D, dfVar_E)

suppfig4B <- ggplot(dfVar) + facet_wrap(~tissue, scales='free') + scale_y_continuous(limits=c(0.0,0.7)) +
  geom_point(aes(x=spc, y=propVar, color=tissue), size=3) +
  scale_x_continuous(breaks=seq(1, 10)) +
  labs(x='Number of SPCs', y='Proportion of\nvariance explained') + 
  theme_custom() + theme(legend.title = element_text(),
                         legend.position = "none",
                         panel.grid.major = element_line(),
                         axis.line=element_line()) #first 3 SPCs explain variance, 
                                                   #but error is not improved much after 2SPCs (see suppfig4A)


# Plot the phase portrait of SPCs over time: Figure 3B
# ----------------------------------------------------
zD <- xD %*% spcResultFinal_D$v[, 1:2]
zE <- xE %*% spcResultFinal_E$v[, 1:2]
colnames(zD) <- c('SPC1', 'SPC2'); colnames(zE) <- colnames(zD)

zD <- data.frame(zD, obs=1:nObs_D, Time=time_D, check.names=FALSE) %>% mutate(tissue="dermis") %>% 
  tibble::rownames_to_column() %>% tidyr::separate(rowname, c("subject","junk"), sep = "_", convert = TRUE) %>% select(-junk)
zE <- data.frame(zE, obs=1:nObs_E, Time=time_E, check.names=FALSE) %>% mutate(tissue="epidermis") %>% 
  tibble::rownames_to_column() %>% tidyr::separate(rowname, c("subject","junk"), sep = "_", convert = TRUE) %>% select(-junk)
z <- rbind(zD, zE)

#zGath_D <- gather(zD, key=SPC, value=Abundance, -obs, -Time, -tissue, -subject)
#zGath_E <- gather(zE, key=SPC, value=Abundance, -obs, -Time, -tissue, -subject)
#zGath_D[,"tissue"] <- "dermis"; zGath_E[,"tissue"] <- "epidermis"
#zGath <- rbind(zGath_D, zGath_E)
#
#ggplot(zGath) + 
#  facet_wrap(tissue~SPC, scales = 'free_y') +
#  geom_point(aes(x = Time, y = Abundance, color = tissue), size = 2, shape = 1) + theme_bw()
##ggplot(zGath) + 
##  facet_grid(SPC~tissue, scales = 'free_y') +
##  geom_point(aes(x = Time, y = Abundance, color = tissue), size = 2, shape = 1) + theme_bw()

data.arrow <- data.frame(SPC1_start = max(z %>% filter(tissue=="dermis") %$% SPC1),
                         SPC2_start = (max(z %>% filter(tissue=="dermis") %$% SPC2)-1),
                         SPC1_end   = (max(z %>% filter(tissue=="dermis") %$% SPC1)-1),
                         SPC2_end   = max(z %>% filter(tissue=="dermis") %$% SPC2),
                         tissue="dermis") %>%
  rbind(data.frame(SPC1_start = max(z %>% filter(tissue=="epidermis") %$% SPC1),
                   SPC2_start = (max(z %>% filter(tissue=="epidermis") %$% SPC2)-1),
                   SPC1_end   = (max(z %>% filter(tissue=="epidermis") %$% SPC1)-1),
                   SPC2_end   = max(z %>% filter(tissue=="epidermis") %$% SPC2),
                   tissue="epidermis"))

fig3B <- ggplot(z) + 
  geom_point(aes(x=SPC1, y=SPC2, color=as.character(Time)), size=2) + theme_custom() + scale_color_viridis(discrete=TRUE) +
  geom_curve(data=data.arrow, aes(x=SPC1_start, y=SPC2_start, xend=SPC1_end, yend=SPC2_end), 
             arrow=arrow(length=unit(0.03, "npc")), lineend="round") + #expand_limits(y=-5) +
  facet_wrap(.~tissue, scales='free', nrow=2) + theme(legend.position="none") +
  scale_y_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) + scale_x_continuous(expand = c(0.04, 0.75, 0.04, 0.75))


# Plot coefficients of the features (time-telling genes) for the SPCs: Figure 3A
# ------------------------------------------------------------------------------
vD <- data.frame(spcResultFinal_D$v[, 1:2])
vE <- data.frame(spcResultFinal_E$v[, 1:2])
colnames(vD) <- c('SPC 1', 'SPC 2'); colnames(vE) <- colnames(vD)

vD <- vD[apply(vD, 1, function(r) any(r != 0)), ]; vE <- vE[apply(vE, 1, function(r) any(r != 0)), ]
vD[vD == 0] <- NA; vE[vE == 0] <- NA
vD <- vD[do.call(order, vD), ]; vE <- vE[do.call(order, vE), ]
vD$feature <- rownames(vD); vE$feature <- rownames(vE)

vD <- inner_join(vD, yave$genes %>% as.data.frame() %>% select(Symbol) %>% mutate(feature=as.character(1:n())))
vE <- inner_join(vE, yave$genes %>% as.data.frame() %>% select(Symbol) %>% mutate(feature=as.character(1:n())))

vGath_D = gather(vD, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vD$feature)),
                Symbol = factor(Symbol, levels = rev(vD$Symbol)), 
                tissue = "dermis")
vGath_E = gather(vE, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vE$feature)),
                Symbol = factor(Symbol, levels = rev(vE$Symbol)), 
                tissue = "epidermis")
vGath = rbind(vGath_D, vGath_E)
vGath <- vGath %>% mutate(sign = ifelse(Coefficient < 0, "-", "+"),
                          Symbol_it = paste0("italic('", Symbol, "')")) %>% filter(!is.na(Coefficient))
vGath$Coefficient <- abs(vGath$Coefficient)

# Check which ZeitZeiger genes are also found as highly time-variant genes through variance partition
vp_time.D <- read.csv("visualize/data/variancePartition_top20_timevariantgenes_dermis.csv") %>% 
  dplyr::select(Symbol, Time) %>% mutate(tissue="dermis")
vp_time.E <- read.csv("visualize/data/variancePartition_top20_timevariantgenes_epidermis.csv") %>% 
  dplyr::select(Symbol, Time) %>% mutate(tissue="epidermis")
vp_time <- rbind(vp_time.D, vp_time.E)

vGath <- vGath %>% left_join(vp_time) %>% mutate(varPart_gene=ifelse(is.na(Time), FALSE, TRUE)) %>% dplyr::select(-Time)


fig3A_1 <- ggplot(vGath %>% filter(tissue=="dermis")) + facet_wrap(~spc, scales="free") +
  geom_label(aes(x=spc, y=feature, label=Symbol_it, size=Coefficient, color=sign, fill=varPart_gene), 
             label.size=NA, parse=TRUE) + 
  scale_color_manual(values=c("#d11141", "steelblue3")) + scale_fill_manual(values=c("transparent", "#ffe62f")) +
  theme_custom() + labs(size="SPC coefficient\n(absolute)") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position="right", 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.title = element_text(face="bold"),
        strip.background = element_rect(fill=alpha("#1B9E77", 0.5)),
        strip.text = element_text(size=16),
        aspect.ratio=2.5) + ggtitle(paste0("dermis, sumabsv=", sumabsv_D)) +
  scale_size(limits = c(NA, NA), range = c(3, 8)) + guides(fill=FALSE)
# https://stackoverflow.com/questions/63393553/color-legend-key-labels-with-r-ggplot2-and-remove-the-keys

fig3A_2 <- ggplot(vGath %>% filter(tissue=="epidermis")) + facet_wrap(~spc, scales="free") +
  geom_label(aes(x=spc, y=feature, label=Symbol_it, size=Coefficient, color=sign, fill=varPart_gene), 
             label.size=NA, parse=TRUE) +
  scale_color_manual(values=c("#d11141", "steelblue3")) + scale_fill_manual(values=c("transparent", "#ffe62f")) + 
  theme_custom() + labs(size="SPC coefficient\n(absolute)") + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position="right", 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.title = element_text(face="bold"),
        strip.background = element_rect(fill=alpha("#D95F02", 0.5)),
        strip.text = element_text(size=16),
        aspect.ratio=2.5) + ggtitle(paste0("epidermis, sumabsv=", sumabsv_E)) +
  scale_size(limits = c(NA, NA), range = c(3, 8)) + guides(fill=FALSE)

fig3A <- ggpubr::ggarrange(fig3A_1, NULL, fig3A_2, nrow=1, ncol=3, common.legend=TRUE, legend="right", widths=c(1.,0.1,1))

## Barplot
#fig3A_1 <- ggplot(vGath %>% filter(tissue=="dermis")) + facet_grid(~ spc) +
#  geom_bar(aes(x = Symbol, y = Coefficient), stat = 'identity', fill="#1B9E77") +
#  labs(x = 'Feature') + coord_flip() +
#  theme_bw() + ggtitle(paste0("sumabsv_dermis=", sumabsv_D)) +
#  theme_custom() + theme(panel.grid.major = element_line(),
#                         axis.text.y = element_text(face="italic"))
#fig3A_2 <- ggplot(vGath %>% filter(tissue=="epidermis")) + facet_grid(~ spc) +
#  geom_bar(aes(x = Symbol, y = Coefficient), stat = 'identity', fill="#D95F02") +
#  labs(x = 'Feature') + coord_flip() +
#  theme_bw() + ggtitle(paste0("sumabsv_epidermis=", sumabsv_E)) +
#  theme_custom() + theme(panel.grid.major = element_line(),
#                         axis.text.y = element_text(face="italic"))
#
#fig3A <- plot_grid(fig3A_1, NULL, fig3A_2, nrow=3, ncol=1, rel_heights=c(1,0.1,1))


# Plot timeseries of time-telling genes: Figure 3C
# ------------------------------------------------
zz.genes_D <- data.frame(Symbol = vGath %>% filter(tissue=="dermis") %$% Symbol %>% unique())
zz.genes_D <- zz.genes_D %>% inner_join(yave$genes %>% as.data.frame()) %>% select(-EntrezID, -ENSEMBL)
zz.genes_E <- data.frame(Symbol = vGath %>% filter(tissue=="epidermis") %$% Symbol %>% unique())
zz.genes_E <- zz.genes_E %>% inner_join(yave$genes %>% as.data.frame()) %>% select(-EntrezID, -ENSEMBL)

yD <- yave$E %>% as.data.frame() %>% filter(rownames(.) %in% zz.genes_D$ProbeName) %>% select(contains("D")) %>%
  tibble::rownames_to_column("ProbeName") %>% full_join(zz.genes_D) %>% select(-ProbeName)
yE <- yave$E %>% as.data.frame() %>% filter(rownames(.) %in% zz.genes_E$ProbeName) %>% select(contains("E"))%>%
  tibble::rownames_to_column("ProbeName") %>% full_join(zz.genes_E) %>% select(-ProbeName)

yGath_D <- yD %>% gather(key, expression, -Symbol) %>%
  tidyr::separate(key, c("tissuetime","subject"), sep = "_", convert = TRUE) %>% 
  tidyr::separate(tissuetime, c("tissue","time"), sep = "(?<=[A-Za-z])(?=[0-9])", convert = TRUE) %>%
  inner_join(experiment %>% select(tissue, time, subject, internal_time))
yGath_E <- yE %>% gather(key, expression, -Symbol) %>%
  tidyr::separate(key, c("tissuetime","subject"), sep = "_", convert = TRUE) %>% 
  tidyr::separate(tissuetime, c("tissue","time"), sep = "(?<=[A-Za-z])(?=[0-9])", convert = TRUE) %>%
  inner_join(experiment %>% select(tissue, time, subject, internal_time))

fig3C_1 <- ggplot(yGath_D) + geom_line(aes(x=internal_time, y=expression, color=subject)) + 
  facet_wrap(~Symbol, scales="free", ncol=5, nrow=5) + xlab("internal time") + ylab(bquote(~log[2]*'expression (normalized)')) +
  theme_custom() + theme(strip.text = element_text(face="bold.italic"),
                         legend.position="right",
                         legend.title=element_text(),
                         strip.background = element_rect(fill=alpha("#1B9E77", 0.5)),
                         panel.border = element_rect(fill="transparent"),
                         axis.line = element_blank(),
                         legend.text=element_blank()) + labs(color="subjects\n1 to 11") +
  scale_color_viridis(discrete=TRUE) + scale_x_continuous(breaks=c(8,20,32)) + expand_limits(x=c(6,34)) 

fig3C_2 <- ggplot(yGath_E) + geom_line(aes(x=internal_time, y=expression, color=subject)) + 
  facet_wrap(~Symbol, scales="free", ncol=5, nrow=5) + xlab("internal time") + ylab(bquote(~log[2]*'expression (normalized)')) +
  theme_custom() + theme(strip.text = element_text(face="bold.italic"),
                         legend.position="right",
                         legend.title=element_text(),
                         strip.background = element_rect(fill=alpha("#D95F02", 0.5)),
                         panel.border = element_rect(fill="transparent"),
                         axis.line = element_blank(),
                         legend.text=element_blank()) + labs(color="subjects\n1 to 11") +
  scale_color_viridis(discrete=TRUE) + scale_x_continuous(breaks=c(8,20,32)) + expand_limits(x=c(6,34)) 

fig3C <- ggpubr::ggarrange(fig3C_1, NULL, fig3C_2, ncol=3, nrow=1, common.legend=TRUE, legend="right", heights=c(1.,0.1,0.6))
#fig3C <- ggpubr::ggarrange(fig3C_1, NULL, fig3C_2, nrow=3, ncol=1, common.legend=TRUE, legend="right", heights=c(1.,0.1,0.6))


###################
###################


# Arrange plots in grid
# ---------------------
fig3_1 <- plot_grid(fig3A, fig3B, labels=c("", "B"), ncol=2, nrow=1, rel_widths=c(2, 1))
fig3 <- plot_grid(fig3_1, NULL, fig3C_1, NULL, fig3C_2, nrow=5, ncol=1, 
                  rel_heights=c(1,0.1,2.,0.1,2.*0.6), labels=c("A", "", "C", ""), align="v", axis="l")
fig3 %>% ggsave('figures/fig3.pdf', ., width = 11, height = 22.5)


sfig4 <- plot_grid(suppfig4A, NULL, suppfig4B, labels=c("A", "", "B"), nrow=1, ncol=3, rel_widths=c(1,0.05,0.85),
                   align="h", axis="b")
sfig4 %>% ggsave('figures/suppfig4.pdf', ., width = 11, height = 3.5)


###################
###################

#TODO 1: supplementary figure with other combinations of #SPCs and sumabsv?
#TODO 2: Zeitzeiger on D and E together?