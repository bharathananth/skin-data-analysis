# go to directory where the renv is located and set it as working directory
# note that 0_preana.R should be run before this file (to pre-process microarray gene expression data)
renv::activate('./renv/') 

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
suppressPackageStartupMessages(library(variancePartition))


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
info_subjects <- read.csv("data/info_subjects_short.csv") %>% dplyr::select(-X) # read info of subjects
experiment <- readRDS("data/experiment.rds") %>% full_join(info_subjects) %>% # read sample details from column names
  dplyr::mutate(layer = ifelse(layer=="D", "dermis", "epidermis"),
                MSF_sc_dec = lubridate::hms(MSF_sc),
                MSF_sc_dec = round((hour(MSF_sc_dec) + minute(MSF_sc_dec) / 60 + lubridate::second(MSF_sc_dec) / 360),2),
                internal_time = time - MSF_sc_dec)
yave <- readRDS("data/rawdata.rds") # read y gene expression data (without outlier removal)


# Remove outliers in yave
PCA_outliers <- "E32_P109" #see PCA analysis in preana.R

ind <- which(colnames(yave) == PCA_outliers)    
yave <- yave[, -ind] 
experiment <- experiment[-ind,]
nrow(experiment) == ncol(yave)   #input-check
wts <- NULL #since we are removing PCA-outliers, weights from lmFit (performed later) are set to NULL 

layer  <- factor(experiment$layer)
subject <- factor(experiment$subject)
sex     <- factor(experiment$sex)
age     <- experiment$age 
time    <- experiment$time
internal_time <- experiment$internal_time


###########################
###########################


# 2. ZEITZEIGER (https://zeitzeiger.hugheylab.org/articles/introduction.html#load-the-necessary-packages-1)
# -------------
# Dataframe preparation for Zeitzeiger: cols = genes, rows = all observations (across all subjs), separately for each layer
xD <- yave$E %>% as.data.frame() %>% dplyr::select(contains("D")) %>% t()
rownames_ord <- paste0(sapply(strsplit(rownames(xD), split="_"), "[", 2), '_', sapply(strsplit(rownames(xD), split="_"), "[", 1)) 
rownames(xD) <- rownames_ord
xD <- xD[order(rownames(xD)), ] %>% as.matrix() #ZeitZeiger takes input in this form

xE <- yave$E %>% as.data.frame() %>% dplyr::select(contains("E")) %>% t()
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
time_D <- (experiment %>% arrange(subject) %>% filter(layer=="dermis") %$% internal_time)/24 #rearrange time
time_D <- ifelse(time_D>1, time_D-1, time_D)
time_E <- (experiment %>% arrange(subject) %>% filter(layer=="epidermis") %$% internal_time)/24 
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


# Plot the error for each set of parameter values: Supplementary Figure 4A
# ------------------------------------------------------------------------
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
cvResultGathGroup_D[,"layer"] <- "dermis"
cvResultGathGroup_E = cvResultGath_E %>% group_by(sumabsv, nSpc) %>%
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()
cvResultGathGroup_E[,"layer"] <- "epidermis"
cvResultGathGroup <- rbind(cvResultGathGroup_D, cvResultGathGroup_E)

suppfig4A <- ggplot(cvResultGathGroup) + facet_wrap(~layer, scales='free') + theme_custom() +
  geom_point(aes(x=nSpc, y=medae*24, shape=sumabsv, color=layer), alpha=0.8, 
             size=3) +
  labs(x='Number of SPCs', y='Median absolute error (h)') + guides(color="none") +
  scale_y_continuous(limits=c(0.03*24, 0.16*24)) + scale_x_continuous(breaks=c(1,2,3,4,5,6)) + expand_limits(x=c(1,6)) +
  theme(legend.title = element_text(),
        legend.position = "right",
        panel.grid.major = element_line(),
        axis.line=element_line(),
        title = element_text()) 


# Train a model on the full dataset
# ---------------------------------
sumabsv_D <- 1.5; sumabsv_E <- 2

fitResultFinal_D <- zeitzeigerFit(xD, time_D)
fitResultFinal_E <- zeitzeigerFit(xE, time_E)
spcResultFinal_D <- zeitzeigerSpc(fitResultFinal_D$xFitMean, fitResultFinal_D$xFitResid, 
                                  sumabsv=sumabsv_D) #sumabsv=2 + 2SPCs gives low MAE in D (doesn't improve much after)
spcResultFinal_E <- zeitzeigerSpc(fitResultFinal_E$xFitMean, fitResultFinal_E$xFitResid, 
                                  sumabsv=sumabsv_E) #sumabsv=2 + 2SPCs gives least MAE in E (doesn't improve much after)

dfVar_D <- data.frame(spc = 1:length(spcResultFinal_D$d), propVar = spcResultFinal_D$d^2 / sum(spcResultFinal_D$d^2))
dfVar_E <- data.frame(spc = 1:length(spcResultFinal_E$d), propVar = spcResultFinal_E$d^2 / sum(spcResultFinal_E$d^2))
dfVar_D[,"layer"] <- "dermis"; dfVar_E[,"layer"] <- "epidermis"
dfVar <- rbind(dfVar_D, dfVar_E)


# Plot the phase portrait of SPCs over time: Figure 3B
# ----------------------------------------------------
zD <- xD %*% spcResultFinal_D$v[, 1:2]
zE <- xE %*% spcResultFinal_E$v[, 1:2]
colnames(zD) <- c('SPC1', 'SPC2'); colnames(zE) <- colnames(zD)

zD <- data.frame(zD, obs=1:nObs_D, Time=time_D, check.names=FALSE) %>% dplyr::mutate(layer="dermis") %>% 
  tibble::rownames_to_column() %>% tidyr::separate(rowname, c("subject","junk"), sep = "_", convert = TRUE) %>% dplyr::select(-junk)
zE <- data.frame(zE, obs=1:nObs_E, Time=time_E, check.names=FALSE) %>% dplyr::mutate(layer="epidermis") %>% 
  tibble::rownames_to_column() %>% tidyr::separate(rowname, c("subject","junk"), sep = "_", convert = TRUE) %>% dplyr::select(-junk)
z <- rbind(zD, zE)

data.arrow <- data.frame(SPC1_start = max(z %>% filter(layer=="dermis") %$% SPC1),
                         SPC2_start = (max(z %>% filter(layer=="dermis") %$% SPC2)-1),
                         SPC1_end   = (max(z %>% filter(layer=="dermis") %$% SPC1)-1),
                         SPC2_end   = max(z %>% filter(layer=="dermis") %$% SPC2),
                         layer="dermis") %>%
  rbind(data.frame(SPC1_start = max(z %>% filter(layer=="epidermis") %$% SPC1),
                   SPC2_start = (max(z %>% filter(layer=="epidermis") %$% SPC2)-1),
                   SPC1_end   = (max(z %>% filter(layer=="epidermis") %$% SPC1)-1),
                   SPC2_end   = max(z %>% filter(layer=="epidermis") %$% SPC2),
                   layer="epidermis"))

fig3B <- ggplot(z) + 
  geom_point(aes(x=SPC1, y=SPC2, color=as.character(Time)), size=2) + theme_custom() + 
  scale_color_viridis(discrete=TRUE, option='E') +
  geom_curve(data=data.arrow, aes(x=SPC1_start, y=SPC2_start, xend=SPC1_end, yend=SPC2_end), 
             arrow=arrow(length=unit(0.03, "npc")), lineend="round") + #expand_limits(y=-5) +
  facet_wrap(layer~., scales='free', ncol=2) + 
  theme(legend.position="none", 
        strip.background = element_rect(fill=alpha("#1B9E77", 0.5), color="white")) +
  scale_y_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) + scale_x_continuous(expand = c(0.04, 0.75, 0.04, 0.75))


# Plot the phase portrait of SPCs over time, but for each subject: supplementary Figure 4D
# -----------------------------------------------------------------------------------------
suppfig4D_1 <- ggplot(z %>% filter(layer=="dermis")) + 
  geom_point(aes(x=SPC1, y=SPC2, color=as.character(Time)), size=2) + theme_custom() + 
  scale_color_viridis(discrete=TRUE, option='E') +
  geom_curve(data=data.arrow %>% filter(layer=="dermis"), aes(x=SPC1_start, y=SPC2_start, xend=SPC1_end, yend=SPC2_end), 
             arrow=arrow(length=unit(0.03, "npc")), lineend="round") + 
  facet_wrap(subject~., scales='free') + theme(legend.position="none") +
  scale_y_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) + scale_x_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) +
  theme(strip.background = element_rect(fill=alpha("#1B9E77", 0.5))) + xlab("\nSPC1") + ylab("SPC2\n") 

suppfig4D_2 <- ggplot(z %>% filter(layer=="epidermis")) + 
  geom_point(aes(x=SPC1, y=SPC2, color=as.character(Time)), size=2) + theme_custom() + 
  scale_color_viridis(discrete=TRUE, option='E') +
  geom_curve(data=data.arrow %>% filter(layer=="epidermis"), aes(x=SPC1_start, y=SPC2_start, xend=SPC1_end, yend=SPC2_end), 
             arrow=arrow(length=unit(0.03, "npc")), lineend="round") + 
  facet_wrap(subject~., scales='free') + theme(legend.position="none") +
  scale_y_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) + scale_x_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) +
  theme(strip.background = element_rect(fill=alpha("#D95F02", 0.5))) + xlab("\nSPC1") + ylab("SPC2\n") 
  
suppfig4D <- plot_grid(NULL, suppfig4D_1, NULL, suppfig4D_2, ncol=4, rel_widths = c(0.1,1,0.1,1))


# Plot coefficients of the features (time-telling genes) for the SPCs: Figure 3A
# ------------------------------------------------------------------------------
vD <- data.frame(spcResultFinal_D$v[, 1:2])
vE <- data.frame(spcResultFinal_E$v[, 1:2])
colnames(vD) <- c('SPC 1', 'SPC 2'); colnames(vE) <- colnames(vD)

vD <- vD[apply(vD, 1, function(r) any(r != 0)), ]; vE <- vE[apply(vE, 1, function(r) any(r != 0)), ]
vD[vD == 0] <- NA; vE[vE == 0] <- NA
vD <- vD[do.call(order, vD), ]; vE <- vE[do.call(order, vE), ]
vD$feature <- rownames(vD); vE$feature <- rownames(vE)

vD <- inner_join(vD, yave$genes %>% as.data.frame() %>% dplyr::select(Symbol) %>% dplyr::mutate(feature=as.character(1:n())))
vE <- inner_join(vE, yave$genes %>% as.data.frame() %>% dplyr::select(Symbol) %>% dplyr::mutate(feature=as.character(1:n())))

vGath_D = gather(vD, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vD$feature)),
                Symbol = factor(Symbol, levels = rev(vD$Symbol)), 
                layer = "dermis")
vGath_E = gather(vE, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(vE$feature)),
                Symbol = factor(Symbol, levels = rev(vE$Symbol)), 
                layer = "epidermis")
vGath = rbind(vGath_D, vGath_E)
vGath <- vGath %>% dplyr::mutate(sign = ifelse(Coefficient < 0, "-", "+"),
                                 Symbol_it = paste0("italic('", Symbol, "')")) %>% filter(!is.na(Coefficient))
vGath$Coefficient <- abs(vGath$Coefficient)

fig3A_1 <- ggplot(vGath %>% filter(layer=="dermis")) + facet_wrap(~spc, scales="free") +
  geom_label(aes(x=spc, y=feature, label=Symbol_it, size=Coefficient, color=sign), 
             label.size=NA, parse=TRUE) + 
  scale_color_manual(values=c("#d11141", "steelblue3")) + 
  scale_fill_manual(values=c("transparent", "#04bedb", "#fa8e9d", "#4B0082")) +
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
        aspect.ratio=1.5) + ggtitle(paste0("dermis, sumabsv=", sumabsv_D)) +
  scale_size(limits = c(NA, NA), range = c(4, 7)) + guides(fill="none")
# https://stackoverflow.com/questions/63393553/color-legend-key-labels-with-r-ggplot2-and-remove-the-keys

fig3A_2 <- ggplot(vGath %>% filter(layer=="epidermis")) + facet_wrap(~spc, scales="free") +
  geom_label(aes(x=spc, y=feature, label=Symbol_it, size=Coefficient, color=sign), 
             label.size=NA, parse=TRUE) +
  scale_color_manual(values=c("#d11141", "steelblue3")) + 
  scale_fill_manual(values=c("transparent", "#04bedb", "#fa8e9d", "#4B0082")) +
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
        aspect.ratio=1.5) + ggtitle(paste0("epidermis, sumabsv=", sumabsv_E)) +
  scale_size(limits = c(NA, NA), range = c(4, 7)) + guides(fill="none")

fig3A <- ggpubr::ggarrange(fig3A_1, NULL, fig3A_2, nrow=1, ncol=3, common.legend=TRUE, legend="right", widths=c(1.,0.1,1))

ZZ_Wu2018.E <- c("ARNTL", "C2CD4B", "TRIM35", "IFFO2", "GALNT11", "NR1D1", "FKBP5", "HLF", "DBP", "NR1D2", "PER3", "TEF", 
                 "PER1", "CIART", 'PHTF2', "DYNC1LI2", "RGS3", "FANCL", "TMEM168", "METTL3", "KIAA0907", "TTC14", "NHLH2", 
                 "BLCAP", "PLLP", "TSC22D3", "RALBP1", "ELOVL6", "ZBTB16")
ZZ_Wu2020.E <- c("ARNTL", "NR1D2", "TEF", "PER3", "PER1", "CIART", "PER2", "NAMPT", "PAQR3", "FAM117A", "DBP", "NR1D1")
ZZ_Wu2020.D <- c("NR1D1", "CIART", "TEF", "NR1D2", "TSC22D3", "PER3", "DBP", 'PER1', "HLF", "BHLHE41", "KLF9", "FKBP5",
                 "ZBTB16", "ARNTL", "IRAK3", "PER2", "PIK3R1", "KLF15", "TACC1", "DSE", "NPAS2")

ourZZgenes_inWu2018.E <- vGath %>% filter(layer=="epidermis") %>% filter(Symbol %in% ZZ_Wu2018.E)
ourZZgenes_inWu2020.E <- vGath %>% filter(layer=="epidermis") %>% filter(Symbol %in% ZZ_Wu2020.E)
ourZZgenes_inWu2020.D <- vGath %>% filter(layer=="dermis") %>% filter(Symbol %in% ZZ_Wu2020.D) 


# Plot timeseries of time-telling genes: Supplementary Figure 4B
# --------------------------------------------------------------
zz.genes_D <- data.frame(Symbol = vGath %>% filter(layer=="dermis") %$% Symbol %>% unique())
zz.genes_D <- zz.genes_D %>% inner_join(yave$genes %>% as.data.frame()) %>% dplyr::select(-EntrezID, -ENSEMBL)
zz.genes_E <- data.frame(Symbol = vGath %>% filter(layer=="epidermis") %$% Symbol %>% unique())
zz.genes_E <- zz.genes_E %>% inner_join(yave$genes %>% as.data.frame()) %>% dplyr::select(-EntrezID, -ENSEMBL)

yD <- yave$E %>% as.data.frame() %>% filter(rownames(.) %in% zz.genes_D$ProbeName) %>% dplyr::select(contains("D")) %>%
  tibble::rownames_to_column("ProbeName") %>% full_join(zz.genes_D) %>% dplyr::select(-ProbeName)
yE <- yave$E %>% as.data.frame() %>% filter(rownames(.) %in% zz.genes_E$ProbeName) %>% dplyr::select(contains("E"))%>%
  tibble::rownames_to_column("ProbeName") %>% full_join(zz.genes_E) %>% dplyr::select(-ProbeName)

yGath_D <- yD %>% gather(key, expression, -Symbol) %>%
  tidyr::separate(key, c("layertime","subject"), sep = "_", convert = TRUE) %>% 
  tidyr::separate(layertime, c("layer","time"), sep = "(?<=[A-Za-z])(?=[0-9])", convert = TRUE) %>%
  dplyr::mutate(layer=ifelse(layer=="D", "dermis", "epidermis")) %>%
  inner_join(experiment %>% dplyr::select(layer, time, subject, internal_time))
yGath_E <- yE %>% gather(key, expression, -Symbol) %>%
  tidyr::separate(key, c("layertime","subject"), sep = "_", convert = TRUE) %>% 
  tidyr::separate(layertime, c("layer","time"), sep = "(?<=[A-Za-z])(?=[0-9])", convert = TRUE) %>%
  dplyr::mutate(layer=ifelse(layer=="D", "dermis", "epidermis")) %>%
  inner_join(experiment %>% dplyr::select(layer, time, subject, internal_time))

suppfig4B_1 <- ggplot(yGath_D) + geom_line(aes(x=internal_time, y=expression, color=subject)) + 
  facet_wrap(~Symbol, scales="free", ncol=4, nrow=5) + 
  xlab(bquote('time after'*~MSF[sc]*' (h)')) + ylab(bquote(~log[2]*'expression (normalized)')) +
  theme_custom() + theme(strip.text = element_text(face="bold.italic"),
                         legend.position="none",
                         legend.title=element_text(),
                         strip.background = element_rect(fill=alpha("#1B9E77", 0.5)),
                         panel.border = element_rect(fill="transparent"),
                         axis.line = element_blank(),
                         legend.text=element_blank()) + labs(color="subjects\n1 to 11") +
  scale_color_viridis(discrete=TRUE) + scale_x_continuous(breaks=c(0,12,24)) + expand_limits(x=c(0,12,24)) 

suppfig4B_2 <- ggplot(yGath_E) + geom_line(aes(x=internal_time, y=expression, color=subject)) + 
  facet_wrap(~Symbol, scales="free", ncol=4, nrow=5) + 
  xlab(bquote('time after'*~MSF[sc]*' (h)')) + ylab(bquote(~log[2]*'expression (normalized)')) +
  theme_custom() + theme(strip.text = element_text(face="bold.italic"),
                         legend.position="none",
                         legend.title=element_text(),
                         strip.background = element_rect(fill=alpha("#D95F02", 0.5)),
                         panel.border = element_rect(fill="transparent"),
                         axis.line = element_blank(),
                         legend.text=element_blank()) + labs(color="subjects\n1 to 11") +
  scale_color_viridis(discrete=TRUE) + scale_x_continuous(breaks=c(0,12,24)) + expand_limits(x=c(0,12,24)) 

suppfig4B <- ggpubr::ggarrange(suppfig4B_1, NULL, suppfig4B_2, ncol=3, nrow=1, 
                               common.legend=TRUE, legend="right", heights=c(1.,0.1,0.6))


# Where are the ZeitZeiger genes in the scatterplot of variability in amplitude/magnitude/phase across subjects?
# --------------------------------------------------------------------------------------------------------------
vp.full <- read.csv("data/variance_rhythmic_parameters_full.csv") %>% dplyr::select(-X) 
vp.full %<>% filter(Amp>.15)

# Calculate cv from variances -> this is from the full vP, not D vs E
variation_A   <- vp.full %>% dplyr::select(ProbeName, Amp, var_A_subject, var_A_layer) %>% 
  gather(variable, variance, -ProbeName, -Amp) %>%
  dplyr::mutate(variable = ifelse(variable=="var_A_subject", "A_S", "A_T"),
                sd = sqrt(variance),
                cv = sd/Amp) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_phi <- vp.full %>% dplyr::select(ProbeName, phase, var_phi_subject, var_phi_layer) %>%
  gather(variable, variance, -ProbeName, -phase) %>%
  dplyr::mutate(variable = ifelse(variable=="var_phi_subject", "phi_S", "phi_T"),
                sd = sqrt(variance),
                cv = sd/24) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))
variation_magn <- vp.full %>% dplyr::select(ProbeName, magn, var_magn_subject, var_magn_layer) %>%
  gather(variable, variance, -ProbeName, -magn) %>%
  dplyr::mutate(variable = ifelse(variable=="var_magn_subject", "magn_S", "magn_T"),
                sd = sqrt(variance),
                cv = sd/magn) %>% arrange(desc(cv)) %>% inner_join(yave$genes %>% dplyr::select(ProbeName, Symbol))

variation_full <- rbind(variation_A %>% dplyr::rename(c("value_fit"="Amp")), 
                        variation_phi %>% dplyr::rename(c("value_fit"="phase")) %>%
                          dplyr::mutate(value_fit = ifelse(value_fit < 0, value_fit + 24, value_fit))) %>% 
  rbind(variation_magn %>% dplyr::rename(c("value_fit"="magn"))) %>%
  tidyr::separate(variable, c("rhythmic_par","variable"), sep = "_", convert = TRUE) %>%
  dplyr::mutate(rhythmic_par=ifelse(rhythmic_par == "A", "amplitude", 
                                    ifelse(rhythmic_par=="phi", "phase", "magnitude")),
                variable=ifelse(variable == "S", "subject", "layer")) 

variation_full$rhythmic_par <- factor(variation_full$rhythmic_par, levels=c("magnitude", "amplitude", "phase"))
variation_full$effect <- ifelse(variation_full$variable=="layer", "layer", "subject")
variation_full$effect <- factor(variation_full$effect, levels=c("layer", "subject"))
variation_full$Symbol_it <- paste0("italic('", variation_full$Symbol, "')")

df_cv_magn_amp <- variation_full %>% 
  dplyr::select(rhythmic_par, variable, cv, Symbol, Symbol_it) %>% 
  filter(variable=="subject", rhythmic_par != "phase") %>% 
  spread(rhythmic_par, cv) 
df_cv_magn_phi <- variation_full %>% 
  dplyr::select(rhythmic_par, variable, cv, Symbol, Symbol_it) %>% 
  filter(variable=="subject", rhythmic_par != "amplitude") %>% 
  spread(rhythmic_par, cv) 

suppfig4C_1 <- ggplot(df_cv_magn_amp) +
    geom_point(aes(x=magnitude, y=amplitude), color="#00798c", alpha=0.6) +
    geom_point(data=filter(df_cv_magn_amp, Symbol %in% zz.genes_D$Symbol | Symbol %in% zz.genes_E$Symbol),
               aes(x=magnitude, y=amplitude), color="red", alpha=1) +
    labs(x='magnitude CV', y='amplitude CV') #BECN1 and FOCAD missing!!
suppfig4C_2 <- ggplot(df_cv_magn_phi) +
    geom_point(aes(x=magnitude, y=phase), color="#00798c", alpha=0.6) +
    geom_point(data=filter(df_cv_magn_phi, Symbol %in% zz.genes_D$Symbol | Symbol %in% zz.genes_E$Symbol),
               aes(x=magnitude, y=phase), color="red", alpha=1) +
    labs(x='magnitude CV', y='phase CV')

suppfig4C <- plot_grid(NULL, suppfig4C_1, NULL, suppfig4C_2, ncol=4, rel_widths = c(0.1,1,0.1,1))

###################
###################


# Arrange plots in grid
# ---------------------
fig3 <- plot_grid(NULL, fig3A, NULL, fig3B, labels=c("A","","B", ""), ncol=4, nrow=1, rel_widths=c(0.01,1.,0.02,.6))
fig3 %>% ggsave('figures/fig3.pdf', ., width = 11, height = 3.)

sfig4 <- plot_grid(plot_grid(suppfig4A, suppfig4C, nrow=1, rel_widths = c(1,1), labels = c("A", "D")), 
                   plot_grid(suppfig4B_1, NULL, suppfig4B_2, nrow=1, rel_widths = c(1,.1,1), labels = c("B", "", "C")),
                   rel_heights=c(0.75,1.), nrow=2) %T>%
        ggsave('figures/suppfig4.pdf', ., width = 11, height = 8.5)


##########
##########

renv::deactivate()