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
x <- yave$E %>% as.data.frame() %>% t()
rownames_ord <- paste0(sapply(strsplit(rownames(x), split="_"), "[", 2), '_', sapply(strsplit(rownames(x), split="_"), "[", 1)) 
rownames(x) <- rownames_ord
x <- x[order(rownames(x)), ] %>% as.matrix() #ZeitZeiger takes input in this form


# Crossvalidation internal time
# -----------------------------
registerDoParallel(cores=2)
sumabsv <- c(1, 1.5, 2, 3)
nSpc <- 1:6
nFolds <- length(unique(subject))*2
nObs <- dim(x)[1]
nFeatures <- dim(x)[2]

foldid <- rep(1:nFolds, each=length(unique(time)))[-126] #E32_P109 is outlier so it is removed 

# Zeitzeiger takes time between 0 and 1
time_D <- (experiment %>% arrange(subject) %>% filter(tissue=="D") %$% internal_time)/24 #rearrange time
time_D <- ifelse(time_D>1, time_D-1, time_D)
time_E <- (experiment %>% arrange(subject) %>% filter(tissue=="E") %$% internal_time)/24 
time_E <- ifelse(time_E>1, time_E-1, time_E)
time <- c(time_D, time_E)

fitResultList <- zeitzeigerFitCv(x, time, foldid)

spcResultList <- list()
for (ii in 1:length(sumabsv)) {
  spcResultList[[ii]] <- zeitzeigerSpcCv(fitResultList, sumabsv=sumabsv[ii])}

predResultList <- list()
for (ii in 1:length(sumabsv)) {
  predResultList[[ii]] <- zeitzeigerPredictCv(x, time, foldid, spcResultList[[ii]], nSpc=nSpc)}


# Plot the error for each set of parameter values
# -----------------------------------------------
# Before plotting, we need to reorganize the output, making a data.frame with the information for each prediction
timePredList <- lapply(predResultList, function(a) a$timePred)

cvResult <- data.frame(do.call(rbind, timePredList),
                       timeObs=rep(time, length(sumabsv)),
                       sumabsv=rep(sumabsv, each = length(time)),
                       obs=rep(1:nObs, length(sumabsv)),
                       stringsAsFactors=FALSE)

cvResultGath           <- gather(cvResult, key=nSpc, value=timePred, -obs, -timeObs, -sumabsv)
cvResultGath$nSpc      <- as.integer(sapply(as.character(cvResultGath$nSpc), function(a) substr(a, 2, nchar(a))))
cvResultGath$sumabsv   <- factor(cvResultGath$sumabsv)
cvResultGath$timeError <- getCircDiff(cvResultGath$timePred, cvResultGath$timeObs)

# Now calculate the median absolute error for each set of parameter values
cvResultGathGroup <- cvResultGath %>% group_by(sumabsv, nSpc) %>% 
  dplyr::summarize(medae = median(abs(timeError))) %>% as.data.frame()

suppfig4A <- ggplot(cvResultGathGroup) + theme_custom() +
  geom_point(aes(x=nSpc, y=medae*24, shape=sumabsv, color=sumabsv), size=3) +
  labs(x='Number of SPCs', y='Median absolute error (h)') + 
  scale_y_continuous(limits=c(0.03*24, 0.16*24)) + scale_x_continuous(breaks=c(1,2,3,4,5,6)) + expand_limits(x=c(1,6)) +
  theme(legend.title = element_text(),
        legend.position = "right",
        panel.grid.major = element_line(),
        axis.line=element_line()) + scale_color_brewer(palette="Set1")


# Train a model on the full dataset
# ---------------------------------
sumabsv <- 3

fitResultFinal <- zeitzeigerFit(x, time)
spcResultFinal <- zeitzeigerSpc(fitResultFinal$xFitMean, fitResultFinal$xFitResid, 
                                sumabsv=sumabsv) #sumabsv=3 + 2SPCs gives low MAE (doesn't improve much after)

dfVar <- data.frame(spc = 1:length(spcResultFinal$d), propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))


## Plot the phase portrait of SPCs over time: Figure 3B
## ----------------------------------------------------
#z <- x %*% spcResultFinal$v[, 1:2]
#colnames(z) <- c('SPC1', 'SPC2')
#
#z <- data.frame(z, obs=1:nObs, Time=time, check.names=FALSE) %>% mutate(tissue="dermis") %>% 
#  tibble::rownames_to_column() %>% tidyr::separate(rowname, c("subject","junk"), sep = "_", convert = TRUE) %>% select(-junk)
#
#data.arrow <- data.frame(SPC1_start = max(z %$% SPC1),
#                         SPC2_start = (max(z %$% SPC2)-1),
#                         SPC1_end   = (max(z %$% SPC1)-1),
#                         SPC2_end   = max(z %$% SPC2))
#
#suppfig7D <- ggplot(z) + 
#  geom_point(aes(x=SPC1, y=SPC2, color=as.character(Time)), size=2) + theme_custom() + 
#  scale_color_viridis(discrete=TRUE, option='E') +
#  geom_curve(data=data.arrow, aes(x=SPC1_start, y=SPC2_start, xend=SPC1_end, yend=SPC2_end), 
#             arrow=arrow(length=unit(0.03, "npc")), lineend="round") + #expand_limits(y=-5) +
#  facet_wrap(.~tissue, scales='free', nrow=2) + theme(legend.position="none") +
#  scale_y_continuous(expand = c(0.04, 0.75, 0.04, 0.75)) + scale_x_continuous(expand = c(0.04, 0.75, 0.04, 0.75))


# Plot coefficients of the features (time-telling genes) for the SPCs: Figure 3A
# ------------------------------------------------------------------------------
v <- data.frame(spcResultFinal$v[, 1:2])
colnames(v) <- c('SPC 1', 'SPC 2')

v <- v[apply(v, 1, function(r) any(r != 0)), ]
v[v == 0] <- NA
v <- v[do.call(order, v), ]
v$feature <- rownames(v)

v <- inner_join(v, yave$genes %>% as.data.frame() %>% select(Symbol) %>% mutate(feature=as.character(1:n())))

vGath <- gather(v, key=spc, value=Coefficient, -feature, -Symbol) %>%
  dplyr::mutate(feature = factor(feature, levels = rev(v$feature)),
                Symbol = factor(Symbol, levels = rev(v$Symbol)), 
                tissue = "dermis",
                sign = ifelse(Coefficient < 0, "-", "+"), Symbol_it = paste0("italic('", Symbol, "')")) %>% 
  filter(!is.na(Coefficient))
vGath$Coefficient <- abs(vGath$Coefficient)

# Check which ZeitZeiger genes are also found as highly time-variant genes through variance partition
vp <- read.csv("visualize/data/variancePartition_full.csv") %>% dplyr::select(-X) %>% rename("timetissue"="time.tissue")
vp_time.full <- vp %>% arrange(desc(time)) %>% dplyr::select(Symbol, time) %>% mutate(tissue="dermis") %>% head(20)
vp_timetissue.full <- vp %>% arrange(desc(timetissue)) %>% dplyr::select(Symbol, timetissue) %>% 
  mutate(tissue="dermis") %>% head(20)

vGath <- vGath %>% 
  left_join(vp_time.full) %>% mutate(varPart_time_gene=ifelse(is.na(time), FALSE, TRUE)) %>% dplyr::select(-time) %>%
  left_join(vp_timetissue.full) %>% mutate(varPart_timetissue_gene=ifelse(is.na(timetissue), FALSE, TRUE)) %>% 
  dplyr::select(-timetissue) 

vGath$varPart_gene <- ifelse(vGath$varPart_time_gene == FALSE & vGath$varPart_timetissue_gene == TRUE, "time:tissue", 
                             ifelse(vGath$varPart_time_gene == TRUE  & vGath$varPart_timetissue_gene == FALSE, "time", 
                                    vGath$varPart_gene <- ifelse(vGath$varPart_time_gene == TRUE  & 
                                                                   vGath$varPart_timetissue_gene == TRUE, 
                                                                 "time & time:tissue", FALSE)))

suppfig4B <- ggplot(vGath) + facet_wrap(~spc, scales="free") +
  geom_label(aes(x=spc, y=feature, label=Symbol_it, size=Coefficient, color=sign, fill=varPart_gene), 
             label.size=NA, parse=TRUE) + 
  scale_color_manual(values=c("#d11141", "steelblue3")) + 
  scale_fill_manual(values=c("transparent", "khaki1", "#F39B7FB2", "#00B159")) +
  theme_custom() + labs(size="SPC coefficient\n(absolute)") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position="right", 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.title = element_text(face="bold"),
        strip.background = element_rect(fill=alpha("grey", 0.5)),
        strip.text = element_text(size=16),
        aspect.ratio=2.5) + ggtitle(paste0("sumabsv=", sumabsv)) +
  scale_size(limits = c(NA, NA), range = c(3, 8)) + guides(fill=FALSE)
# https://stackoverflow.com/questions/63393553/color-legend-key-labels-with-r-ggplot2-and-remove-the-keys


###################
###################


# Arrange plots in grid
# ---------------------
sfig4 <- plot_grid(suppfig4A, NULL, suppfig4B, labels=c("A","B",""), ncol=3, rel_widths=c(1,0.1,1.5))
sfig4 %>% ggsave('figures/suppfig4.pdf', ., width = 11, height = 7)


###################
###################

