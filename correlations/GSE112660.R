#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112660

"
I haven't removed background -> do I have to?
How does BA get the RDS files from GEO?
"

setwd("~/Documents/WORK/POSTDOC/projects/skin-data-analysis/correlations/")
library(dplyr) 
library(tidyr)
library(ggplot2)
library(tibble)
library(limma)
library(magrittr)
library(hgu219.db) #library(hgug4112a.db)
library(stringr)
library(annotate)
library(limma)
library(corrplot)
library(GGally)
library(Hmisc)

###############################
###############################

# MAKE CHANGES HERE!!!
save_fitdata <- FALSE
fdr_cutoff <- 0.05
amp_cutoff <- log2(1 + 0.2) 

###############################
###############################

# Read GSE112660 (Wu Hogenesch epidermis PNAS 2018) data and organize data frame
data <- read.csv("../data/GSE112660/GSE112660_timelabelled.csv")

subjects <- sapply(strsplit(colnames(data)[2:80], split="_"), "[", 2)
time <- sapply(strsplit(colnames(data)[2:80], split="_"), "[", 4) %>% as.numeric() %/% 100
colnames(data)[2:80] <- paste(paste0("P",subjects), time, sep="_")
colnames(data)[1] <- "ProbeName"
data <- data[-c(1,2,49389),]
row.names(data) <- data$ProbeName
data %<>% dplyr::select(-ProbeName) 
##data[, c(1:79)] <- sapply(data[, c(1:79)], as.character)
##data[, c(1:79)] <- sapply(data[, c(1:79)], as.numeric)
#data$Symbol <- mapIds(hgu219.db,  as.character(data$ProbeName), keytype = "PROBEID", column = "SYMBOL")

# Associate ProbeIDs to Gene Symbols
genes <- AnnotationDbi::select(hgu219.db, as.character(data %>% tibble::rownames_to_column("ProbeName") %$% ProbeName), 
                               c("PROBEID", "SYMBOL","ENTREZID", "ENSEMBL"))

# Remove duplicates, remove genes lacking IDs...
NoID    <- is.na(genes$ENSEMBL)
#IsExpr  <- rowSums(y$other$gIsWellAboveBG>0) >= 77 #!!! READ
genes <- genes[!NoID,]

data2 <- data %>% mutate(PROBEID = rownames(.)) %>% # Data (expressed, identified) with gene annotation
  inner_join(genes) 

data2 <- avereps(data2, data2[, "ENSEMBL"]) %>% as.data.frame() # Averaging probes mapping to the same gene
data2 %<>% dplyr::select(-SYMBOL, -ENTREZID, -PROBEID)
data2 <- data2 %>% tibble::column_to_rownames("ENSEMBL")
data <- data2
data[, c(1:79)] <- sapply(data[, c(1:79)], as.character)
data[, c(1:79)] <- sapply(data[, c(1:79)], as.numeric)

# Extract sample details from column names
experiment <- data.frame(subject = character(), time = integer()) %>%
  {strcapture("(\\d+)_(\\w+)", colnames(data), ., perl = TRUE)} 

# -----------------------

# Fit linear model (assuming there are no inter-individual differences, ~Fig1)
# Prepare sample details for design matrix
time    <- experiment$time
subject <- factor(experiment$subject)

inphase  <- cos(2*pi*time/24)
outphase <- sin(2*pi*time/24) #a*cos(wt) + cos(wt-phi) == A*cos(wt - phi), where A=sqrt(a**2 + b**2), phi = atan2(b,a)

# design matrix + fits
design <- model.matrix(~ 0 + subject + inphase + outphase) #H0: rhythms are same across subjects

fit  <- limma::lmFit(data, design) 
fit2 <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

if (save_fitdata == TRUE){
  rhy_indices <- which(grepl("phase",colnames(design)))
  results <- limma::topTable(fit2, coef = rhy_indices, number = Inf, sort.by = "none") %>% 
    set_colnames(gsub("\\.","_", colnames(.))) 
  saveRDS(results, file = "../visualize/data/results_GSE112660.rds")
}  

# -----------------------

# Results: rhythmic and expressed genes
clock_genes <- c("PER1","PER2","PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "ARNTL", "ARNTL2", "CLOCK", 
                 "NPAS2","RORA","RORB","RORC", "CSNK1D", "CSNK1E", "DBP")

results <- readRDS("../visualize/data/results_GSE112660.rds")
results %<>% 
  dplyr::mutate(A = sqrt(inphase^2 + outphase^2), #Amp = cos**2 + sin**2 (note the log2 values)
                phase = atan2(outphase, inphase)*12/pi) #atan2 takes two arguments (y,x), atan takes the angle
# phase = 0 means that the respective gene peaks at 12:00 (time of first sampling)

rhy_results <- results %>% filter(adj_P_Val < fdr_cutoff & A > amp_cutoff) %>% 
  tibble::rownames_to_column("ENSEMBL") %>% inner_join( genes %>% dplyr::select(ENSEMBL, SYMBOL) )
rhy_results <- rhy_results[!duplicated(rhy_results$SYMBOL), ]
# -----------------------

# Spearman correlations
df <- data %>% rownames_to_column("ENSEMBL") %>% inner_join( genes %>% dplyr::select(ENSEMBL, SYMBOL) )
df <- df[!duplicated(df$SYMBOL), ]
timeseries_cg <- df[which (df$SYMBOL %in% clock_genes),] %>% 
  dplyr::select(-ENSEMBL) %>%
  gather(time_subj, value, -SYMBOL) %>%
  #dplyr::select(-ProbeName) %>%
  tidyr::separate(time_subj, c("subject","time"), sep = "_", convert = TRUE) 

timeseries_cg_raw <- timeseries_cg
timeseries_cg %<>%
  group_by(subject, SYMBOL) %>%
  mutate(value=value-mean(value)) %>% as.data.frame()

data_corrmat <- timeseries_cg %>% #mutate(gene_tissue = paste(Symbol, tissue, "_")) %>% 
  spread(SYMBOL, value)
data_corrmat_raw <- timeseries_cg_raw %>% #mutate(gene_tissue = paste(Symbol, tissue, "_")) %>% 
  spread(SYMBOL, value)

corrmat <- data_corrmat %>% mutate(subjtime = paste(subject, time, sep="_")) %>%
  dplyr::select(-time, -subject) %>% tibble::column_to_rownames('subjtime') %>% cor(., method="spearman") %>% as.data.frame() %>%
  tibble::rownames_to_column() %>% 
  gather(key, value, -rowname)

corrmat$key <- factor(corrmat$key, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "RORB", "RORC", "CRY1", 
                                            "CRY2", "NR1D1", "NR1D2", "CSNK1D", "CSNK1E", "PER1", "PER2", "PER3", "DBP"))
corrmat$rowname<- factor(corrmat$rowname, levels=c("ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA","RORB","RORC", "CRY1", 
                                                   "CRY2", "NR1D1", "NR1D2", "CSNK1D", "CSNK1E", "PER1", "PER2", "PER3", "DBP"))

# All subjects together, mean substracted correlations
# No significance info
my.lines <- data.frame(x=c(.5,8.5), y=c(8.5,.5), 
                       xend=c(17.5,8.5), yend=c(8.5,17.5))
fig_corr1 <- ggplot(data = corrmat, aes(x=key, y=rowname, fill=value)) + 
  geom_tile() + theme_bw() + xlab("") + ylab("") + labs(fill=expression("Spearman's"~rho)) +
  scale_fill_distiller(palette = "RdBu", limits=c(-1.05,1.05), breaks = c(1, 0.5, 0, -0.5, -1)) + 
  #geom_hline(yintercept = "CRY1"+1) +
  theme(aspect.ratio=1,
        axis.text.x = element_text(angle = 90, vjust = 0.5, face="italic", hjust=1),
        axis.text.y = element_text(face="italic"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        #strip.text = element_text(face="bold"))
        strip.text = element_text(face="bold"),
        legend.title = element_text(vjust = 0.8),
        legend.position="bottom") + 
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F) +
  ggtitle('Spearman correlations time-labelled data GSE112660')

if (!file.exists("../figures/GSE112660_fig1_corr.pdf")){ 
  fig_corr1 %>% ggsave("../figures/GSE112660_fig1_corr.pdf", .) 
}

# With significance info
# DERMIS
data_corrmat <- data_corrmat[,c("subject", "time", "ARNTL", "ARNTL2", "NPAS2", "CLOCK", "RORA", "RORB", "RORC", "CRY1", 
                                        "CRY2", "NR1D1", "NR1D2", "CSNK1D", "CSNK1E", "PER1", "PER2", "PER3", "DBP")] %>%
  dplyr::select(-subject, -time)
col2 = colorRampPalette(c('#053061', '#2166AC', '#4393C3',  '#92C5DE', '#D1E5F0', '#FFFFFF',
                          '#FDDBC7','#F4A582',  '#D6604D', '#B2182B','#67001F'))

signif <- cor.mtest(data_corrmat, conf.level = 0.95)
corr   <- cor(data_corrmat, method='spearman')

pdf("../figures/GSE112660_fig2_corr.pdf")
fig_corr2 <- corrplot(corr, tl.col="black", tl.cex=0.8,tl.srt=45, p.mat = signif$p,  col=col2(200), method='color',
         insig = "blank", sig.level = 0.1, cl.pos="b", cl.cex=0.8)
fig_corr2 <- corrRect(fig_corr2, namesMat = rbind(c('ARNTL', 'CRY2', 'CRY1', 'DBP'), c('CRY2', 'ARNTL', 'DBP', 'CRY1')))
dev.off()

### Subject-by-subject doesn't really make sense (n=4...)

# Scatterplots after mean-centering
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, method = "spearman", ...){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method = method, use = "pair")               # correlation coef
  p <- cor.test(x, y, method = method)$p.val                  # p-value
  txt <- paste0(prefix, format(c(r, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  txt <- format(r, digits = 2)                                # Format r-value
  txt1 <- format(p, digits = 2)                                 # Format p-value
  txt2 <- paste0("r= ", txt, '\n', "p= ", txt1) # Make panel text
  text(0.5, 0.5, txt2, cex = 0.4 + cex.cor * (abs(r)*0.5), ...) # Resize the text by level of correlation
}

pdf(width=11, height=11.5, file="../figures/GSE112660_fig3_corr.pdf")
par(mfrow=c(4,1), cex=0.5)
ggpairs(data_corrmat,  
        upper = list(continuous = wrap(ggally_cor, method='spearman', size = 2, align_percent = 1)), 
        lower = list(continuous = wrap("points", alpha = 0.7, size=0.5, color='orange'),
                     combo = wrap("dot", alpha = 0.7, size=0.5, color='orange')),
        diag = "blank") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 4.5, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 4.5),
        strip.text = element_text(size=7, face="italic"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=10, face='bold')) +
  ggtitle(paste0("Scatterplots clock gene expression GSE112660 (PROCESSED DATA downloaded)\n",
                 "***, **, *, . -> pval < 0.001, 0.01, 0.05, 0.1 respectively"))
dev.off()


