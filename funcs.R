library(magrittr)
library(tibble)
library(tidyr)
library(dplyr) 
library(ggplot2)
library(data.table)
library(tibble)

plot_gene <- function(gene_symbol, tissue){
  # This function plots the expression values as well as the cosine fits from any specific gene across all subjects
  # It takes as arguments the gene that is to be plotted and the tissue of interest (D or E)
  
  probe <- yave$genes[which(yave$genes$Symbol == gene_symbol),]$ProbeName
  
  # extended dataframe of raw expression values
  expr <- yave$E %>% as.data.frame %>% rownames_to_column() %>% filter(rowname == probe) %>%
    select(contains(tissue), -rowname) %>% 
    gather(tissuetimesubj, expr_raw) %>% 
    separate(tissuetimesubj, c("tissue","timesubj"), convert = TRUE, sep = 1) %>% 
    select(-tissue) %>% 
    separate(timesubj, c("time","subject"), sep = "_", convert = TRUE)
  
  # extended dataframe of fitted values
  fitted_vals <- fitted_values %>% as.data.frame %>% rownames_to_column() %>% filter(rowname == probe) %>%
    select(contains(tissue), -rowname) %>% 
    gather(tissuetimesubj, fitted_vals) %>% 
    separate(tissuetimesubj, c("tissue","timesubj"), convert = TRUE, sep = 1) %>% 
    select(-tissue) %>% 
    separate(timesubj, c("time","subject"), sep = "_", convert = TRUE)  
  
  # extended dataframe of fitted_values across all the cosine curve (fitted_curve)
  amp <- results %>% filter(ProbeName == probe) %>% 
    select(contains("A_")) %>% select(contains(tissue)) %>% 
    transpose() %>% rename("amp" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  phi <- results %>% filter(ProbeName == probe) %>% 
    select(contains(paste0("phase", tissue))) %>% select(contains(tissue)) %>%
    transpose() %>% rename("phi" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  baseline <- fit$coefficients[,1:11] %>% as.data.frame() %>% 
    rownames_to_column %>% filter(rowname == probe) %>% select(-rowname) %>% transpose() %>%
    rename("baseline" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  
  if(tissue=="E"){
    diff_baseline_E <- fit$coefficients[,"tissueE"] %>% as.data.frame() %>% 
      rownames_to_column %>% filter(rowname == probe) %$% . %>% as.numeric()
    baseline <- baseline %>% mutate(baseline_E = baseline + diff_baseline_E) %>% #!!! Is this correct?
      select(-baseline) %>% rename("baseline" = "baseline_E")
  } 
  
  t <- seq(min(expr$time), max(expr$time), 0.01)
  fitted_curve <- data.frame(time=t) 
  
  for (i in 1:dim(amp)[1]){
    subj_i <- amp$subject[i]
    amp_i <- amp$amp[i]
    phi_i <- phi$phi[i] #in h
    bli_i <- baseline$baseline[i]
    
    yvals_i <- matrix(0, nrow=length(t)) %>% as.data.frame()
    colnames(yvals_i) <- subj_i
    yvals_i[,1] <- amp_i*cos(2*pi*t/24 - 2*pi*phi_i/24) + bli_i
    
    fitted_curve <- cbind(fitted_curve, yvals_i)
  }
  fitted_curve <- fitted_curve %>% 
    gather(subject, fitted_vals, -time)  
   
  # for an input_check for myself
  fitted_vals_check <- format(round(fitted_vals$fitted_vals, digits=6), nsmall = 6) 
  fitted_curv_check <- format(round(fitted_curve$fitted_vals, digits=6), nsmall = 6) 
  all(fitted_vals_check %in% fitted_curv_check) 
  
  # annotations
  tis <- ifelse(tissue=="D", "D", "E")
  annots <- rhythms_weighted %>% filter(ProbeName == probe & tissue==tis) %>% #filter(tissue == tissue) %>% 
    select(subject, rss, modulus_rhy, modulus_rhywt) 
  #annots[,2:4] <- round(annots[,2:4], 3)
  annots$Label <- paste(paste0("mod=", round(annots$modulus_rhy,3)), paste0("rss=", round(annots$rss,3)), 
                        paste0("mod_w=", round(annots$modulus_rhywt,3)), sep="\n")
  annots$threshold_surpass <- ifelse(annots$modulus_rhywt >= amp_cutoff & annots$modulus_rhy >= amp_cutoff, "a", "b")
  
  # merge all dataframes
  df <- full_join(expr, fitted_vals) %>% full_join(fitted_curve) %>%
    gather(key, value, -time, -subject) %>% full_join(annots)
  
  plot <- ggplot(df) +
    geom_line(data = df %>% filter(key=="fitted_vals"), aes(x=time, y=value), 
              color = ifelse(tissue=="D", "#F8766D","#00AFBB")) + 
    geom_point(data = df %>% filter(key=="expr_raw"), aes(x=time, y=value), color="black", size=2) +
    facet_wrap(~subject) +
    geom_text(data=annots, 
              aes(label = Label, x=Inf, y=Inf, color=threshold_surpass), 
              vjust=1, hjust=1) +
    scale_color_manual(values=c("black", "red")) + 
    theme_bw() + xlab("time (h)") + theme(legend.position = "none") + 
    ylab("expression") + ggtitle(paste0(gene_symbol, " in ", tissue))
  
  return(plot)  
}


circ.mean <- function (x) {
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  circmean
}

plot_allfits <- function(gene_symbol, tissue) {
  probe <- yave$genes[which(yave$genes$Symbol == gene_symbol),]$ProbeName
  
  expr <- yave$E %>% as.data.frame %>% rownames_to_column() %>% filter(rowname == probe) %>%
    select(contains(tissue), -rowname) %>% 
    gather(tissuetimesubj, expr_raw) %>% 
    separate(tissuetimesubj, c("tissue","timesubj"), convert = TRUE, sep = 1) %>% 
    select(-tissue) %>% 
    separate(timesubj, c("time","subject"), sep = "_", convert = TRUE)
  
  # extended dataframe of fitted values
  fitted_vals <- fitted_values %>% as.data.frame %>% rownames_to_column() %>% filter(rowname == probe) %>%
    select(contains(tissue), -rowname) %>% 
    gather(tissuetimesubj, fitted_vals) %>% 
    separate(tissuetimesubj, c("tissue","timesubj"), convert = TRUE, sep = 1) %>% 
    select(-tissue) %>% 
    separate(timesubj, c("time","subject"), sep = "_", convert = TRUE)  
  
  # extended dataframe of fitted_values across all the cosine curve (fitted_curve)
  amp <- results %>% filter(ProbeName == probe) %>% 
    select(contains("A_")) %>% select(contains(tissue)) %>% 
    transpose() %>% rename("amp" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  phi <- results %>% filter(ProbeName == probe) %>% 
    select(contains(paste0("phase", tissue))) %>% select(contains(tissue)) %>%
    transpose() %>% rename("phi" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  baseline <- fit$coefficients[,1:11] %>% as.data.frame() %>% 
    rownames_to_column %>% filter(rowname == probe) %>% select(-rowname) %>% transpose() %>%
    rename("baseline" = "V1") %>% mutate(subject=fitted_vals$subject %>% unique)
  
  if(tissue=="E"){
    diff_baseline_E <- fit$coefficients[,"tissueE"] %>% as.data.frame() %>% 
      rownames_to_column %>% filter(rowname == probe) %$% . %>% as.numeric()
    baseline <- baseline %>% mutate(baseline_E = baseline + diff_baseline_E) %>% #!!! Is this correct?
      select(-baseline) %>% rename("baseline" = "baseline_E")
  } 
  
  t <- seq(min(expr$time), max(expr$time)*2, 0.01)
  
  fitted_curve <- data.frame(time=t) 
  
  for (i in 1:4){#dim(amp)[1]){
    subj_i <- amp$subject[i]
    amp_i <- 1#amp$amp[i]
    phi_i <- phi$phi[i] #in h
    bli_i <- baseline$baseline[i]
    
    yvals_i <- matrix(0, nrow=length(t)) %>% as.data.frame()
    colnames(yvals_i) <- subj_i
    yvals_i[,1] <- amp_i*cos(2*pi*t/24 - 2*pi*phi_i/24) #+ bli_i
    
    fitted_curve <- cbind(fitted_curve, yvals_i)
  }
  fitted_curve <- fitted_curve %>% 
    gather(subject, fitted_vals, -time)  
  
  plot <- ggplot(fitted_curve) + geom_line(aes(x=time, y=fitted_vals, color=subject)) + theme_bw()
  return(plot)
}


plot_gene_rhyU_polar <- function(gene_symbol, tissue){
  if(tissue=="D"){
    X <- XD
    U <- U %>% filter(tissue=="D")
  } else if(tissue=="E"){
    X <- XE
    U <- U %>% filter(tissue=="E")
  }
  genes <- results %>% select(ProbeName, Symbol)
  gene_ofinterest <- genes[which(genes$Symbol == gene_symbol),]
  idx <- which(colnames(XD) == gene_ofinterest$ProbeName)
  
  p <- ggplot(data = data.frame(colX = X[,idx], U = U$U, XcolUh = c(X[,idx]*t(Conj(U$U)))) %>% tibble::rownames_to_column() %>%
           gather(key, value, -rowname) %>% 
           mutate(modulus=Mod(value), argument=Arg(value)),# %>% filter(rowname=="P111_D" | rowname=="P115_D"), 
         aes(x=argument, y=modulus, color=rowname)) + geom_point() + facet_grid(~key) +
    coord_polar(start = pi/2, direction=-1) + theme_bw() +
    geom_segment(aes(x=argument, y=0, xend=argument, yend=modulus)) + 
    scale_x_continuous(limits=c(-pi,pi), labels = math_format(.x * pi, format = function(x) x / pi), 
                       breaks = c(0, pi/2, pi, -pi/2)) +
    ggtitle(paste0("How rhythms of ", gene_symbol, " in ", tissue,
                   " look like in the 11 subjects\n"))
  return(p)
}
