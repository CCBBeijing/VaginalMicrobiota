rm(list = ls())
library(nlme)
library(drc)
library(aomisc)
library(pbapply)
library(ggplot2)
library(patchwork)
library(gridExtra)
get_p1_base <- function(){
  dat = read.table("otu_table.p.absolute.xls.txt",row.names = 1, header = T)[1:240]
  
  dat_av = dat[,1:80]
  dat_n = dat[,81:240]
  
  exp_index_save = log2(colSums(dat)+1)
  exp_index_av = log2(colSums(dat_av)+1)
  exp_index_n = log2(colSums(dat_n)+1)
  
  dat_av = log2(dat_av+1)
  dat_n = log2(dat_n+1)
  dat = log2(dat+1)
  #exp_index_save = colSums(dat)
  
  power_fit <- function(x,y){
    x <- as.numeric(x)
    y <- as.numeric(y)
    model <- try(nls(y~a*x^b,start = list(a =0.01, b = 1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =0.3, b = 0.1),
                       control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =3, b = 0.1),
                       control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =0.1, b = 0.1),
                       control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a = -0.1, b = 0.1),
                       control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a = 0.8, b = -0.1),
                       control = nls.control(maxiter = 100000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =-0.1, b = -0.1),
                       control = nls.control(maxiter = 1000000,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      result <- c(NA,NA)
    }
    else{
      result <- coef(model)}
    return(result)
  } 
  
  power_equation <- function(x, dat_par){ t(sapply(1:nrow(dat_par),function(c) dat_par[c,1]*x^dat_par[c,2] ) )}
  
  selection = which(rowSums(dat)>50)
  dat = dat[selection,]
  dat_par_av = t(sapply(1:nrow(dat),function(c) power_fit(exp_index_av, dat_av[c,])))
  dat_par_n = t(sapply(1:nrow(dat),function(c) power_fit(exp_index_n, dat_n[c,])))
  
  
  #selection = which(rowSums(dat)>50)
  dat_av = dat_av[selection,]
  #dat_par_av = dat_par_av[selection,]
  
  dat_n = dat_n[selection,]
  #dat_par_n = dat_par_n[selection,]
  
  #dat = dat[selection,]
  
  library(reshape2)
  library(grid)
  
  colnames(dat_av) = exp_index_av
  df_av = melt(as.matrix(dat_av))
  df_av$group = "AV"
  
  colnames(dat_n) = exp_index_n
  df_n = melt(as.matrix(dat_n))
  df_n$group = "N"
  
  df = rbind(df_av,df_n)
  
  ###########3
  times = seq(min(exp_index_save),max(exp_index_save),length = 100)
  dat_fit_av = power_equation(times, dat_par_av)
  rownames(dat_fit_av) = rownames(dat)
  colnames(dat_fit_av) = times
  
  df.fit_av = melt(as.matrix(dat_fit_av))
  df.fit_av$group = "AV"
  
  dat_fit_n = power_equation(times, dat_par_n)
  rownames(dat_fit_n) = rownames(dat)
  colnames(dat_fit_n) = times
  
  df.fit_n = melt(as.matrix(dat_fit_n))
  df.fit_n$group = "N"
  
  df.fit = rbind(df.fit_av,df.fit_n)
  
  xlabels = expression(2^14,2^15,2^16,2^17,2^18)
  ylabels = expression(0,2^5,2^10,2^15,2^20)
  
  p <- ggplot() + 
    geom_point(df, mapping = aes(x = Var2, y = value,colour = group),  
               size = 1.25, show.legend = F, alpha = 0.25) +
    geom_line(df.fit, mapping = aes(x = Var2, y = value,colour = group), size = 1.15, show.legend = F)  + 
    scale_color_manual(values = c("AV" = "#FF7BA9", "N" = "#65C18C")) + 
    facet_wrap(~Var1, nrow = 3) + 
    geom_text(df, mapping = aes(label = Var1), 
              x = 16.3, y = 18.5, check_overlap = TRUE, size = 3) + 
    xlab("Habitat Index") + ylab("Abundance of Individual Phyla") + theme(axis.title=element_text(size=18)) +
    scale_x_continuous(labels = xlabels) + scale_y_continuous(labels = ylabels, limits = c(0,19.5))
  
  p1 <- p + 
    theme_bw() + #theme(panel.grid =element_blank()) +
    theme(panel.spacing = unit(0.0, "lines")) +
    theme(plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(strip.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(margin=margin(0,0,0,0,"pt")))+
    theme(axis.text.y = element_text(hjust = 0))
  
  #ggsave("Fig1.png",p1,width = 10, height = 4)
  return(p1)
}

get_p2_base <- function(){
  
  dat =read.delim("C:/Users/ddaa/Desktop/wq_microbe/otu_table.g.absolute.xls.txt", header=T)
  dat = aggregate(dat[,c(-1,-242)], list(gene=dat[,1]), FUN = sum)
  rownames(dat) = dat[,1]
  dat = dat[,-1]
  
  dat_av = dat[,1:80]
  dat_n = dat[,81:240]
  
  exp_index_save = log2(colSums(dat)+1)
  exp_index_av = log2(colSums(dat_av)+1)
  exp_index_n = log2(colSums(dat_n)+1)
  
  dat_av = log2(dat_av+1)
  dat_n = log2(dat_n+1)
  dat = log2(dat+1)
  #exp_index_save = colSums(dat)
  
  power_fit <- function(x,y){
    #tmp <- data.frame(x=as.numeric(x),y=as.numeric(y))  
    #model <- try(drm(y ~ x, fct = DRC.powerCurve(),
    #                 data = tmp))
    x <- as.numeric(x)
    y <- as.numeric(y)
    model <- try(nls(y~a*x^b,start = list(a =0.01, b = 1),
                     control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =0.3, b = 0.1),
                       control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =3, b = 0.1),
                       control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =0.1, b = 0.1),
                       control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a = -0.1, b = 0.1),
                       control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a = 0.8, b = -0.1),
                       control = nls.control(maxiter = 1e4,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      model <- try(nls(y~a*x^b,start = list(a =-0.1, b = -0.1),
                       control = nls.control(maxiter = 1e40,minFactor = 1e-200)))
    }
    if ('try-error' %in% class(model)) {
      tmp <- data.frame(x=as.numeric(x),y=as.numeric(y))  
      model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                       data = tmp))
    }
    if ('try-error' %in% class(model)) {
      result <- c(NA,NA)
    }
    else{
      result <- coef(model)}
    return(result)
  } 
  
  power_equation <- function(x, dat_par){ t(sapply(1:nrow(dat_par),function(c) dat_par[c,1]*x^dat_par[c,2] ) )}
  
  
  selection = c("Lactobacillus", "Streptococcus", "Aerococcus", "Gardnerella",
                  "Atopobium", "Prevotella", "Pseudomonas","Ureaplasma",
                  "Prevotella","Staphylococcus", "Escherichia" ,"Others")
  select_name = c("Lactobacillus", "Streptococcus", "Aerococcus", "Gardnerella",
                "Atopobium", "Prevotella", "Pseudomonas","Ureaplasma",
                "Prevotella","Staphylococcus", "Escherichia" ,"Others")
  dat_av = dat_av[selection,]
  dat_par_av = t(pbsapply(1:nrow(dat_av),function(c) power_fit(exp_index_av, dat_av[c,])))
  
  #tmp = rowSums(dat)[order(rowSums(dat))]
  #tmp
  dat_n = dat_n[selection,]
  dat_par_n = t(pbsapply(1:nrow(dat_n),function(c) power_fit(exp_index_n, dat_n[c,])))
  
  dat = dat[selection,]
  
  
  library(reshape2)
  library(grid)
  dat_av2 = dat_av[match(select_name,rownames(dat)),]
  dat_n2 = dat_n[match(select_name,rownames(dat)),]
  
  colnames(dat_av2) = exp_index_av
  df_av = melt(as.matrix(dat_av2))
  df_av$group = "AV"
  
  colnames(dat_n2) = exp_index_n
  df_n = melt(as.matrix(dat_n2))
  df_n$group = "N"
  
  df = rbind(df_av,df_n)
  
  ###########3
  times = seq(min(exp_index_save),max(exp_index_save),length = 100)
  dat_fit_av = power_equation(times, dat_par_av[match(select_name,rownames(dat)),])
  rownames(dat_fit_av) = rownames(dat_av2)
  colnames(dat_fit_av) = times
  
  df.fit_av = melt(as.matrix(dat_fit_av))
  df.fit_av$group = "AV"
  
  dat_fit_n = power_equation(times, dat_par_n[match(select_name,rownames(dat)),])
  rownames(dat_fit_n) = rownames(dat_n2)
  colnames(dat_fit_n) = times
  
  df.fit_n = melt(as.matrix(dat_fit_n))
  df.fit_n$group = "N"
  
  df.fit = rbind(df.fit_av,df.fit_n)
  
  xlabels = expression(2^14,2^15,2^16,2^17,2^18)
  ylabels = expression(0,2^5,2^10,2^15,2^20)
  
  p <- ggplot() + 
    geom_point(df, mapping = aes(x = Var2, y = value,colour = group),  
               size = 1.25, show.legend = F, alpha = 0.25) +
    geom_line(df.fit, mapping = aes(x = Var2, y = value,colour = group), size = 1.15, show.legend = F)  + 
    scale_color_manual(values = c("AV" = "#FF7BA9", "N" = "#65C18C")) + 
    facet_wrap(~Var1, nrow = 2) + 
    geom_text(df, mapping = aes(label = Var1), 
              x = 16.3, y = 18.5, check_overlap = TRUE, size = 3) + 
    xlab("Habitat Index") + ylab("Abundance of Individual Genera") + theme(axis.title=element_text(size=18)) +
    scale_x_continuous(labels = xlabels) + scale_y_continuous(labels = ylabels, limits = c(0,19.5))
  
  p1 <- p + 
    theme_bw() + #theme(panel.grid =element_blank()) +
    theme(panel.spacing = unit(0.0, "lines")) +
    theme(plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(strip.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(margin=margin(0,0,0,0,"pt"))) +
    theme(axis.text.y = element_text(hjust = 0))
  
  #ggsave("Fig2.png",p1,width = 10, height = 4)
  return(p1)
}

p1 = get_p1_base()
p2 = get_p2_base()


pp = p1/p2 + plot_annotation(tag_levels = 'A') + plot_layout(heights =  c(3, 2))

ggsave(filename = "Fig2.pdf",pp, width = 13, height = 9)

