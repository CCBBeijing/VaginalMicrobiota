rm(list = ls())
library(ggplot2)
library(mvtnorm)
library(nlme)
library(drc)
library(aomisc)
library(reshape2)
library(scales)
load("microbe_g.Rdata")
df = data.frame(cbind(dat_n, dat_av))

get_SAD1_covmatrix <- function(par,n){
  phi <- par[1]; gamma <- par[2];
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  #sigma[sigma<0.25] = 0
  return(gamma^2*sigma)
}

#EM
mle <- function(par, prob, X, k){
  df_n = X[,1:160]; df_av = X[,161:240];
  n = dim(X)[1]; d1 = dim(df_n)[2]; d2 = dim(df_av)[2];
  times = c(exp_index_n,exp_index_av)
  par.cov <- par[1:4]
  par.mu <- matrix(par[-c(1:4)], nrow = k, ncol = 4)
  cov_n = get_SAD1_covmatrix(par.cov[1:2], n = d1)
  cov_av = get_SAD1_covmatrix(par.cov[3:4], n = d2)
  mu_n <- power_equation(times[1:160], par.mu[1:k,1:2])
  mu_av <- power_equation(times[161:240], par.mu[1:k,3:4])
  mvn_n <- sapply(1:k, function(c) dmvnorm(df_n, mu_n[c,], cov_n))
  mvn_av <- sapply(1:k, function(c) dmvnorm(df_av, mu_av[c,], cov_av))
  mvn = sapply(1:k,function(c) mvn_n[,c]* prob[c]*mvn_av[,c])
  mvn[mvn==0] = 1e-250
  LL <- sum(-log(rowSums(mvn)))
  return(LL)
}

power_fit <- function(y,times){
  tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
  model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                   data = tmp))
  #x <- as.numeric(times)
  #y <- as.numeric(y)
  #model <- try(nls(y~a*x^b,start = list(a =5, b = 0.5),
  #                control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.3, b = 0.1),
                     control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =3, b = 0.1),
                     control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.1, b = 0.1),
                     control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = -0.1, b = 0.1),
                     control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = 0.8, b = -0.1),
                     control = nls.control(maxiter = 10000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =-0.1, b = -0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    result <- c(NA,NA)
  }
  else{
    result <- coef(model)}
  return(result)
}

power_equation <- function(x, dat_par){ t(sapply(1:nrow(dat_par),function(c) dat_par[c,1]*x^dat_par[c,2] ) )}

get.par.int <- function(X, k, times){
  n = dim(X)[1]; d = dim(X)[2]; 
  
  cov.int = c(0.7,1.3,0.8,1.6)
  #use k-means to acquire initial parameters
  init.cluster <- kmeans(X,centers = k)
  prob <- table(init.cluster$cluster)/n 
  
  mu.par.int.n <- t(sapply(1:k,function(c) power_fit(init.cluster$centers[c,1:160],times[1:160])))
  mu.par.int.av <- t(sapply(1:k,function(c) power_fit(init.cluster$centers[c,161:240],times[161:240])))
  mu.par.int = cbind(mu.par.int.n, mu.par.int.av)
  colnames(mu.par.int) = c("par.a_n","par.b_n","par.a_av","par.b_av")
  
  
  return_obj <- list(initial_cov_params = cov.int,
                     initial_mu_params = mu.par.int,
                     initial_probibality = prob, 
                     model = init.cluster
  )
  return(return_obj)
}


fun.clu <- function(X, k, inv.cov = NULL, eplison = 1, iter.max = 1e1){
  #attribute
  df_n = X[,1:160]; df_av = X[,161:240];times = c(exp_index_n,exp_index_av);
  
  n = dim(X)[1]; d1 = dim(df_n)[2]; d2 = dim(df_av)[2]; iter = 0;
  # initial pars
  initial.pars = get.par.int(X, k, times)
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob <-  initial.pars$initial_probibality
  #prob = rep(1/k,k)
  
  while( abs(eplison) > 1e-1 && iter <= iter.max ){
    #E step
    par.mu <- matrix(par.int[-c(1:4)],nrow = k,ncol = 4 )
    par.cov = par.int[1:4]
    cov_n = get_SAD1_covmatrix(par.cov[1:2], n = d1)
    cov_av = get_SAD1_covmatrix(par.cov[3:4], n = d2)
    mu_n <- power_equation(times[1:160], par.mu[1:k,1:2])
    mu_av <- power_equation(times[161:240], par.mu[1:k,3:4])
    
    mvn_n <- sapply(1:k, function(c) dmvnorm(df_n, mu_n[c,], cov_n))
    mvn_av <- sapply(1:k, function(c) dmvnorm(df_av, mu_av[c,], cov_av))
    #mean(mvn_av)
    prob = prob 
    mvn = sapply(1:k,function(c) mvn_n[,c]*prob[c]*mvn_av[,c])
    mvn[mvn==0]=1e-250
    
    omega <- mvn/rowSums(mvn)
    table(apply(omega,1,which.max))
    LL_mem <-  mle(par.int, prob, X, k)
    
    #M stepx
    prob <- colSums(omega)/sum(omega)

    m.maximization <- try(optim(par = par.int, mle, 
                                prob = prob, 
                                X = X,
                                k = k,
                                method = "Nelder-Mead",
                                #method = "BFGS",
                                #lower = -1e3,
                                #upper = 1e3,
                                control = list(trace = TRUE,maxit = 1e3
                                )))
    
    par.hat <- m.maximization$par
    par.int = par.hat
    
    LL_next <-  mle(par.hat, prob, X, k)
    eplison <-  LL_next - LL_mem
    LL_mem <- LL_next
    
    iter = iter + 1
    
    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL_next, "\n")
    
  }
  BIC <- 2*(LL_next) + log(n)*length(par.hat)
  
  #tmp = which(sapply(1:n, function(c) all(omega[c,]== 1/k))==T)
  X.clustered <- data.frame(X, apply(omega,1,which.max))
  #X.clustered$apply.omega..1..which.max.[tmp] = initial.pars$cluster[tmp]
  table(X.clustered$apply.omega..1..which.max.)
  return_obj <- list(cluster_number = k,
                     Log_likelihodd = LL_next,
                     BIC_value = BIC,
                     cov_par = par.cov,
                     mu_par = par.mu,
                     probibality = prob,
                     omega = omega,
                     cluster = X.clustered)
  return(return_obj)
}

a = fun.clu(df,k = 3, iter.max = 10)

fun.clu.plot <- function(result){
  k = result$cluster_number
  times = c(exp_index_n, exp_index_av)
  par.mu = result$mu_par
  mu.fit = cbind(power_equation(times[1:160], par.mu[1:k,1:2]), power_equation(times[161:240], par.mu[1:k,3:4]))
  mu.fit_n = data.frame(mu.fit[,1:160])
  mu.fit_n$cluster = 1:k
  mu.fit_n = melt(mu.fit_n,id=c("cluster"))
  mu.fit_n$time = rep(exp_index_n, each = k)
  mu.fit_n$color = paste0("N_",mu.fit_n$cluster)
  mu.fit_n$type = "N"
  
  mu.fit_av = data.frame(mu.fit[,161:240])
  mu.fit_av$cluster = 1:k
  mu.fit_av = melt(mu.fit_av,id=c("cluster"))
  mu.fit_av$time = rep(exp_index_av, each = k)
  mu.fit_av$color = paste0("AV_",mu.fit_av$cluster)
  mu.fit_av$type = "AV"
  
  plot.df_n = result$cluster[,c(1:160,241)]
  #plot.df_n$apply.omega..1..which.max. = paste0("N_",plot.df$apply.omega..1..which.max.)
  plot.df_n$name = rownames(plot.df_n)
  plot.df_n = melt(plot.df_n,id=c("apply.omega..1..which.max.", "name"))
  plot.df_n = plot.df_n[order(plot.df_n$name),]
  plot.df_n$time = rep(exp_index_n,106)
  colnames(plot.df_n)[1] = c("cluster")
  plot.df_n$type = "N"
  plot.df_n$color = paste0("N_",plot.df_n$cluster)
  
  plot.df_av = result$cluster[,c(161:241)]
  #plot.df_av$apply.omega..1..which.max. = paste0("AV_",plot.df$apply.omega..1..which.max.)
  plot.df_av$name = rownames(plot.df_av)
  plot.df_av = melt(plot.df_av,id=c("apply.omega..1..which.max.", "name"))
  plot.df_av = plot.df_av[order(plot.df_av$name),]
  plot.df_av$time = rep(exp_index_av,106)
  colnames(plot.df_av)[1] = c("cluster")
  plot.df_av$type = "AV"
  plot.df_av$color = paste0("AV_",plot.df_av$cluster)
  
  name.df = data.frame(label = paste0("M_",1:k), x = 16.25, y = 20, cluster = 1:k)
  
  p = ggplot() + geom_point(plot.df_n, mapping = aes(x = time, y = value, colour = type), alpha = 0.35) +
    geom_point(plot.df_av, mapping = aes(x = time, y = value, colour = type), alpha = 0.35)+
    geom_line(mu.fit_n, mapping = aes(x = time, y = value, colour = type),size=1.25)+
    geom_line(mu.fit_av, mapping = aes(x = time, y = value, colour = type),size=1.25)+
    facet_wrap(~cluster
               #,scales = "free_y"
    ) + 
    scale_color_manual(values = c("AV" = "#FF7BA9", "N" = "#65C18C")) + 
    xlab("Habitat Index") + ylab("Module Expression") + theme(axis.title=element_text(size=18)) +
    scale_x_continuous(labels = math_format(expr = 2^.x)) + 
    scale_y_continuous(labels = math_format(expr = 2^.x)) +
    theme(legend.position="none")+ 
    geom_text(name.df, mapping = aes(x = x, y = y,label = label), check_overlap = TRUE, size = 5) 
  
  p1 <- p + 
    theme_bw() + #theme(panel.grid =element_blank()) +
    theme(panel.spacing = unit(0.0, "lines")) +
    theme(plot.margin = unit(c(1,1,1,1), "lines")) +
    theme(strip.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(margin=margin(0,0,0,0,"pt")))
  
  #ggsave("Fig_cluster.png",p1,width = 10, height = 4)
  return(p1)
}

fun.clu.plot(a)

