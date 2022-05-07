rm(list = ls())

library(ggplot2)

load("microbe_g.Rdata")


df = data.frame(cbind(dat_fit_n, dat_fit_av))
times = c(times_n,times_av)

get_SAD1_covmatrix <- function(par,n){
  phi <- par[1]; gamma <- par[2];
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  return(gamma^2*sigma)
}

get_biSAD1 <- function(par,n1,n2){
  sig1 <- get_SAD1_covmatrix(par[1:2],n1)
  sig2 <- get_SAD1_covmatrix(par[3:4],n2)
  sig12 <- array(0, dim=c(n1,n2))
  sig21 <- array(0, dim=c(n2,n1))
  sigma1 <- cbind(sig1,sig12)
  sigma2 <- cbind(sig21,sig2)
  sigma <- rbind(sigma1,sigma2)
  return(sigma)
}

#EM
Q.function <- function(par, omega, X, k, times){
  n = dim(X)[1]; d = dim(X)[2]
  par.cov <- par[1:4]
  par.mu <-  matrix(par[-c(1:4)], nrow = k, ncol = 4)
  mu <- cbind(power_equation(times[1:30],par.mu[1:k,1:2]), power_equation(times[31:60], par.mu[1:k,3:4]))
  cov <- get_biSAD1(par.cov, n1=30,n2=30)
  mvn.log <- sapply(1:k, function(c) dmvnorm(X, mu[c,], cov, log = TRUE))
  #mvn.log[mvn.log==-Inf] = log(1e-250)
  Q.LL <- -sum(omega * mvn.log)
  return(Q.LL)
}


mle <- function(par, prob, X, k){
  n = dim(X)[1]; d = dim(X)[2]
  par.cov <- par[1:4]
  par.mu <- matrix(par[-c(1:4)], nrow = k, ncol = 4)
  mu <- cbind(power_equation(times[1:30],par.mu[1:k,1:2]), power_equation(times[31:60], par.mu[1:k,3:4]))
  cov <- get_biSAD1(par.cov, n1=30,n2=30)
  mvn <- sapply(1:k, function(c) dmvnorm(X, mu[c,], cov)*prob[c] )
  #for (i in 1:n) {
  #if (any(mvn[i,]==0)) {
  #mvn[i,which(mvn[i,]==0)] = 1e-250 * prob[which(mvn[i,]==0)]
  #}
  #}
  mvn[mvn==0] = 1e-250
  LL <- sum(-log(rowSums(mvn)))
  return(LL)
}

power_fit <- function(y,times){
  tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
  model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                   data = tmp))
  x <- as.numeric(times)
  y <- as.numeric(y)
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

get.par.int <- function(X, k, times){
  n = dim(X)[1]; d = dim(X)[2]; 
  #use k-means to acquire initial parameters
  init.cluster <- kmeans(X,centers = k,iter.max = 1e3)
  prob <- table(init.cluster$cluster)/n 
  cov.par.int <- c(0.6,1.2,0.8,1.2)
  
  mu.par.int.n <- t(sapply(1:k,function(c) power_fit(init.cluster$centers[c,1:30],times[1:30])))
  mu.par.int.av <- t(sapply(1:k,function(c) power_fit(init.cluster$centers[c,31:60],times[31:60])))
  mu.par.int = cbind(mu.par.int.n, mu.par.int.av)
  colnames(mu.par.int) = c("par.a_n","par.b_n","par.a_av","par.b_av")
  #check
  #plot(times, init.cluster$centers[2,], type = 'l')
  #lines(times, legendre_fit(mu.par.int[4,], times),col=2 )
  
  return_obj <- list(initial_cov_params = cov.par.int,
                     initial_mu_params = mu.par.int,
                     initial_probibality = prob)
  return(return_obj)
}
#initial.pars = get.par.int(df, 6, c(exp_index_n,exp_index_av))
#X = df;k=5;prob = a$probibality;par.int = c(a$cov_par,a$mu_par)
fun.clu <- function(X, k, inv.cov = NULL, eplison = 1, iter.max = 1e1){
  #attribute
  n = dim(X)[1]; d = dim(X)[2]; times = c(exp_index_n,exp_index_av); iter = 0;
  # initial pars
  initial.pars = get.par.int(X, k, times)
  
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob <-  initial.pars$initial_probibality
  
  while( abs(eplison) > 1e-1 && iter <= iter.max ){
    #E step
    par.cov <- par.int[1:4]
    par.mu <- matrix(par.int[-c(1:4)],nrow = k,ncol = 4 )
    mu <- cbind(power_equation(times[1:30], par.mu[1:k,1:2]), power_equation(times[31:60], par.mu[1:k,3:4]))
    cov <- get_biSAD1(par.cov, n1=30, n2 = 30)
    mvn <- sapply(1:k, function(c) dmvnorm(X, mu[c,], cov)* prob[c])

    
    mvn[mvn==0] = 1e-250
    omega <- mvn/rowSums(mvn)
    
    LL_mem <-  mle(par.int, prob, X, k)
    
    #M step
    prob <- colSums(omega)/sum(omega)
    mle.maximization <- try(optim(par = par.int, mle, 
                                prob = prob, 
                                X = X,
                                k = k,
                                method = "Nelder-Mead",
                                #method = "L-BFGS-B",
                                #lower = -1e3,
                                #upper = 1e3,
                                control = list(trace = TRUE,maxit = 1000
                                )))
    
    par.hat <- mle.maximization$par  
    par.int = par.hat
    
    LL_next <-  mle(par.hat, prob, X, k)
    eplison <-  LL_next - LL_mem
    LL_mem <- LL_next
    
    iter = iter + 1
    
    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL_next, "\n")
    
  }
  BIC <- 2*(LL_next) + log(n)*length(par.hat)
  BIC2 =  LL_next + log(n*d)*length(par.hat)
  X.clustered <- data.frame(X, apply(omega,1,which.max))
  
  return_obj <- list(cluster_number = k,
                     Log_likelihodd = LL_next,
                     BIC_value = BIC,
                     BIC2 = BIC2 ,
                     cov_par = par.cov,
                     mu_par = par.mu,
                     probibality = prob,
                     cluster = X.clustered)
  return(return_obj)
}

##################################################33
#load("funclu_k_6.Rdata")

a = fun.clu(df,k = 3, iter.max = 10)


table(a$cluster$apply.omega..1..which.max.)


result_all = replicate(10, fun.clu(df,k = 15, iter.max = 100))

result=a
fun.clu.plot2 <- function(result){
  k = result$cluster_number
  times = c(times_n, times_av)
  par.mu = result$mu_par
  mu.fit = cbind(power_equation(times[1:30], par.mu[1:k,1:2]), power_equation(times[31:60], par.mu[1:k,3:4]))
  mu.fit_n = data.frame(mu.fit[,1:30])
  mu.fit_n$cluster = 1:k
  mu.fit_n = melt(mu.fit_n,id=c("cluster"))
  mu.fit_n$time = rep(times_n, each = k)
  mu.fit_n$color = paste0("N_",mu.fit_n$cluster)
  mu.fit_n$type = "N"
  
  mu.fit_av = data.frame(mu.fit[,31:60])
  mu.fit_av$cluster = 1:k
  mu.fit_av = melt(mu.fit_av,id=c("cluster"))
  mu.fit_av$time = rep(times_av, each = k)
  mu.fit_av$color = paste0("AV_",mu.fit_av$cluster)
  mu.fit_av$type = "AV"
  
  plot.df_n = result$cluster[,c(1:30,61)]
  #plot.df_n$apply.omega..1..which.max. = paste0("N_",plot.df$apply.omega..1..which.max.)
  plot.df_n$name = rownames(plot.df_n)
  plot.df_n = melt(plot.df_n,id=c("apply.omega..1..which.max.", "name"))
  plot.df_n = plot.df_n[order(plot.df_n$name),]
  plot.df_n$time = rep(times_n,106)
  colnames(plot.df_n)[1] = c("cluster")
  plot.df_n$type = "N"
  plot.df_n$color = paste0("N_",plot.df_n$cluster)
  
  plot.df_av = result$cluster[,c(31:61)]
  #plot.df_av$apply.omega..1..which.max. = paste0("AV_",plot.df$apply.omega..1..which.max.)
  plot.df_av$name = rownames(plot.df_av)
  plot.df_av = melt(plot.df_av,id=c("apply.omega..1..which.max.", "name"))
  plot.df_av = plot.df_av[order(plot.df_av$name),]
  plot.df_av$time = rep(times_av,106)
  colnames(plot.df_av)[1] = c("cluster")
  plot.df_av$type = "AV"
  plot.df_av$color = paste0("AV_",plot.df_av$cluster)
  
  name.df = data.frame(label = paste0("M_",1:k), x = 16.25, y = 20, cluster = 1:k)
  
  p = ggplot() + geom_point(plot.df_n, mapping = aes(x = time, y = value, colour = type), alpha = 0.35) +
    geom_point(plot.df_av, mapping = aes(x = time, y = value, colour = type), alpha = 0.35)+
    geom_line(mu.fit_n, mapping = aes(x = time, y = value, colour = type),size=1.25)+
    geom_line(mu.fit_av, mapping = aes(x = time, y = value, colour = type),size=1.25)+
    facet_wrap(~cluster) + scale_color_manual(values = c("AV" = "#FF7BA9", "N" = "#65C18C")) + 
    xlab("Habitat Index") + ylab("Module Expression") + theme(axis.title=element_text(size=18)) +
    scale_x_continuous(labels = math_format(expr = 2^.x)) + 
    scale_y_continuous(labels = math_format(expr = 2^.x)) +
    theme(legend.position="none") + 
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
fun.clu.plot2(a)











