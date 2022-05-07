rm(list = ls())
library(orthopolynom)
library(igraph)
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

times_av = seq(min(exp_index_av),max(exp_index_av),length = 30)
dat_fit_av = power_equation(times_av, dat_par_av)
rownames(dat_fit_av) = rownames(dat)

times_n = seq(min(exp_index_n),max(exp_index_n),length = 30)
dat_fit_n = power_equation(times_n, dat_par_n)
rownames(dat_fit_n) = rownames(dat)

#selection = which(rowSums(dat)>50)
dat_av = dat_av[selection,]
#dat_par_av = dat_par_av[selection,]

dat_n = dat_n[selection,]

################
get_legendre_par <- function(order,exp_index,times) {
  get_interaction <- function(data,col){
    n <- nrow(data)
    clean_data <- data
    gene_list <- list()
    m <- clean_data[,col]
    M <- clean_data[,-col]
    x_matrix <- M
    x_matrix <- as.matrix(x_matrix)
    vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
    x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
    name <- colnames(clean_data)
    ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", 
                           family="gaussian",nfold = 10,alpha = 0)
    best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
    
    fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                         nfold = 10,alpha = 1,
                         penalty.factor = 1/best_ridge_coef,
                         keep = TRUE)
    best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
    
    gene_list_one <- list()
    gene_list_one[[1]] <- name[col]
    gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
    gene_list_one[[3]] <- best_alasso_coef1@x[-1]
    gene_list[[col]] <- gene_list_one
    
    if( identical(gene_list_one[[2]], character(0))==T ){
      gene_list_one[[2]] <- colnames(x_matrix)
      gene_list_one[[3]] <- rep(1,length(colnames(x_matrix)))
    }
    return(gene_list_one)
  }
  
  cluster_mean <- exp_index
  
  module_relationship <- pblapply(1:nrow(cluster_mean),function(c) get_interaction(t(cluster_mean),c))
  #----------------------
  get_effect <- function(pars,effect,times,order){
    if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
    LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
    d_LOP_fit <-  sapply(1:length(pars),function(c)
      pars[c]*polynomial.derivatives(LOP)[[c+1]])
    h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
    #rk4 for legendre with step=h
    LOP_rk4 <- function(x0,y0){
      f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
      k1 <- f(x0,y0) 
      k2 <- f(x0+h/2,y0+h/2*k1)
      k3 <- f(x0+h/2,y0+h/2*k2)
      k4 <- f(x0+h,y0+h*k3)
      y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
      return(y)
    }
    #dy_LOP, the increasment of each step
    dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],0))
    #dy_LOP*y= main or sub effect
    dy_fit <- effect*c(0,dy[1:(length(times)-1)])
    return(cumsum(dy_fit))
  }
  
  ode_optimize <- function(pars,ind,dep,times,data,order){
    ind_pars <- matrix(pars,ncol=order)[1,]
    dep_pars <- matrix(pars,ncol=order)[-1,]
    ind_effect <- get_effect(ind_pars,data[,ind],times,order)
    if ( is.null(nrow(dep_pars)) ) {
      dep_effect <- get_effect(dep_pars,data[,dep],times,order)
      y <- ind_effect+dep_effect+data[,ind][1]
    }else{
      dep_effect <- sapply(1:length(dep), function(c)
        get_effect(dep_pars[c,],data[,dep[c]],times,order))
      y <- ind_effect+rowSums(dep_effect)+data[,ind][1]
    }
    ssr <- sum((data[,ind]-y)^2)
    return(ssr)
  }
  
  get_value <- function(effect,data,times,order){
    #input
    ind <- data[[1]]
    dep <- data[[2]]
    ind_no <- as.numeric(which(colnames(effect)==ind))
    dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
    init_pars <- rep(0.001,(length(ind_no)+length(dep_no))*order)
    result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                    times=times,data=effect,order=order,
                    method = "Nelder-Mead", control=list(maxit=1e4,trace=T))
    par_after <- matrix(result$par,length(ind)+length(dep),order)
    return(par_after)
  }
  
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(orthopolynom)})
  clusterExport(cl, c("get_value","ode_optimize","get_effect","get_interaction",
                      "cluster_mean","module_relationship","order","times_av","times_n")
                ,envir=environment()
  )
  lop_par <- pblapply(1:nrow(cluster_mean),function(c)get_value(t(cluster_mean),
                                                                module_relationship[[c]],
                                                                times,
                                                                order),cl=cl)
  stopCluster(cl)
  return(list(lop_par,module_relationship))
}
#a <- get_legendre_par(1,4)
legendre_order=3
all_lop_par1 <- get_legendre_par(order=legendre_order,dat_fit_n, times_n)
all_lop_par2 <- get_legendre_par(order=legendre_order,dat_fit_av, times_av)
#all_lop_par[[3]] <- get_legendre_par(3,order=legendre_order)

#2.output for result---------------------------------------------------------------------------------------------
get_output <- function(relationship,par,effect,times,order){
  get_effect <- function(pars,effect,times,order){
    if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
    LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
    d_LOP_fit <-  sapply(1:length(pars),function(c)
      pars[c]*polynomial.derivatives(LOP)[[c+1]])
    h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
    #rk4 for legendre with step=h
    LOP_rk4 <- function(x0,y0){
      f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
      k1 <- f(x0,y0) 
      k2 <- f(x0+h/2,y0+h/2*k1)
      k3 <- f(x0+h/2,y0+h/2*k2)
      k4 <- f(x0+h,y0+h*k3)
      y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
      return(y)
    }
    #dy_LOP, the increasment of each step
    dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],0))
    #dy_LOP*y= main or sub effect
    dy_fit <- effect*c(0,dy[1:(length(times)-1)])
    return(cumsum(dy_fit))
  }
  output <- list()
  output[[1]] <- relationship[[1]]  
  output[[2]] <- relationship[[2]]  
  output[[3]] <- par[1,]
  output[[4]] <- par[2:nrow(par),]
  ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
  dep_no <- as.numeric(sapply(1:length(output[[2]]), 
                              function(c) which(colnames(effect)==output[[2]][c])))
  inital_value <- effect[,ind_no][1]/(length(ind_no)+length(dep_no))
  ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order)+inital_value
  if (length(dep_no)==1) {
    dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order)+inital_value
  }else{
    dep_effect <- sapply(1:length(dep_no), function(c)
      get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order))+inital_value
    colnames(dep_effect) <- dep_no
  }
  #------------
  all_effect <- cbind(ind_effect,dep_effect)
  #effect_mean <- all_effect[5,]
  effect_mean <- apply(all_effect,2,mean)
  output[[5]] <- effect_mean
  output[[6]] <- all_effect
  return(output)
}

get_net <- function(exp_index,all_lop_par){
  cluster_mean <- exp_index
  module_relationship <- all_lop_par[[2]]
  net <- pblapply(1:nrow(exp_index),function(c)
    get_output(module_relationship[[c]],all_lop_par[[1]][[c]],
               t(cluster_mean),times=seq(1,ncol(cluster_mean),length=30),order=legendre_order))
  return(net)
}
all_net <- list()
all_net[[1]] <- get_net(dat_fit_n,all_lop_par1)
all_net[[2]] <- get_net(dat_fit_av,all_lop_par2)
get_after <- function(i){
  temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
  temp[,1] <- i[[2]]
  temp[,2] <- i[[1]]
  temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
  
  colnames(temp) <- c('from','to','dep_effect')
  temp <- data.frame(temp)
  temp[,3] <- as.numeric(as.character(temp[,3]))
  return(temp)
}

get_max_effect <- function(k){
  after <- do.call(rbind,lapply(k, get_after))
  max_dep_effect <- max(abs(after$dep_effect))
  
  temp <- aggregate(dep_effect ~ to, data = after, sum)
  all_dep_effect <- max(abs(temp$dep_effect))
  return(c(max_dep_effect,all_dep_effect))
}

max_effect <- t(sapply(1:2,function(c)get_max_effect(all_net[[c]])))
#3.plot---------------------------------------------------------------------------------------------------------

network_plot <- function(k,title){
  
  #extra <- as.numeric(table(df$cluster))
  get_extra <- function(i){
    temp <- i[[5]][1]
    return(temp)
  }
  extra <- sapply(k,get_extra)
  
  after <- do.call(rbind,lapply(k, get_after))
  
  colfunc <- colorRampPalette(c("#619CFF", #ggplot blue
                                "#ffdead", 
                                "#F8766D"))#ggplot red
  #unchange-weight-data
  #links
  #colour_edge
  edge_number <- round(max(max_effect[,1]))
  edge_col <- data.frame(colfunc(2*edge_number+1),seq(-edge_number,edge_number,1))
  
  get_edge_colour <- function(i){
    temp <-  round(i)
    temp2 <- adjustcolor(edge_col[which(edge_col[,2]==temp),1], alpha=0.5)
    return(temp2)
  }
  links <- after
  colnames(links) <- c("from","to","weight")
  links$edge.colour <- pbsapply(links$weight,get_edge_colour) #add colour for links
  
  #nodes
  node_number <- round(max(max_effect[,2]))
  node_col <- data.frame(colfunc(2*node_number+1),seq(-node_number,node_number,1))
  get_vertex_colour <- function(i){
    temp <-  round(i)
    temp2 <- adjustcolor(node_col[which(node_col[,2]==temp),1], alpha=1)
    return(temp2)
  }
  
  nodes <- data.frame(unique(links[,2]),row.names(dat_av),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  nodes <- nodes[order(nodes[,1]),]
  nodes$influence <- aggregate(dep_effect ~ to, data = after, sum)[,2]
  nodes$colour <- pbsapply(nodes$influence,get_vertex_colour) #add colour for links
  
  #normalization
  normalization <- function(x){(x-min(x))/(max(x)-min(x))*1.5+1} 
  
  #final plot
  links[,3] <- normalization(abs(links[,3]))
  nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  net <- graph_from_data_frame( d=links,vertices = nodes,directed = T ) 
  
  #layout
  set.seed(1)
  l <- layout_on_sphere(net)
  
  plot.igraph(net,
               vertex.label=V(net)$name,
               vertex.label.color="black",
               vertex.shape="circle", 
               vertex.label.cex=V(net)$ind_effect,
               vertex.size=V(net)$ind_effect*15,
               edge.curved=0.05,
               edge.color=E(net)$edge.colour,
               edge.frame.color=E(net)$edge.colour,
               edge.width=E(net)$weight,
               vertex.color=V(net)$colour,
               layout=l,
               main=title,
               margin=c(-.05,-.05,-.05,-.05)
  )
  #extra
  #mark.groups=list(c(1,5,21), c(10,23)), 
  #mark.col=c(grDevices::adjustcolor("#F8766D", alpha=0.5),grDevices::adjustcolor("#619CFF", alpha=0.5)),mark.border=NA
}

network_plot(all_net[[1]], title = "")


