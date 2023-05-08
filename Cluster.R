rm(list = ls())
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)
library("patchwork")

FunMap_par<- allpar
marker_data<-read.csv("geno.csv")
marker_data<-marker_data[,-1]
get_miu <- function(par,t){
  t <- 1:8
  y <- par[1]/(1+(par[2]*exp(-par[3]*t)))
  return(y)
}

get_u = function(par,t){
  get_miu(par[1:3])-get_miu(par[4:6])
}

diff_vg <- c() 

for (a in 1:n) {
  
  AA <- as.numeric(which(marker_data[a,]==1))
  aa <- as.numeric(which(marker_data[a,]==0))
  
  NAA <- length(AA)
  Naa <- length(aa)
  
  p1 <- (NAA*2)/((NAA+Naa)*2) #A????????
  p0 <- (Naa*2)/((NAA+Naa)*2) #a????????
  
  mean_AA <- get_u( as.numeric(FunMap_par[a,2:7]),t=1:8)
  mean_aa <- get_u( as.numeric(FunMap_par[a,8:13]),t=1:8)
  AE <- (mean_AA - mean_aa)/2 
  
  Vg <- 2*p1*p0*(AE^2)  
  diff_vg <- rbind(diff_vg,Vg)
  cat(a,"finished","\n")
  
}

sd<-sqrt(diff_vg)
data<-sd


get_init_par <- function(data,k,legendre_order){
  get_legendre_par <- function(y,legendre_order,x){
    #lm_method
    get_legendre_matrix <- function(x,legendre_order){
      legendre_coef <- legendre.polynomials(n=legendre_order, normalized=F)
      legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
        polynomials=legendre_coef,x=scaleX(x, u=-1, v=1))))
      colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
      return(legendre_matrix[,2:(legendre_order+1)])
    }
    legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
    return(legendre_par)
  }
  #get initial pars based on k-means
  
  init_cluster <- kmeans(data,centers = k,iter.max = 100)
  cuM <- init_cluster$centers 
  
  init_curve_par <- t(sapply(1:k,function(c)get_legendre_par(cuM[c,1:8],legendre_order,1:8)))
  init_SAD_par <- c(1.06,0.25)
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(init_SAD_par,init_curve_par,init_pro)
  names(return_object)<-c("init_SAD_par","init_curve_par","init_pro")
  return(return_object)
}

get_cluster <- function(data,k,input,legendre_order){
  Delta <- 100; iter <- 1; itermax <- 100;
  get_SAD1_covmatrix <- function(par,d){
    phi <- par[1]; gamma <- par[2]; 
    sigma <- array(dim=c(d,d))
    #formula 1, diag element
    diag(sigma) <- sapply(1:d, function(c)(1-phi^(2*c))/(1-phi^2) )
    #formula 2, non-diag element
    sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(d-1),function(c)phi^seq(1:(d-c))*diag(sigma)[c]))
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
    return(gamma^2*sigma)
  }
  
  legendre_fit <- function(par){
    x <- 1:ncol(data)
    fit <- sapply(1:length(par),function(c)par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
    legendre_fit <- as.matrix(as.data.frame(polynomial.values(polynomials=fit,x=scaleX(x, u=-1, v=1))))
    x_interpolation <- rowSums(legendre_fit)
    return(x_interpolation)
  }
  
  mle <- function(par,data,prob){
    par1 <- par[1:2]
    par2 <- matrix(par[-c(1:2)],nrow = k,ncol = (legendre_order+1))
    miu <- t(sapply(1:k, function(c)c(legendre_fit(par2[c,1:(legendre_order+1)]))))
    temp_S <- sapply(1:k,function(c)dmvnorm(data,miu[c,],get_SAD1_covmatrix(par1,d=8))*prob[c])
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  
  cat(paste0("Start biFunClu Calculation ","\n","Cluster_number=",k," Legendre_order=", legendre_order))
  while ( Delta > 1 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <- input[[1]]
      init_curve_par <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(par=old_par,data,prob=pro)
    miu <- t(sapply(1:k, function(c)c(legendre_fit(init_curve_par[c,1:(legendre_order+1)]))))
    
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,miu[c,],get_SAD1_covmatrix(init_SAD_par,d=8))*pro[c] )
    
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:2]
    init_curve_par <- matrix(new_par$par[-c(1:2)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 20000)
      break
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  cat("Finish biFunClu Calculation")
  #plot-----------
  cluster <- apply(omega,1,which.max)
  clustered<- data.frame(row.names(data),data[,1:8],cluster)
  clustered_data <- clustered
  get_plot <- function(clustered_data){
    colnames(clustered_data) <- c("marker",1:8,"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") + xlab("Time")+ylab("generic_effect")
    return(p)
  }
  
  p1 <- get_plot(clustered)
  clustered<- clustered[,-1]
  return_object <- list(init_SAD_par,init_curve_par,pro,LL_mem,BIC,clustered_ck,p1) #cluster????
  cat("Finish biFunClu Calculation")
  names(return_object)<-c("SAD_par", "curve_par", "pro", "LL", 
                          "BIC", "clustered_ck","plot1")
  return(return_object)
  
}


set.seed(500) 
biFunClu_intial_pars<- get_init_par(data=data,k=4,legendre_order=4)
biFunClu_results<- get_cluster(data=data,k=4,input=biFunClu_intial_pars,legendre_order=4)



