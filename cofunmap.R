rm(list = ls())
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)



all_geno<-read.csv("genoH1.csv")
phenoA<-read.csv("H1u.csv")
phenoB<-read.csv("L1u.csv")
phenoB<-phenoB[,-1]
pheno<-phenoA-phenoB


f<- function(i){
  geno=all_geno[i,]
  pheno0<-pheno[which(geno==0),]
  pheno1<-pheno[which(geno==1),]
  pheno9<-pheno[which(geno==9),]
  X <- rbind(pheno0,pheno1)
  
  #Logistic model
  get_miu <- function(par,t){
    t <- 1:8
    y <- par[1]/(1+(par[2]*exp(-par[3]*t)))
    return(y)
  }
  #SAD model
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
  
  get_u = function(par,t){
    get_miu(par[1:3])-get_miu(par[4:6])
  }
  get_u = function(par,t){
    get_SAD1_covmatrix(par[7:8])+get_SAD1_covmatrix(par[9:10])
  }
  
  get_initial_par<- function(pheno,t){
    y <- apply(pheno[,],2,mean)  
    c(max(y),
      max((y[-1]-y[-length(y)])/(t[-1]-t[-length(t)])),
      t[which.max(((y[-1]-y[-length(y)])/(t[-1]-t[-length(t)])))]-y[which.max(((y[-1]-y[-length(y)])/(t[-1]-t[-length(t)])))]/max((y[-1]-y[-length(y)])/(t[-1]-t[-length(t)])))  
  }
  
  get_r2=function(par,t,y){
    sum((y-get_miu(par,t))^2)
  }
  
  get_initial_par2 <- function(pheno,t){
    mcurve =colMeans(pheno)
    init<- optim(par=get_initial_par(phenoA,t),get_r2,t,y=mcurve)
    return(init$par)
  }
  
  L0 <- function(par){
    miu1 <- get_u(par[1:6],t=1:8)
    sigma1 <- get_SAD1_covmatrix(par[7:10],d=8)
    L0 <- -sum(dmvnorm(X,miu1,sigma1,log = T))
    L0
  } 
  #lower=c(11,15,0.5,15,10,1,0.25,0.75),upper = c(15,25,3,25,30,2,0.35,0.85)
  #rm(c,d,dat1,c1,c2)
  inpar <- c(get_initial_par2(phenoA,t),get_initial_par2(phenoB,t),0.5,2)
  NH_0 <- optim(inpar,L0,control=list(maxit=10000))
  NH_0$par
  
  
  L1 <- function(par){
    miu2 <- get_u(par[1:6])
    miu3 <- get_u(par[7:12])
    sigma2 <- get_SAD1_covmatrix(par[13:16],d=8)
    L_0 <- -sum(dmvnorm(pheno0,miu2,sigma2,log = T))
    L_1 <- -sum(dmvnorm(pheno1,miu3,sigma2,log = T))
    LL <- L_0 + L_1
    return(LL)
  }
  
  #rm(c,d,dat1,c1,c2)
  h1_pars <- c(NH_0$par[1:6],NH_0$par[1:6],NH_0$par[7:10])
  NH_1 <- optim(h1_pars,L1,control=list(maxit=10000))
  
  
  LR <- 2*(NH_0$value - NH_1$value)
  
  
  allpar<-c(LR,NH_1$par)
  
  allpar
}