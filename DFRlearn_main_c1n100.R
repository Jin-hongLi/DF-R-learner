######################################################################
#########  Simulation code For DF R-learner              #############
######################################################################
rm(list = ls())
setwd(getwd())
######################################################################
#########  Load packages                                 #############
######################################################################
library(Matrix)
library(stringr)
library(mclust) # Adjusted Rand Index
library(glmnet)
library(ggplot2)
######################################################################
#########  parameters                                    #############
######################################################################
c=1 #Level of heterogeneity
n=100 #Sample size in each study
p=4 #Number of covariates (including intercept)
K=10 #Number of studies
nk=rep(n,K) #List of all sample sizes
nk.test=rep(100,K) #List of all sample sizes for testing set
f=5 #Folds of cross fitting
random=F #Random trial indicator
discrete=T #Discrete grouping indicator
N=sum(nk) #Total sample size
rho=1 #Parameter for ADMM 
alpha=1 #Parameter for fusion penalty, not for glmnet
alpha.glmnet=1 #Parameter for glmnet
r=0.5 #Parameter for adaptive weight
s=0.5 #Parameter for group-wise weight
zeta=1 # parameter for EBIC
threshold = 1e-8 #Convergence threshold
n_digit=3L #integer indicating the number of decimal places to round coefficients
maxiter=2500 #Maximum number of iterations

vk<-rep(1,K) #List of standard deviation for random error in each study

lambda <- 10^(seq(-6,0,length.out=100)) 
date<-"250701"
######################################################################
#########             Essential Functions                #############
######################################################################
# generate coefficients for m(), e() and tau()
gencoef<-function(c=1,K=10,discrete=F,p=4,seed=NA){
  if (!is.na(seed)){
    set.seed(seed)
  }
  beta_tau_null<-matrix(nrow = p,ncol = K)
  beta_tau_null[,]<-t(c(1.5,2,-3,1))
  beta_m_null<-matrix(nrow = p,ncol = K)
  beta_m_null[,]<-t(c(1,0.5,-1.5,2))
  beta_e_null<-matrix(nrow = p,ncol = K)
  beta_e_null[,]<-t(c(0.2,0.1,-0.3,0.2))
  g1k<-rbinom(n=K, size=1, prob=0.5)#grouping for beta1
  g2k<-rbinom(n=K, size=1, prob=0.5)#grouping for beta2
  if (discrete){
    beta_tau_null[2,]<-beta_tau_null[2,]+c*g1k
    beta_tau_null[3,]<-beta_tau_null[3,]+c*g2k
  }else {
    beta_tau_null[2,]<-beta_tau_null[2,]+c*g1k+runif(K,-0.1,0.1)
    beta_tau_null[3,]<-beta_tau_null[3,]+c*g2k+runif(K,-0.1,0.1)
  }
  beta_m_null[2,]<-beta_m_null[2,]+g1k #c=1 for m(x)
  beta_m_null[3,]<-beta_m_null[3,]+g2k
  beta_e_null[2,]<-beta_e_null[2,]+0.2*g1k #c=0.2 for e(x)
  beta_e_null[3,]<-beta_e_null[3,]+0.2*g2k
  beta_tau_null<-t(beta_tau_null)
  beta_m_null<-t(beta_m_null)
  beta_e_null<-t(beta_e_null)
  result <- list("beta_m_null"=beta_m_null, "beta_e_null"=beta_e_null,
                 "beta_tau_null"=beta_tau_null) 
  return(result)
}

# generate simulation data with defined coefficients
generatedata<-function(nk,vk,p,K=10,random=F,beta_m,beta_e,beta_tau){
  K=K
  for (k in 1:K){
    x<-c()
    prob<-c()
    ecolumn<-c()
    n<-nk[k]
    v<-vk[k]
    x_temp<-rep(1,n) #intercept, x_0
    j=0
    while (j<(p-1)) {
      x_temp<-cbind(x_temp,rnorm(n))
      j=j+1
    } 
    if (random){
      prob<-rep(0.5,n)
    } else {
      odds<-x_temp[,1:p]%*%beta_e[k,]
      prob<-exp(odds)/(1+exp(odds))
    }
    for (e in 1:n) {
      ecolumn[e]<-sample(0:1,1,replace=T,prob = c(1-prob[e],prob[e]))
    }
    x_temp<-cbind(x_temp,ecolumn)
    x<-as.matrix(x_temp)
    colnames(x)<-c(paste(paste("x",k,sep = ""),0:(p-1),sep = "_"),paste("trt",k,sep = "_"))
    #generate y by Robinson's transformation
    y<-x[,1:p]%*%beta_m[k,]+diag(as.list(x[,p+1]-prob))%*%(x[,1:p]%*%beta_tau[k,])+rnorm(nrow(x),sd=v)
    data<-cbind(y,x)
    colnames(data)[1]<-"y"
    data<-data.frame(data)
    data.sim[[k]]<-data
  }
  return(data.sim)
}

# IPD version of DF-R-learner (FR-IPD)
FRIPD <- function(X,Y,rho,r,s,lambda,alpha=0,zeta=1, nk, vk, n_digit=3L,
                  K,maxiter=2500,threshold=1e-8){
  p=ncol(X)/K
  N<-nrow(X)
  #initialize beta
  beta <- solve(crossprod(X)) %*% crossprod(X,Y)
  rownames(beta)<-colnames(X)
  #generate the D matrix
  if (K==1){
    #error message
    stop("K should be at least 2")
    D <- diag(p)
  } else{
    #generate D matrix for multiple studies
    D <- matrix(0, nrow=(K)*p, ncol=(K)*p) 
    for (j in 1:p){
      beta.j<-beta[c((p)*(c(0:(K-1)))+j),1]
      orderU<-order(beta.j)
      orderV<-order(abs(beta.j))
      for (i in 1:(K-1)){
        D[1+(j-1)*K,(orderV[1]-1)*p+j] <- 1 #the first line indicates smallest absolute value
        D[i+1+(j-1)*K,(orderU[i]-1)*p+j] <- (-1)
        D[i+1+(j-1)*K,(orderU[i+1]-1)*p+j] <- 1
      }
    }
  }
  #initialize theta
  theta<-as.vector(D%*%beta) #for calculating W
  #generate W matrix for weights 
  W<-rep(0,(K)*p)
  for (j in 1:p){
    W[1+(j-1)*K]<-alpha*abs(theta[1+(j-1)*K])^(-r)
    thetarange<-sum(theta[(2+(j-1)*K):(j*K)])
    for (i in 2:K){
      W[i+(j-1)*K]<-abs(thetarange)^(-s)*abs(theta[i+(j-1)*K])^(-r)
    }
  }
  W<-replace(W,W==Inf,1e8)
  Vinv<-c()
  for(i in 1:K){
    Vinv<-c(Vinv,rep(vk[i],nk[i])^-1)
  }
  Vinv<-diag(Vinv) #Lambda^(-1)
  
  theta <- matrix(0, ncol=1, nrow=nrow(D)) #initial value
  kappa <- matrix(0, ncol=1, nrow=nrow(D)) #initial value
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  converge<-0
  p1<-solve(t(X)%*%Vinv%*%X+N*rho*t(D)%*%D) #avoid repeat computation
  p2<-t(X)%*%Vinv%*%Y #avoid repeat computation
  # ADMM updates
  for (t in 2:maxiter){
    #Update beta
    beta.new<-p1 %*% (p2 + N*rho*t(D)%*%(theta-kappa))
    #Update theta (soft threshold)
    theta.new<-sign(D%*%beta.new+kappa) * 
      pmax(abs(D%*%beta.new+kappa) - lambda*W/rho, 0)
    #Update gamma
    kappa.new=kappa+D%*%beta.new-theta.new
    if ((sqrt(sum((beta.new-beta)^2))/ncol(X)) <= threshold){
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      converge<-1
      break
    }
    else{
      beta=beta.new
      theta=theta.new
      kappa=kappa.new
    }
  }
  
  #final result of beta
  beta.hat=beta.new  
  beta.hat.round<-round(beta.hat,n_digit)
  #Elements for BIC
  residual <- as.vector(Y - X %*% beta.hat)
  stp<-1
  n2ll<-c()
  df<-c()
  combin<-0
  for (i in 1:K){
    n2ll[i]<-nk[i]*log(sum(residual[stp:(stp+nk[i]-1)]^2)/(nk[i]))
    stp<-stp+nk[i]
  }
  for (j in 1:p){
    beta.j<-beta.hat.round[c((p)*(c(0:(K-1)))+j)]
    df[j]<-length(unique(beta.j))
    combin<-combin+choose(K,df[j])
  }
  #EBIC
  bic<-sum(n2ll)+sum(df)*log(N)+2*zeta*log(combin)
  #AIC
  aic<-sum(n2ll)+2*sum(df)
  
  #reformat beta 
  beta.hat.matrix<-data.frame(matrix(beta.hat,ncol=p,byrow = T))
  rownames(beta.hat.matrix)<-paste("study",1:K,sep = "_")
  colnames(beta.hat.matrix)<-c("Intercept",paste0("x",1:(p-1),sep=""))
  beta.hat.matrix<-data.matrix(beta.hat.matrix,rownames.force=T)
  result <- list("beta.hat.matrix"=beta.hat.matrix, "beta.hat"=beta.hat,
                 "conv_message"=message, "iter"=t,"BIC"=bic,"AIC"=aic,
                 "converge"=converge)   
  return(result)
}

# CD version of DF-R-learner (DFR-CD)
DFRCD <- function(sigma_n_inv,b,rho,r,s,lambda, n_digit=3L, nk,alpha=0,
                 zeta=1,
                 K,maxiter=2500,threshold=1e-8){
  p=ncol(sigma_n_inv)/K
  N=sum(nk)
  #initialize beta matrix
  beta <- matrix(b, ncol=1, nrow=ncol(sigma_n_inv)) 
  rownames(beta)<-colnames(sigma_n_inv)
  #generate the D matrix
  if (K==1){
    #error message
    stop("K should be at least 2")
    D<-diag(p)
  } else{
    #generate D matrix for multiple studies 
    D <- matrix(0, nrow=(K)*p, ncol=(K)*p) 
    for (j in 1:p){
      beta.j<-beta[c((p)*(c(0:(K-1)))+j),1]
      orderU<-order(beta.j)
      orderV<-order(abs(beta.j))
      for (i in 1:(K-1)){
        D[1+(j-1)*K,(orderV[1]-1)*p+j]<-1
        D[i+1+(j-1)*K,(orderU[i]-1)*p+j]<-(-1)
        D[i+1+(j-1)*K,(orderU[i+1]-1)*p+j]<-1
      }
    }
  }
  #initialize theta
  theta<-as.vector(D%*%beta) # for calculating W
  #generate W matrix for weights
  W<-rep(0,(K)*p)
  for (j in 1:p){
    W[1+(j-1)*K]<-alpha*abs(theta[1+(j-1)*K])^(-r)
    thetran<-sum(theta[(2+(j-1)*K):(j*K)])
    for (i in 2:K){
      W[i+(j-1)*K]<-abs(thetran)^(-s)*abs(theta[i+(j-1)*K])^(-r)
    }
  }
  W<-replace(W,W==Inf,1e8)
  theta <- matrix(0, ncol=1, nrow=nrow(D)) #initial value
  kappa <- matrix(0, ncol=1, nrow=nrow(D)) #initial value
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  converge<-0
  p1<-solve(sigma_n_inv+N*rho*crossprod(D))
  p2<-sigma_n_inv%*%b
  # ADMM updates
  for (t in 2:maxiter){
    #Update beta
    beta.new<-p1 %*% (p2 + N*rho*t(D)%*%(theta-kappa))
    #Update theta
    theta.new<-sign(D%*%beta.new+kappa) * 
      pmax(abs(D%*%beta.new+kappa) - lambda*W/rho, 0)
    #Update kappa
    kappa.new=kappa+D%*%beta.new-theta.new
    # Check convergence
    # Define Euclidean (l2) norm of a vector
    if ((sqrt(sum((beta.new-beta)^2))/ncol(sigma_n_inv)) <= threshold ){
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      converge<-1
      break
    }
    else{
      beta=beta.new
      theta=theta.new
      kappa=kappa.new
    }
  }
  
  #final result of beta
  beta.hat=beta.new   
  beta.hat.round<-round(beta.hat,n_digit)
  #elements for BIC
  df<-c()
  combin<-0
  for (j in 1:p){
    beta.j<-beta.hat.round[c((p)*(c(0:(K-1)))+j)]
    df[j]<-length(unique(beta.j))
    combin<-combin+choose(K,df[j])
  }
  #EBIC
  bic<-t(b-beta.hat)%*%sigma_n_inv%*%(b-beta.hat)+sum(df)*log(N)+2*zeta*log(combin)
  #AIC
  aic<-t(b-beta.hat)%*%sigma_n_inv%*%(b-beta.hat)+2*sum(df)
  #reformat beta 
  beta.hat.matrix<-data.frame(matrix(beta.hat,ncol=p,byrow = T))
  rownames(beta.hat.matrix)<-paste0("study",1:K,sep="")
  colnames(beta.hat.matrix)<-c("Intercept",paste0("x",1:(p-1),sep=""))
  beta.hat.matrix<-data.matrix(beta.hat.matrix,rownames.force=T)
  result <- list("beta.hat.matrix"=beta.hat.matrix, "beta.hat"=beta.hat,
                 "conv_message"=message, "iter"=t,"BIC"=bic,"AIC"=aic,
                 "converge"=converge)  
  return(result)
}

# Multi Study R learner 
# (https://github.com/cathyshyr/multi-study-r-learner/blob/main/Utils.R) 
ms_rlearner <- function(X, Y, alpha, K, p, p_mat, newX){
  ms_rlearner_mod <- glmnet::cv.glmnet(x = as.matrix(X),
                                   y = as.matrix(Y),
                                   alpha = alpha,
                                   lambda = NULL,
                                   standardize = TRUE,
                                   intercept = FALSE,
                                   standardize.response = TRUE)
  
  # Remove the first entry that corresponds to the intercept (which is empty)
  tau_beta <- as.vector(t(coef(ms_rlearner_mod, s = "lambda.min")[-1]))
  # p + 1 coefficients per study (p covariates and 1 intercept)
  tau_beta_k <- split(tau_beta, ceiling(seq_along(tau_beta)/(p + 1)))
  
  ms_rlearner_pred_list <- vector("list", length = K)
  for(k in 1:K){
    ms_rlearner_pred_list[[k]] <- newX %*% as.matrix(tau_beta_k[[k]]) * p_mat[, k]
  }
  ms_rlearner_pred_mat <- do.call("cbind", ms_rlearner_pred_list)
  ms_rlearner_pred <- apply(ms_rlearner_pred_mat, 1, sum)
  return(ms_rlearner_pred)
}

# Generate labels by rounding values and assigning integers
generatelabel <- function(x) {
  x <- round(x, 3)
  labels <- match(x, unique(x))
  return(labels)
}

# Extract beta indices for covariate (1-based indexing)
get_covariate_indices <- function(covariate_index, K, p) {
  return(seq(from = covariate_index, to = K * p, by = p))
}

# calculate comparing metrics
eva_perform<-function(true, est,K,p){
  bias_beta <- mean(abs(true - est))  # MAB
  mse_beta <- mean((true - est)^2)    # MSE
  ari_by_cov=c()
  for (pp in 1:p) {
    idx <- get_covariate_indices(pp, K, p)
    if (discrete){
      beta_true_pp <- true[idx]
    } else {
      beta_true_pp<-round(true[idx]*2,0)/2 # for non-discrete grouping
    }
    
    beta_est_pp <- est[idx]
    
    label_true <- generatelabel(beta_true_pp)
    label_est <- generatelabel(beta_est_pp)
    ari_by_cov=c(ari_by_cov,adjustedRandIndex(label_true, label_est))
  }
  ari_beta=mean(ari_by_cov)
  return(list(ari_beta=ari_beta,
              bias_beta=bias_beta,
              mse_beta=mse_beta))
}
######################################################################
#########         loop over seed                         #############
######################################################################
# number 1000 can be changed
seedlist<-(1:1000)
result=c() #store results for all methods
result.data.all=c()
beta.est.all=c()
for (seed in seedlist){
 
  tryCatch(expr={
    set.seed(seed)
    print(seed)
    ######################################################################
    #########  store true and fitted coefficients            #############
    ######################################################################
    beta.est<-c()
    ######################################################################
    #########  generate coefficients                         #############
    ######################################################################
    simcoef<-gencoef(c=c,K=K,discrete=discrete,p=p,seed=NA) #seed is pre-specified
    beta_m=simcoef$beta_m_null
    beta_e=simcoef$beta_e_null
    beta_tau=simcoef$beta_tau_null
    beta.est<-matrix(t(beta_tau),nrow=1)#first line is true value
    ######################################################################
    #########  generate data                                 #############
    ######################################################################
    data.sim<-c()
    data.sim<-generatedata(nk=nk,vk=vk,p=p,K=K,
                           random=random,beta_e=beta_e,beta_m=beta_m,beta_tau=beta_tau)
    data.test<-generatedata(nk=nk.test,vk=vk,p=p,K=K,
                            random=random,beta_e=beta_e,beta_m=beta_m,beta_tau=beta_tau)
    ######################################################################
    #########          prepare data                          #############
    ######################################################################
    #For testing data
    tau.test.true=c()
    m_hat_test_all<-c()
    e_hat_test_all<-c()
    for (k in 1:K){
      data.temp=as.matrix(data.test[[k]])
      tau.test.true=c(tau.test.true,data.temp[,2:(p+1)]%*%(as.vector(beta_tau[k,])))
      }
    
    # Oracle Setting
    {start_time <- Sys.time()
      X.or<-c()
      Y.or<-c()
      X_tilde_true<-c() #X_tilde
      Y_tilde_true<-c() #Y_tilde
      Colnames<-c()
      v.or<-c()
      for (k in 1:K){
        n<-nk[k]
        data<-as.matrix(data.sim[[k]])
        y<-data[,1]
        x<-data[,-1]
        if (random){
          prob<-rep(0.5,n)
        } else {
          odds<-x[,1:(p)]%*%beta_e[k,]
          prob<-exp(odds)/(1+exp(odds))
        }
        y_tilde_true<-y-x[,1:(p)]%*%beta_m[k,]
        x_tilde_true<-diag(as.list(x[,p+1]-prob))%*%x[,1:(p)]
        Y_tilde_true[[k]]<-y_tilde_true
        X_tilde_true[[k]]<-x_tilde_true
        Colnames<-c(Colnames,colnames(X_tilde_true[[k]])) #store all colnames
        Y.or<-c(Y.or,y_tilde_true)
        fit.or <- glm(y_tilde_true~x_tilde_true-1,family="gaussian")
        v.or<-c(v.or,sum(fit.or$residuals^2)/n)
      }
      X.or<-as.matrix(bdiag(X_tilde_true)) #block diagonal matrix to store all X's
      colnames(X.or)<-Colnames
      end_time <- Sys.time()
      prepare.time.or<-end_time - start_time}
    
    # IPD Setting
    start_time <- Sys.time()
    X_tilde_est<-c() #List of X_tilde
    Y_tilde_est<-c() #List of Y_tilde
    X.ipd<-c() #store all X_tilde in a matrix
    Y.ipd<-c() #store all Y_tilde in a vector
    m_hat_all<-c()
    e_hat_all<-c()
    Colnames<-c()
    v.ipd<-c()
    for (k in 1:K){
      n<-nk[k]
      data.temp<-as.matrix(data.sim[[k]])
      y<-data.temp[,1]
      x<-data.temp[,-1]
      # fold ID for cross-validation
      foldid = sample(rep(seq(f), length = length(y)))
      
      m_hat<-c()
      e_hat<-c()
      #######cross fitting m
      m_fit = glmnet::cv.glmnet(x[,-c(1,p+1)], y,
                                foldid = foldid,
                                intercept=T,
                                keep = TRUE,
                                lambda = NULL,
                                alpha = alpha.glmnet,
                                penalty.factor = rep(1,p-1))
      m_lambda_min = m_fit$lambda[which.min(m_fit$cvm[!is.na(colSums(m_fit$fit.preval))])]
      m_hat = m_fit$fit.preval[,!is.na(colSums(m_fit$fit.preval))][, m_fit$lambda[!is.na(colSums(m_fit$fit.preval))] == m_lambda_min]
      m_hat_all<-c(m_hat_all,m_hat)
      #######cross fitting e
      if (random){
        e_hat=rep(0.5,length(y))
      } else {
        e_fit = glmnet::cv.glmnet(x[,-c(1,p+1)], x[,(p+1)],
                                  foldid = foldid,
                                  intercept=T,
                                  lambda = NULL,
                                  keep = TRUE,
                                  alpha = alpha,
                                  penalty.factor = rep(1,p-1))
        
        e_lambda_min = e_fit$lambda[which.min(e_fit$cvm[!is.na(colSums(e_fit$fit.preval))])]
        e_hat = e_fit$fit.preval[,!is.na(colSums(e_fit$fit.preval))][, e_fit$lambda[!is.na(colSums(e_fit$fit.preval))] == e_lambda_min]
      }
      e_hat_all<-c(e_hat_all,e_hat)
      
      y_tilde_est<-y-m_hat
      x_tilde_est<-diag(x[,p+1]-e_hat)%*%x[,1:(p)] #exclude trt
      Y.ipd<-c(Y.ipd,y_tilde_est)
      X_tilde_est[[k]]<-x_tilde_est
      Y_tilde_est[[k]]<-y_tilde_est
      Colnames<-c(Colnames,colnames(x_tilde_est))
      fit.ipd <- glm(y_tilde_est~x_tilde_est-1,family="gaussian")
      v.ipd<-c(v.ipd,sum(fit.ipd$residuals^2)/n)
    }
    X.ipd<-as.matrix(bdiag(X_tilde_est))
    colnames(X.ipd)<-Colnames
    end_time <- Sys.time()
    prepare.time.ipd<-(end_time - start_time)
    
    #CD
    start_time <- Sys.time()
    sigma_n_inv_est<-c() #sigma_n inverse matrix
    beta.cd<-c()
    Colnames<-c()
    for (k in 1:K){
      n<-nk[k]
      y_tilde_est<-Y_tilde_est[[k]]
      x_tilde_est<-X_tilde_est[[k]]
      x_data_est<-cbind(y_tilde_est,x_tilde_est)
      colnames(x_data_est)[1]<-"y"
      x_data_est<-data.frame(x_data_est)
      fit.cd <- glm(y~.-1,data=x_data_est,family="gaussian")
      beta.cd<-c(beta.cd,fit.cd$coefficients)
      v<-solve(solve(t(x_tilde_est)%*%x_tilde_est)*sum(fit.cd$residuals^2)/n)#sigma_n_inverse
      Colnames<-c(Colnames,colnames(x_tilde_est))
      sigma_n_inv_est[[k]]<-v
    }
    sigma_n_inv_est.cd<-as.matrix(bdiag(sigma_n_inv_est))
    colnames(sigma_n_inv_est.cd)<-Colnames
    end_time <- Sys.time()
    prepare.time.cd<-(end_time - start_time)+prepare.time.ipd
    
    ######################################################################
    #########    estimate local R-learner LM (LOC)        #############
    ######################################################################
    beta.LOC=beta.cd
    beta.est <- rbind(beta.est,beta.LOC)
    colnames(beta.est)<-Colnames
    
    tau.LOC=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.LOC=c(tau.LOC,as.matrix(data.temp[,2:(p+1)])%*%beta.LOC[((k-1)*p+1):(k*p)])
    }
    
     performance=eva_perform(beta.est[1,], beta.LOC,K,p)
    
     result.LOC=data.frame("seed"=seed,"method"="LOC",
                           "ARI"=performance$ari_beta,
                           "MSE"=performance$mse_beta,
                           "MSE.TAU"=mean((tau.LOC - tau.test.true)^2),
                           "BIAS"=performance$bias_beta)
    result=rbind(result,result.LOC)
   
    ######################################################################
    #########    estimate oracle                             #############
    ######################################################################
    betamat <- matrix(rep(NA, length(lambda)*K*p), ncol=length(lambda))
    bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    bicmat[1,]  <- lambda
    convergedmat[1,]  <- lambda
    itermat[1,]  <- lambda
    start_time <- Sys.time()
    for (i in 1:length(lambda)) {
      dfmodel <- FRIPD(X=X.or,Y=Y.or,rho=rho,r=r,s=s,lambda=lambda[i], K=K, 
                       threshold=threshold, nk=nk,vk=v.or, alpha=alpha, 
                       zeta=zeta, n_digit=n_digit,maxiter = maxiter)
      betamat[,i]  <- dfmodel$beta.hat #store estimated beta
      bicmat[2,i]  <- dfmodel$BIC #store EBIC
      convergedmat[2,i]  <- dfmodel$converge #store convergence
      itermat[2,i]  <- dfmodel$iter #store iteration
    }
    end_time <- Sys.time()
    model.time.or<-(end_time - start_time)/length(lambda)
    rownames(betamat)<-Colnames
    index.lambda.best<-which.min(bicmat[2,])
    beta.or<-betamat[,index.lambda.best]
    beta.est<-rbind(beta.est,beta.or)
    lambda.best.or<-lambda[index.lambda.best]
    converged.or<-convergedmat[2,index.lambda.best]
    iter.or<-itermat[2,index.lambda.best]

    tau.or=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.or=c(tau.or,as.matrix(data.temp[,2:(p+1)])%*%beta.or[((k-1)*p+1):(k*p)])
    }
    
    performance=eva_perform(beta.est[1,], beta.or,K,p)
    
    result.OR=data.frame("seed"=seed,"method"="OR",
                             "ARI"=performance$ari_beta,
                             "MSE"=performance$mse_beta,
                             "MSE.TAU"=mean((tau.or - tau.test.true)^2),
                             "BIAS"=performance$bias_beta)
    result=rbind(result,result.OR)
    ######################################################################
    #########    estimate FR-IPD                                #############
    ######################################################################
    betamat <- matrix(rep(NA, length(lambda)*K*p), ncol=length(lambda))
    bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    bicmat[1,]  <- lambda
    convergedmat[1,]  <- lambda
    itermat[1,]  <- lambda
    start_time <- Sys.time()
    for (i in 1:length(lambda)) {
      dfmodel <- FRIPD(X=X.ipd,Y=Y.ipd,rho=rho,r=r,s=s,lambda=lambda[i], K=K, 
                       threshold=threshold, nk=nk,vk=v.ipd, alpha=alpha,
                       zeta=zeta, n_digit=n_digit,maxiter = maxiter)
      betamat[,i]  <- dfmodel$beta.hat #store estimated beta
      bicmat[2,i]  <- dfmodel$BIC #store EBIC
      convergedmat[2,i]  <- dfmodel$converge #store convergence
      itermat[2,i]  <- dfmodel$iter #store iteration
    }
    end_time <- Sys.time()
    model.time.IPD<-(end_time - start_time)/length(lambda)
    rownames(betamat)<-Colnames
    index.lambda.best<-which.min(bicmat[2,])
    beta.IPD<-betamat[,index.lambda.best]
    beta.est<-rbind(beta.est,beta.IPD)
    lambda.best.IPD<-lambda[index.lambda.best]
    converged.IPD<-convergedmat[2,index.lambda.best]
    iter.IPD<-itermat[2,index.lambda.best]
    
    tau.IPD=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.IPD=c(tau.IPD,as.matrix(data.temp[,2:(p+1)])%*%beta.IPD[((k-1)*p+1):(k*p)])
    }
    
    performance=eva_perform(beta.est[1,], beta.IPD,K,p)
    
    result.FRIPD=data.frame("seed"=seed,"method"="FR-IPD",
                         "ARI"=performance$ari_beta,
                         "MSE"=performance$mse_beta,
                         "MSE.TAU"=mean((tau.IPD - tau.test.true)^2),
                         "BIAS"=performance$bias_beta)
    result=rbind(result,result.FRIPD)
    
    ######################################################################
    #########    estimate  CD                                #############
    ######################################################################
    betamat <- matrix(rep(NA, length(lambda)*K*p), ncol=length(lambda))
    bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
    bicmat[1,]  <- lambda
    convergedmat[1,]  <- lambda
    itermat[1,]  <- lambda
    start_time <- Sys.time()
    for (i in 1:length(lambda)) {
      dfmodel <- DFRCD(sigma_n_inv=sigma_n_inv_est.cd,b=beta.cd,rho=rho,r=r,s=s,nk=nk,lambda=lambda[i], K=K, 
                      threshold=threshold, alpha=alpha,
                      zeta=zeta, n_digit=n_digit,maxiter = maxiter)
      betamat[,i]  <- dfmodel$beta.hat #store estimated beta
      bicmat[2,i]  <- dfmodel$BIC #store EBIC
      convergedmat[2,i]  <- dfmodel$converge #store convergence
      itermat[2,i]  <- dfmodel$iter #store iteration
    }
    end_time <- Sys.time()
    model.time.CD<-(end_time - start_time)/length(lambda)
    rownames(betamat)<-Colnames
    index.lambda.best<-which.min(bicmat[2,])
    beta.CD<-betamat[,index.lambda.best]
    beta.est<-rbind(beta.est,beta.CD)
    lambda.best.CD<-lambda[index.lambda.best]
    converged.CD<-convergedmat[2,index.lambda.best]
    iter.CD<-itermat[2,index.lambda.best]

    tau.CD=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.CD=c(tau.CD,as.matrix(data.temp[,2:(p+1)])%*%beta.CD[((k-1)*p+1):(k*p)])
    }
    performance=eva_perform(beta.est[1,], beta.CD,K,p)
    
    result.DFRCD=data.frame("seed"=seed,"method"="DFR-CD",
                          "ARI"=performance$ari_beta,
                          "MSE"=performance$mse_beta,
                          "MSE.TAU"=mean((tau.CD - tau.test.true)^2),
                          "BIAS"=performance$bias_beta)
    result=rbind(result,result.DFRCD)
    ######################################################################
    #########    estimate K-Means                            #############
    ######################################################################
    betalist<-c()
    for (j in 1:p){
      beta.j<-beta.LOC[p*((1:K)-1)+j] 
      #calculate center=1 first
      betamat<-rep(kmeans(beta.j, centers = 1, nstart=10)$centers,K)
      # calculate center=2 to 9, because center=10==LOC
      for (i in 2:(K-1)){
        km<-kmeans(beta.j, centers = i, nstart=10)
        km.cluster<-km$cluster
        km.mean<-km$centers
        returnmean<-function(x){return(km.mean[x])}#given cluster, give mean
        km.beta<-sapply(km.cluster,returnmean)
        betamat<-rbind(betamat,km.beta)
      }
      #here, include center=10
      betamat<-rbind(betamat,beta.j)
      betalist[[j]]<-betamat
    }
    #matrix of all combination of number of centers (1-10)
    combnmat<-expand.grid(rep(list(1:K), p))
    calculatebic<-function(x){
      comb<-combnmat[x,]
      beta<-c()
      for (j in 1:p){
        beta<-c(beta,betalist[[j]][unlist(comb[j]),])
      }
      beta<-matrix(t(matrix(beta,nrow=K)),ncol =1)
      beta.round<-round(beta,n_digit)
      df<-c()
      combin<-0
      for (j in 1:p){
        beta.j<-beta.round[c((p)*(c(0:(K-1)))+j)]
        df[j]<-length(unique(beta.j))
        combin<-combin+choose(K,df[j])
      }
      bic<-t(beta-beta.LOC)%*%sigma_n_inv_est.cd%*%(beta-beta.LOC)+sum(df)*log(N)+2*zeta*log(combin)
      return(bic)
    }
    biclist<-sapply(1:nrow(combnmat),calculatebic)
    comb<-combnmat[which.min(biclist),]
    beta.KM<-c()
    for (j in 1:p){
      beta.KM<-c(beta.KM,betalist[[j]][unlist(comb[j]),])
    }
    beta.KM<-t(matrix(t(matrix(beta.KM,nrow=K)),ncol = 1))
    beta.est<-rbind(beta.est,beta.KM)
    
    tau.KM=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.KM=c(tau.KM,as.matrix(data.temp[,2:(p+1)])%*%beta.KM[((k-1)*p+1):(k*p)])
    }
    performance=eva_perform(beta.est[1,], beta.KM,K,p)
    
    result.KM=data.frame("seed"=seed,"method"="KM",
                         "ARI"=performance$ari_beta,
                         "MSE"=performance$mse_beta,
                         "MSE.TAU"=mean((tau.KM - tau.test.true)^2),
                         "BIAS"=performance$bias_beta)
    result=rbind(result,result.KM)
    
    ######################################################################
    #########    estimate MA                                 #############
    ######################################################################
    beta.MA<-c()
    for (j in 1:p){
      beta.j<-beta.LOC[c((p)*(c(0:(K-1)))+j)]
      beta.MA<-c(beta.MA,rep(weighted.mean(beta.j, nk), K))
    }
    beta.MA<-t(matrix(t(matrix(beta.MA,nrow=K)),ncol = 1))
    beta.est<-rbind(beta.est,beta.MA)
    
    tau.MA=c()
    for (k in 1:K){
      data.temp=data.test[[k]]
      tau.MA=c(tau.MA,as.matrix(data.temp[,2:(p+1)])%*%beta.MA[((k-1)*p+1):(k*p)])
    }
    performance=eva_perform(beta.est[1,], beta.MA,K,p)
    
    result.MA=data.frame("seed"=seed,"method"="MA",
                         "ARI"=performance$ari_beta,
                         "MSE"=performance$mse_beta,
                         "MSE.TAU"=mean((tau.MA - tau.test.true)^2),
                         "BIAS"=performance$bias_beta)
    result=rbind(result,result.MA)
    
    ######################################################################
    #########    Multi-Study R-learner                       #############
    ######################################################################
    data.ms=c()
    data.test.ms=c()
    colnames.ms=sapply(0:(p-1), function(x) paste("x", x, sep = ""))
    for (k in 1:K){
      data.temp=data.sim[[k]]
      data.temp$S=k
      colnames(data.temp)[2:(p+1)]<-colnames.ms
      colnames(data.temp)[p+2]<-"trt"
      data.ms[[k]]=data.temp
      #test set
      data.test.temp=data.test[[k]]
      colnames(data.test.temp)[2:(p+1)]<-colnames.ms
      colnames(data.test.temp)[p+2]<-"trt"
      data.test.ms[[k]]=data.test.temp
    }
    data.ms.all=as.data.frame(do.call("rbind", data.ms))
    x.ms.all=as.matrix(data.ms.all[, 2:(p+1)])
    #test set
    data.test.ms.all=as.matrix(do.call("rbind", data.test.ms))
    x.test.ms.all=as.matrix(data.test.ms.all[, 2:(p+1)])
    # Estimate ascertainment probability 
    p_mod <- glmnet::cv.glmnet(x = x.ms.all,
                               y = as.matrix(data.ms.all$S),
                               intercept=F,
                               alpha = alpha.glmnet,
                               lambda = NULL,
                               standardize = TRUE,
                               keep = TRUE,
                               family = "multinomial")
    p_hat_lambda_min = p_mod$lambda[which.min(p_mod$cvm[!is.na(colSums(p_mod$fit.preval))])]
    p_pred_mat <- as.data.frame(predict(p_mod, newx = as.matrix(data.ms.all[, 2:(p+1)]), s = p_hat_lambda_min, type = "response"))
    p_pred_mat_test <- as.data.frame(predict(p_mod, newx = as.matrix(data.test.ms.all[, 2:(p+1)]), s = p_hat_lambda_min, type = "response"))
    
    X_res_mat_merge <- lapply(1:K, function(k){
      sweep(x.ms.all, 1, (data.ms.all$trt - e_hat_all) * p_pred_mat[, k], "*")
    })
    X_res_mat_merge_all <- as.data.frame(do.call("cbind", X_res_mat_merge))
    names(X_res_mat_merge_all) <- sapply(1:ncol(X_res_mat_merge_all), function(x) paste("X.", x, sep = ""))
    
    {mod_ms <- glmnet::cv.glmnet(x = as.matrix(X_res_mat_merge_all),
                                 y = as.matrix(Y.ipd),
                                 alpha = alpha.glmnet,
                                 lambda = NULL,
                                 standardize = F,
                                 intercept = FALSE,
                                 standardize.response = F)
      
      # Remove the first entry that corresponds to the intercept (which is empty)
      tau_beta <- as.vector(t(coef(mod_ms, s = "lambda.min")[-1]))
      # Store coefficient estimates in matrix
      tau_beta_matrix=matrix(tau_beta,nrow=p)
      beta.ms.all=as.matrix(p_pred_mat)%*%t(tau_beta_matrix)
      
      beta.ms.all.list=c()
      id_list <- lapply(seq_along(nk), function(i) {
        start <- if (i == 1) 1 else sum(nk[1:(i-1)]) + 1
        end <- sum(nk[1:i])
        seq(start, end)
      })
      for (k in 1:K){
        beta.ms.all.list[[k]]=beta.ms.all[id_list[[k]],]
      }
      # Find the maximum number of rows
      max_rows <- max(sapply(beta.ms.all.list, nrow))
      # Function to pad a matrix to max_rows
      pad_matrix <- function(mat, max_rows) {
        n <- nrow(mat)
        if (n < max_rows) {
          # Pad with NA rows
          pad <- matrix(NA, nrow = max_rows - n, ncol = ncol(mat))
          mat <- rbind(mat, pad)
        }
        return(mat)
      }
      
      # Apply the padding
      beta.ms.padded <- lapply(beta.ms.all.list, pad_matrix, max_rows = max_rows)
      # Combine the padded matrices
      beta.ms <- do.call(cbind, beta.ms.padded)
    }
    
    tau.ms <- ms_rlearner(X = X_res_mat_merge_all, 
                          Y = Y.ipd, 
                          alpha = alpha.glmnet, K = K, p = p-1, 
                          p_mat = p_pred_mat_test,
                          newX = x.test.ms.all)
    
    MSE.ms = mean((tau.ms - tau.test.true)^2)
    
    result.MS=data.frame("seed"=seed,"method"="MS",
                         "ARI"=NA,
                         "MSE"=NA,
                         "MSE.TAU"=MSE.ms,
                         "BIAS"=NA)
    result=rbind(result,result.MS)
    ######################################################################
    #########  save result                                   #############
    ######################################################################
    #save beta
    rownames(beta.est)[1:7]<-c("true","LOC","OR","FR-IPD","DFR-CD","KM","MA")

    #save other results
    result.data<-data.frame("seed"=seed,
                            "n"=n, "c"=c, "p"=p, "K"=K,
                            "random"=random, "discrete"=discrete,"rho"=rho,
                            "alpha"=alpha,"alpha.glmnet"=alpha.glmnet, "r"=r, "s"=s, "threshold" = threshold,
                            "n_digit"=n_digit,
                            "prepare_time.OR"=prepare.time.or,
                            "prepare_time.FRIPD"=prepare.time.ipd,
                            "prepare_time.DFRCD"=prepare.time.cd,
                            "model_time.OR"=model.time.or,
                            "model_time.FRIPD"=model.time.IPD,
                            "model_time.DFRCD"=model.time.CD,
                            "converge.OR"=converged.or,
                            "converged.FRIPD"=converged.IPD, 
                            "converged.DFRCD"=converged.CD,
                            "iter.OR"=iter.or,
                            "iter.FRIPD"=iter.IPD, 
                            "iter.DFRCD"=iter.CD)
    
    result.data.all=rbind(result.data.all,result.data)
    beta.est.all=rbind(beta.est.all,beta.est)
    })
}#end of whole loop

#average ARI
result.ARI<-aggregate(x=result$ARI,
                   by=list(result$method),
                   FUN = mean) 

#average MSE
result.MSE<-aggregate(x=result$MSE,
                      by=list(result$method),
                      FUN = mean) 

#average prediction error
result.Prediction<-aggregate(x=result$MSE.TAU,
                      by=list(result$method),
                      FUN = mean) 

#average mean absolute bias (MAB)
result.BIAS<-aggregate(x=result$BIAS,
                             by=list(result$method),
                             FUN = mean) 

# calculate SE below
{method.list <- c("LOC", "OR", "FR-IPD", "DFR-CD", "KM", "MA")
n_method <- length(method.list)
n_rep <- length(seedlist)
true.ind <- (0:(n_rep - 1)) * (n_method + 1) + 1

# Create a named list of method indices
method_inds <- setNames(
  lapply(1:n_method, function(i) (0:(n_rep - 1)) * (n_method + 1) + (i + 1)),
  method.list
)

# Initialize a list to store SEs by method
se_method_list <- setNames(vector("list", n_method), method.list)

# Function to compute weighted SE for one method and covariate
get_weighted_se <- function(estimates, true_vals) {
  unique_vals <- unique(true_vals)
  weighted_ses <- sapply(unique_vals, function(val) {
    group_vals <- estimates[true_vals == val]
    rep_count <- length(group_vals)
    sd(group_vals, na.rm = TRUE) / sqrt(rep_count)
  })
  weights <- as.numeric(table(true_vals))
  weighted.mean(weighted_ses, w = weights)
}

# Main loop
for (j in 1:p) {
  index.j <- get_covariate_indices(j, K, p)
  true.j <- c(beta.est.all[true.ind, index.j])
  
  for (method in method.list) {
    est.j <- c(beta.est.all[method_inds[[method]], index.j])
    wse <- get_weighted_se(estimates = est.j, true_vals = true.j)
    se_method_list[[method]] <- c(se_method_list[[method]], wse)
  }
}

# Get mean SE across covariates for each method
result.se <- sapply(se_method_list, mean)
}

#save results
save.image(file = "my_workspace.RData")
