######################################################################
#########  Simulation code for testing                   #############
######################################################################
rm(list = ls())
setwd(getwd())
######################################################################
#########  Load packages                                 #############
######################################################################
library(Matrix)
library(stringr)
library(ggplot2)
######################################################################
#########  parameters                                    #############
######################################################################
c=0 #level of heterogeneity
n=100 #sample size for each site
p=4 #number of covariates
number_studies=10
K=number_studies
nk=rep(n,K)
# nk=sample(c(200,500),10,replace = T)
f=5 #folds of cross fitting
M=2 #number of treatment
random=F #random treatment assignment or not
discrete=T #discrete grouping or not
N=sum(nk)
rho=1
alpha=1
r=0.5
s=0.5
threshold = 1e-8
n_digit=3
maxiter=2500
vk<-rep(1,K)
plotname=""
lambda <- 10^(seq(-6,0,length.out=100))
date<-"240228"

######################################################################
#########             functions                          #############
######################################################################

## generate true coefficients for all sites
gencoef<-function(c=1,number_studies=10,discrete=T,p=4,seed=NA){
  if (!is.na(seed)){
    set.seed(seed)
  }
  beta_tau_null<-matrix(nrow = p,ncol = number_studies)
  beta_tau_null[,]<-t(c(1.5,2,-3,1))
  beta_m_null<-matrix(nrow = p,ncol = number_studies)
  beta_m_null[,]<-t(c(1,0.5,-1.5,2))
  beta_e_null<-matrix(nrow = p,ncol = number_studies)
  beta_e_null[,]<-t(c(0.2,0.1,-0.3,0.2))
  g1k<-rbinom(n=number_studies, size=1, prob=0.5)#grouping for beta1
  g2k<-rbinom(n=number_studies, size=1, prob=0.5)#grouping for beta2
  if (discrete){
    beta_tau_null[2,]<-beta_tau_null[2,]+c*g1k
    beta_tau_null[3,]<-beta_tau_null[3,]+c*g2k
  }else {
    beta_tau_null[2,]<-beta_tau_null[2,]+c*rbinom(n=number_studies, size=1, prob=0.5)+runif(number_studies,-0.1,0.1)
    beta_tau_null[3,]<-beta_tau_null[3,]+c*rbinom(n=number_studies, size=1, prob=0.5)+runif(number_studies,-0.1,0.1)
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

## generate simulation data for all sites
generatedata<-function(nk,vk,p,number_studies,M,random=F,beta_m,beta_e,beta_tau){
  K=number_studies
  for (k in 1:K){
    x<-c()
    prob<-c()
    ecolumn<-c()
    n<-nk[k]
    v<-vk[k]
    x_temp<-rep(1,n)#intercept, x_0
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
      ecolumn[e]<-sample(0:(M-1),1,replace=T,prob = c(1-prob[e],prob[e]))
    }
    x_temp<-cbind(x_temp,ecolumn)
    x<-as.matrix(x_temp)
    colnames(x)<-c(paste(paste("x",k,sep = ""),0:(p-1),sep = "_"),paste("trt",k,sep = "_"))
    #y<-log(1+exp(x[,2]+x[,3]+x[,4]))*2+diag(as.list(x[,p+2]-prob))%*%(x[,1:(p+1)]%*%beta_tau[k,])+rnorm(nrow(x),sd=sigma) 
    #for additional simulation where m(x) can not be estimated well
    y<-x[,1:p]%*%beta_m[k,]+diag(as.list(x[,p+1]-prob))%*%(x[,1:p]%*%beta_tau[k,])+rnorm(nrow(x),sd=v)
    data<-cbind(y,x)
    colnames(data)[1]<-"y"
    data<-data.frame(data)
    data.sim[[k]]<-data
  }
  return(data.sim)
}

## generate ids for cross fitting or cross validation
cf_index<-function(n,f){
  ids<-c(1:n)
  indexes<-c()
  for (i in 1:(f-1)) {
    index_for_ids<-sample(1:length(ids),size = round(n/f))
    indexes[[i]]<-ids[index_for_ids]
    ids<-ids[-index_for_ids]
  }
  indexes[[f]]<-ids
  return(indexes)
}

#DF R-learner for IPD
DFIPD <- function(X,Y,rho,r,s,lambda,alpha=0, nk, vk, n_digit=2,
                  number_studies,number_treatments,maxiter=5000,threshold=1e-6){
  p=ncol(X)/number_studies
  N<-nrow(X)
  #initialize beta
  beta <- solve(crossprod(X)) %*% crossprod(X,Y)
  rownames(beta)<-colnames(X)
  #generate the D matrix
  if (number_studies==1){
    #error message
    stop("number_studies should be at least 2")
    D<-diag(p)
  } else{
    #generate D matrix for multiple studies
    D <- matrix(0, nrow=(number_studies)*p, ncol=(number_studies)*p) 
    for (j in 1:p){
      beta.j<-beta[c((p)*(c(0:(number_studies-1)))+j),1]
      orderU<-order(beta.j)
      orderV<-order(abs(beta.j))
      for (i in 1:(number_studies-1)){
        D[1+(j-1)*number_studies,(orderV[1]-1)*p+j]<-1 #the first line indicates smallest 
        D[i+1+(j-1)*number_studies,(orderU[i]-1)*p+j]<-(-1)
        D[i+1+(j-1)*number_studies,(orderU[i+1]-1)*p+j]<-1
      }
    }
  }
  #initialize theta
  theta<-as.vector(D%*%beta) #only for calculating W
  #generate W matrix for weights 
  W<-rep(0,(number_studies)*p)
  for (j in 1:p){
    W[1+(j-1)*number_studies]<-alpha*abs(theta[1+(j-1)*number_studies])^(-r)
    thetarange<-sum(theta[(2+(j-1)*number_studies):(j*number_studies)])
    for (i in 2:number_studies){
      W[i+(j-1)*number_studies]<-abs(thetarange)^(-s)*abs(theta[i+(j-1)*number_studies])^(-r)
    }
  }
  W<-replace(W,W==Inf,1e6)
  #W<-diag(W)
  Vinv<-c()
  for(i in 1:number_studies){
    Vinv<-c(Vinv,rep(vk[i],nk[i])^-1)
  }
  Vinv<-diag(Vinv)
  
  theta <- matrix(0, ncol=1, nrow=nrow(D)) #don't use D*beta as initial value. 
                                           #Other wise wont update
  kappa <- matrix(0, ncol=1, nrow=nrow(D)) 

  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  converge<-0
  p1<-solve(t(X)%*%Vinv%*%X+N*rho*t(D)%*%D)
  p2<-t(X)%*%Vinv%*%Y
  # ADMM updates
  for (t in 2:maxiter){
    # beta.new<-t(solve(t(X)%*%Vinv%*%X+N*rho*t(D)%*%D) %*%
    #               (t(X)%*%Vinv%*%Y + N*rho*t(D)%*%(t(theta)-t(kappa))))
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
  for (i in 1:number_studies){
    n2ll[i]<-nk[i]*log(sum(residual[stp:(stp+nk[i]-1)]^2)/(nk[i]))
    stp<-stp+nk[i]
  }
  for (j in 1:p){
    beta.j<-beta.hat.round[c((p)*(c(0:(number_studies-1)))+j)]
    df[j]<-length(unique(beta.j))
    combin<-combin+choose(number_studies,df[j])
  }
  #EBIC
  bic<-sum(n2ll)+sum(df)*log(N)+2*log(combin)
  #AIC
  aic<-sum(n2ll)+2*sum(df)
  
  #transform beta in 
  beta.hat.matrix<-data.frame(matrix(beta.hat,ncol=p,byrow = T))
  rownames(beta.hat.matrix)<-paste("study",1:number_studies,sep = "_")
  colnames(beta.hat.matrix)<-c("Intercept",paste0("x",1:(p-1),sep=""))
  beta.hat.matrix<-data.matrix(beta.hat.matrix,rownames.force=T)
  result <- list("beta.hat.matrix"=beta.hat.matrix, "beta.hat"=beta.hat,
                 "conv_message"=message, "iter"=t,"BIC"=bic,"AIC"=aic,
                 "converge"=converge)   
  return(result)
}

#DF R-learner for CD
DFCD <- function(sigma_n_inv,b,rho,r,s,lambda, n_digit=2, nk,alpha=0,
                 number_studies,maxiter=5000,threshold=1e-12){
  p=ncol(sigma_n_inv)/number_studies
  N=sum(nk)
  #initialize beta matrix
  beta <- matrix(b, ncol=1, nrow=ncol(sigma_n_inv)) 
  rownames(beta)<-colnames(sigma_n_inv)
  #generate the D matrix
  if (number_studies==1){
    #error message
    stop("number_studies should be at least 2")
    D<-diag(p)
  } else{
    #generate D matrix for multiple studies 
    D <- matrix(0, nrow=(number_studies)*p, ncol=(number_studies)*p) 
    for (j in 1:p){
      beta.j<-beta[c((p)*(c(0:(number_studies-1)))+j),1]
      orderU<-order(beta.j)
      orderV<-order(abs(beta.j))
      for (i in 1:(number_studies-1)){
        D[1+(j-1)*number_studies,(orderV[1]-1)*p+j]<-1
        D[i+1+(j-1)*number_studies,(orderU[i]-1)*p+j]<-(-1)
        D[i+1+(j-1)*number_studies,(orderU[i+1]-1)*p+j]<-1
      }
    }
  }
  #initialize theta
  theta<-as.vector(D%*%beta) # only for calculating W
  #generate W matrix for weights
  W<-rep(0,(number_studies)*p)
  for (j in 1:p){
    W[1+(j-1)*number_studies]<-alpha*abs(theta[1+(j-1)*number_studies])^(-r)
    thetran<-sum(theta[(2+(j-1)*number_studies):(j*number_studies)])
    for (i in 2:number_studies){
      W[i+(j-1)*number_studies]<-abs(thetran)^(-s)*abs(theta[i+(j-1)*number_studies])^(-r)
    }
  }
  W<-replace(W,W==Inf,1e6)
  #W<-diag(W)
  theta <- matrix(0, ncol=1, nrow=nrow(D)) #don't use D*beta as initial value. 
                                           #Other wise wont update
  kappa <- matrix(0, ncol=1, nrow=nrow(D))  
  
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
  beta.hat=beta.new   #changed
  beta.hat.round<-round(beta.hat,n_digit)
  #elements for BIC
  df<-c()
  combin<-0
  for (j in 1:p){
    beta.j<-beta.hat.round[c((p)*(c(0:(number_studies-1)))+j)]
    df[j]<-length(unique(beta.j))
    combin<-combin+choose(K,df[j])
  }
  #EBIC
  bic<-t(b-beta.hat)%*%sigma_n_inv%*%(b-beta.hat)+sum(df)*log(N)+2*log(combin)
  #AIC
  aic<-t(b-beta.hat)%*%sigma_n_inv%*%(b-beta.hat)+2*sum(df)
  #transform beta in 
  beta.hat.matrix<-data.frame(matrix(beta.hat,ncol=p,byrow = T))
  rownames(beta.hat.matrix)<-paste0("study",1:number_studies,sep="")
  colnames(beta.hat.matrix)<-c("Intercept",paste0("x",1:(p-1),sep=""))
  beta.hat.matrix<-data.matrix(beta.hat.matrix,rownames.force=T)
  result <- list("beta.hat.matrix"=beta.hat.matrix, "beta.hat"=beta.hat,
                 "conv_message"=message, "iter"=t,"BIC"=bic,"AIC"=aic,
                 "converge"=converge)  
  return(result)
}

plotpath<-function(betamat,lambda,name){
  data.gam<-as.data.frame(t(betamat))
  data.gam$lambda<-lambda
  data.gam <- reshape(data.gam,              
                      direction = "long",
                      varying = names(data.gam)[1:(ncol(data.gam)-1)],
                      v.names = "y",
                      idvar = c("lambda"),
                      timevar = "Variable",
                      times = names(data.gam)[1:(ncol(data.gam)-1)])
  a<-unlist(str_split(data.gam$Variable, "_"))
  # data.gam$study<-a[sequence(nvec = nrow(data.gam),from=1,by=2)]
  data.gam$study<-a[seq(from=1,to=nrow(data.gam)*2-1,by=2)]
  data.gam$study<-sub('.', 'study', data.gam$study)
  # data.gam$covariate<-a[sequence(nvec = nrow(data.gam),from=2,by=2)]
  data.gam$covariate<-a[seq(from=2,to=nrow(data.gam)*2,by=2)]
  data.gam$covariate<-paste("x",data.gam$covariate,sep = "")
  
  #name.p<-paste(paste("plot",paste("n",n,sep = ""),setting,penalty,seed,date,sep = "_"),"pdf",sep = ".")
  # name.p<-paste(paste("plot",paste("n",n,sep = ""),setting,penalty,seed,date,sep = "_"),"pdf",sep = ".")
  pdf(file=name)
  solutionpath<-ggplot(data.gam, aes(x = lambda, y = y)) +
    geom_line(aes(group = Variable, color = study,linetype = covariate)) + # CHANGED
    scale_x_log10() +
    geom_vline(xintercept = lambda[which.min(bicmat[2,])], linetype="dotted", color = "red", size=1.5) +
    #scale_x_log10() +
    xlab("Lambda") + 
    ylab("Coefficient") +
    theme_bw() +
    theme(legend.key.width = unit(3,"lines"))
  print(solutionpath)
  dev.off()
}

######################################################################
#########         loop over seed                         #############
######################################################################
# seed.seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
seed.seed <- 1
## number 100 can be changed
seedlist<-(1+(seed.seed-1)*100):(seed.seed*100)
for (seed in seedlist){
  set.seed(seed)
  ######################################################################
  #########  store the data, output                        #############
  ######################################################################
  # beta.est<-c()#store fitted and true beta for all estimators
  # converged<-c()
  # iteration.number<-c()
  # prepare.time<-NA
  # model.time<-NA
  # lambda.best<-c()

  
  ######################################################################
  #########  generate coef                                 #############
  ######################################################################
  simcoef<-gencoef(c=c,number_studies=number_studies,discrete=discrete,p=p,seed=NA)
  beta_m=simcoef$beta_m_null
  beta_e=simcoef$beta_e_null
  beta_tau=simcoef$beta_tau_null
  beta.est<-matrix(t(beta_tau),nrow=1)#first line is true value
  ######################################################################
  #########  generate data                                 #############
  ######################################################################
  data.sim<-c()
  data.sim<-generatedata(nk=nk,vk=vk,p=p,number_studies=number_studies,M=M,
                         random=F,beta_e=beta_e,beta_m=beta_m,beta_tau=beta_tau)
  
  ######################################################################
  #########          prepare data                          #############
  ######################################################################
  ## Oracle 
  start_time <- Sys.time()
  X.s1<-c()
  Y.s1<-c()
  X_tilde_true<-c() #X_tilde
  Y_tilde_true<-c() #Y_tilde
  Colnames<-c()
  v.s1<-c()
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
    #Y_tilde_true[[k]]<-y_tilde_true
    X_tilde_true[[k]]<-x_tilde_true
    Colnames<-c(Colnames,colnames(X_tilde_true[[k]]))
    Y.s1<-c(Y.s1,y_tilde_true)
    fit.s1 <- glm(y_tilde_true~x_tilde_true-1,family="gaussian")
    v.s1<-c(v.s1,sum(fit.s1$residuals^2)/n)
  }
  X.s1<-as.matrix(bdiag(X_tilde_true))
  colnames(X.s1)<-Colnames
  end_time <- Sys.time()
  prepare.time.s1<-end_time - start_time
  
  ## IPD
  start_time <- Sys.time()
  X_tilde_est<-c() #X_tilde
  Y_tilde_est<-c() #Y_tilde
  X.s3<-c()
  Y.s3<-c()
  Colnames<-c()
  v.s3<-c()
  for (k in 1:K){
    n<-nk[k]
    data<-as.matrix(data.sim[[k]])
    y<-data[,1]
    x<-data[,-1]
    indexes<-cf_index(n=n,f=f)
    m_hat<-c()
    e_hat<-c()
    for (i in 1:f) {
      index.all<-c(1:n)
      index.test<-indexes[[i]]
      index.train<-index.all[-index.test]
      data.train<-data.frame(data[index.train,])
      data.test<-data.frame(data[index.test,])
      #######cross fitting m
      model.m<-glm(y~.-1,data=data.train[,-ncol(data.train)],family="gaussian")
      m_hat[index.test]<-predict(model.m, data.test, type="response")
      #######cross fitting e
      model.e<-glm(as.formula(paste(colnames(data.train)[p+2], "~ .-1")),
                   data=data.train[,-1],family="binomial")
      e_hat[index.test]<-predict(model.e, data.test, type="response")
    }
    y_tilde_est<-y-m_hat
    x_tilde_est<-diag(x[,p+1]-e_hat)%*%x[,1:(p)]
    Y.s3<-c(Y.s3,y_tilde_est)
    X_tilde_est[[k]]<-x_tilde_est
    Colnames<-c(Colnames,colnames(x_tilde_est))
    fit.s3 <- glm(y_tilde_est~x_tilde_est-1,family="gaussian")
    v.s3<-c(v.s3,sum(fit.s3$residuals^2)/n)
  }
  X.s3<-as.matrix(bdiag(X_tilde_est))
  colnames(X.s3)<-Colnames
  end_time <- Sys.time()
  prepare.time.s3<-(end_time - start_time)
  
  ## CD
  start_time <- Sys.time()
  sigma_n_inv_est<-c() #sigma_n inverse matrix
  beta.s4<-c()
  Colnames<-c()
  for (k in 1:K){
    n<-nk[k]
    data<-as.matrix(data.sim[[k]])
    y<-data[,1]
    x<-data[,-1]
    indexes<-cf_index(n=n,f=f)
    m_hat<-c()
    e_hat<-c()
    for (i in 1:f) {
      index.all<-c(1:n)
      index.test<-indexes[[i]]
      index.train<-index.all[-index.test]
      data.train<-data.frame(data[index.train,])
      data.test<-data.frame(data[index.test,])
      #######cross fitting m
      model.m<-glm(y~.-1,data=data.train[,-ncol(data.train)],family="gaussian")
      m_hat[index.test]<-predict(model.m, data.test, type="response")
      #######cross fitting e
      model.e<-glm(as.formula(paste(colnames(data.train)[p+2], "~ .-1")),data=data.train[,-1],family="binomial")
      e_hat[index.test]<-predict(model.e, data.test, type="response")
    }
    y_tilde_est<-y-m_hat
    x_tilde_est<-diag(x[,p+1]-e_hat)%*%x[,1:(p)]
    x_data_est<-cbind(y_tilde_est,x_tilde_est)
    colnames(x_data_est)[1]<-"y"
    x_data_est<-data.frame(x_data_est)
    fit.s4 <- glm(y~.-1,data=x_data_est,family="gaussian")
    beta.s4<-c(beta.s4,fit.s4$coefficients)
    v<-solve(solve(t(x_tilde_est)%*%x_tilde_est)*sum(fit.s4$residuals^2)/n)#sigma_n_inverse
    Colnames<-c(Colnames,colnames(x_tilde_est))
    sigma_n_inv_est[[k]]<-v
  }
  sigma_n_inv_est.s4<-as.matrix(bdiag(sigma_n_inv_est))
  colnames(sigma_n_inv_est.s4)<-Colnames
  end_time <- Sys.time()
  prepare.time.s4<-(end_time - start_time)
  
  ######################################################################
  #########    estimate LOC                                #############
  ######################################################################
  beta.local<-summary(lm(Y.s3~X.s3-1))$coefficients[,1]
  beta.est <- rbind(beta.est,beta.local)
  # names(beta.local)<-Colnames
  colnames(beta.est)<-Colnames
  ######################################################################
  #########    estimate oracle                             #############
  ######################################################################
  betamat <- matrix(rep(NA, length(lambda)*number_studies*p), ncol=length(lambda))
  bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  bicmat[1,]  <- lambda
  convergedmat[1,]  <- lambda
  itermat[1,]  <- lambda
  start_time <- Sys.time()
  for (i in 1:length(lambda)) {
    dfmodel <- DFIPD(X=X.s1,Y=Y.s1,rho=rho,r=r,s=s,lambda=lambda[i], number_studies=number_studies, 
                     threshold=threshold, nk=nk,vk=v.s1, alpha=alpha, n_digit=n_digit,maxiter = maxiter)
    betamat[,i]  <- dfmodel$beta.hat #store estimated beta
    bicmat[2,i]  <- dfmodel$BIC #store EBIC
    convergedmat[2,i]  <- dfmodel$converge #store convergence
    itermat[2,i]  <- dfmodel$iter #store iteration
  }
  end_time <- Sys.time()
  model.time.oracle<-(end_time - start_time)/length(lambda)
  rownames(betamat)<-Colnames
  index.lambda.best<-which.min(bicmat[2,])
  beta.oracle<-betamat[,index.lambda.best]
  beta.est<-rbind(beta.est,beta.oracle)
  lambda.best.oracle<-lambda[index.lambda.best]
  converged.oracle<-convergedmat[2,index.lambda.best]
  iter.oracle<-itermat[2,index.lambda.best]
  #plotpath(betamat=betamat,lambda=lambda,name="XX.pdf")
  ######################################################################
  #########    estimate IPD                                #############
  ######################################################################
  betamat <- matrix(rep(NA, length(lambda)*number_studies*p), ncol=length(lambda))
  bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  bicmat[1,]  <- lambda
  convergedmat[1,]  <- lambda
  itermat[1,]  <- lambda
  start_time <- Sys.time()
  for (i in 1:length(lambda)) {
    dfmodel <- DFIPD(X=X.s3,Y=Y.s3,rho=rho,r=r,s=s,lambda=lambda[i], number_studies=number_studies, 
                     threshold=threshold, nk=nk,vk=v.s3, alpha=alpha, n_digit=n_digit,maxiter = maxiter)
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
  #plotpath(betamat=betamat,lambda=lambda,name="solution2.pdf")
  
  ######################################################################
  #########    estimate  CD                                #############
  ######################################################################
  betamat <- matrix(rep(NA, length(lambda)*number_studies*p), ncol=length(lambda))
  bicmat <- matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  convergedmat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  itermat<-matrix(rep(NA, length(lambda)*2), ncol=length(lambda))
  bicmat[1,]  <- lambda
  convergedmat[1,]  <- lambda
  itermat[1,]  <- lambda
  start_time <- Sys.time()
  for (i in 1:length(lambda)) {
    dfmodel <- DFCD(sigma_n_inv=sigma_n_inv_est.s4,b=beta.s4,rho=rho,r=r,s=s,nk=nk,lambda=lambda[i], number_studies=number_studies, 
                    threshold=threshold, alpha=alpha, n_digit=n_digit,maxiter = maxiter)
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
  #plotpath(betamat=betamat,lambda=lambda,name="solution3.pdf")
  
  ######################################################################
  #########    estimate K-Means                            #############
  ######################################################################
  betalist<-c()
  for (j in 1:p){
    beta.j<-beta.local[p*((1:number_studies)-1)+j] 
    #calculate center=1 first
    betamat<-rep(kmeans(beta.j, centers = 1, nstart=10)$centers,number_studies)
    # calculate center=2 to 9, because center=10==LOC
    for (i in 2:(number_studies-1)){
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
  combnmat<-expand.grid(rep(list(1:number_studies), p))
  calculatebic<-function(x){
    comb<-combnmat[x,]
    beta<-c()
    for (j in 1:p){
      # temp<-unlist(betalist[[j]])
      beta<-c(beta,betalist[[j]][unlist(comb[j]),])
    }
    beta<-matrix(t(matrix(beta,nrow=number_studies)),ncol =1)
    beta.round<-round(beta,n_digit)
    df<-c()
    combin<-0
    for (j in 1:p){
      beta.j<-beta.round[c((p)*(c(0:(number_studies-1)))+j)]
      df[j]<-length(unique(beta.j))
      combin<-combin+choose(number_studies,df[j])
    }
    bic<-t(beta-beta.local)%*%sigma_n_inv_est.s4%*%(beta-beta.local)+sum(df)*log(N)+2*log(combin)
    return(bic)
  }
  biclist<-sapply(1:nrow(combnmat),calculatebic)
  comb<-combnmat[which.min(biclist),]
  beta.KM<-c()
  for (j in 1:p){
    # temp<-unlist(betalist[[j]])
    beta.KM<-c(beta.KM,betalist[[j]][unlist(comb[j]),])
  }
  beta.KM<-t(matrix(t(matrix(beta.KM,nrow=number_studies)),ncol = 1))
  beta.est<-rbind(beta.est,beta.KM)
  
  ######################################################################
  #########    estimate MA                                 #############
  ######################################################################
  beta.MA<-c()
  for (j in 1:p){
    beta.j<-beta.local[c((p)*(c(0:(number_studies-1)))+j)]
    beta.MA<-c(beta.MA,rep(mean(beta.j),number_studies))
  }
  beta.MA<-t(matrix(t(matrix(beta.MA,nrow=number_studies)),ncol = 1))
  beta.est<-rbind(beta.est,beta.MA)
  ######################################################################
  #########  save result                                   #############
  ######################################################################
  #save beta
  rownames(beta.est)<-c("true","local","oracle","IPD","CD","KM","MA")
  # beta.data<-data.frame("true"=matrix(t(beta_tau),ncol=1),"local"=beta.local,"oracle"=beta.oracle,
  #                     "IPD"=beta.IPD,"CD"=beta.CD, "KM"=t(beta.KM),"MA"=t(beta.MA))
  name1<-paste(paste("beta",paste("n",n,sep = ""),seed,date,sep = "_"),"csv",sep = ".")
  # name1<-paste(paste("beta",paste("n","mix",sep = ""),setting,penalty,seed,date,sep = "_"),"csv",sep = ".")
  write.csv(beta.est,file=name1)
  #save time
  time.data<-data.frame("seed"=seed,"prepare_time.oracle"=prepare.time.s1,"prepare_time.IPD"=prepare.time.s3,
                        "prepare_time.CD"=prepare.time.s4,"model_time.oracle"=model.time.oracle,
                        "model_time.IPD"=model.time.IPD,"model_time.CD"=model.time.CD)
  name2<-paste(paste("time",paste("n",n,sep = ""),seed,date,sep = "_"),"csv",sep = ".")
  write.csv(time.data,file=name2,row.names = F)
  #save convergency
  converge.data<-data.frame("seed"=seed,"converge.oracle"=converged.oracle,
                            "converged.IPD"=converged.IPD, 
                            "converged.CD"=converged.CD)
  name3<-paste(paste("conv",paste("n",n,sep = ""),seed,date,sep = "_"),"csv",sep = ".")
  write.csv(converge.data,file=name3,row.names = F)
  #save iteration
  iter.data<-data.frame("seed"=seed,"iter.oracle"=iter.oracle,
                            "iter.IPD"=iter.IPD, 
                            "iter.CD"=iter.CD)
  name4<-paste(paste("iter",paste("n",n,sep = ""),seed,date,sep = "_"),"csv",sep = ".")
  write.csv(iter.data,file=name4,row.names = F)
  #save parameters
  para.data<-data.frame("seed"=seed,"n"=n, "c"=c, "M"=M, "p"=p, 
                        "number_studies"=number_studies,
                        "random"=random, "discrete"=discrete,"rho"=rho,
                        "alpha"=alpha, "r"=r, "s"=s, "threshold" = threshold,
                        "n_digit"=n_digit)
  name5<-paste(paste("para",paste("n",n,sep = ""),seed,date,sep = "_"),"csv",sep = ".")
  write.csv(para.data,file=name5,row.names = F)
}#end of whole loop