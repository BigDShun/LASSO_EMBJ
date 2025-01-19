


library(glmnet)
library(relliptical)
library(caret) 
library(latex2exp) 
library(MASS)
#library("NormalBetaPrime")  

EM.NormalCens <- function(cc, x, y, lambda=0.1, alpha=1, error = 0.0001, iter.max = 100){
  ################################################################################
  ## cc is a vector nx1 of left-censoring 0 = uncensoring or 1 = censoring
  ## x is the design matrix of dimension nxp
  ## y vector of responses nx1
  ## alpha=1 Lasso
  ## alpha=0 Ridge
  ## alpha \in (0,1) elastic-net
  ################################################################################
  
  p <- ncol(x)
  n <- nrow(x)
  
  ## Initial values
  
  reg <-  glmnet(x, y, lambda=lambda,family="gaussian",intercept = F, alpha=alpha)
  beta<- reg$beta[,1] 
  mu<- x%*%beta
  sigma2 <- sum((y-mu)^2)/n                                             
  
  
  ################################################################################
  ###                                     Normal 
  ################################################################################
  
  cont <- 0
  criterio <- 1
  lkante   <- 1
  
  while(criterio > error){
    
    cont <- (cont+1)
    
    E01<-y
    E02<-y^2
    
    if(sum(cc)>0)
    { 
      mu1<- mu[cc==1]
      y1<- y[cc==1]
      np<- length(mu1)
      aux1MomW<- matrix(0,np,2)
      
      for(j in 1:np){
        A1a<- mvtelliptical(-Inf, y1[j], mu1[j], sigma2, "Normal", n=1e6)
        aux1MomW[j,]<-c(A1a$EY,A1a$EYY)
      }
      
      
      
      E01[cc==1]<- aux1MomW[,1]
      E02[cc==1]<- aux1MomW[,2]
    }
    
    
    lambda1<-lambda*sigma2
    la.eq<-  glmnet(x, E01, lambda=lambda1,family="gaussian", intercept = F, alpha=alpha)
    beta<-la.eq$beta[,1] 
    mu<-x%*%beta
    sigma2<-sum(E02-2*E01*mu+mu^2)/n		
    
    
    
    ################################################################################
    auxpdf<-dnorm(y[cc==0], mu[cc==0], sqrt(sigma2),log = TRUE)
    auxcdf<- pnorm((y[cc==1]- mu[cc==1])/sqrt (sigma2),log.p = TRUE)
    lk<-sum(auxpdf)+sum(auxcdf)
    logver<-lk
    ################################################################################
    
    criterio <- abs((lk/lkante-1))
    
    lkante <- lk
    
    if (cont==iter.max){
      criterio <- error/10
    }
    
    
  }
  #print(cont)
  ychap<-E01 
  teta_novo<-matrix(c(beta,sigma2),ncol=1)
  
  nf<-(p-length(which(beta ==0)))+1
  
  
  bic<-  (-2*lk + log(n)*nf)/n
  
  return(list(theta=teta_novo, iter=cont,logver=logver,BIC=bic,yhat=ychap))	
}


BJ.NormalCens <- function(cc, x, y, lambda=0.1, alpha=1, error = 0.0001, iter.max = 100){
  ################################################################################
  ## cc is a vector nx1 of left-censoring 0 = uncensoring or 1 = censoring
  ## x is the design matrix of dimension nxp
  ## y vector of responses nx1
  ## alpha=1 Lasso
  ## alpha=0 Ridge
  ## alpha \in (0,1) elastic-net
  ################################################################################
  
  p <- ncol(x)
  n <- nrow(x)
  
  ## Initial values
  
  reg <- glmnet(x, y, lambda=lambda,family="gaussian",intercept = F, alpha=alpha)
  beta<- reg$beta[,1] 
  mu<- x%*%beta
  sigma2 <- sum((y-mu)^2)/n                                             
  
  ################################################################################
  ###                                     Normal 
  ################################################################################
  
  cont <- 0
  criterio <- 1
  lkante   <- 1
  
  while(criterio > error){
    
    cont <- (cont+1)
    
    E01<-y
    
    if(sum(cc)>0)
    { 
      mu1<- mu[cc==1]
      y1<- y[cc==1]
      np<- length(mu1)
      aux1MomW<- matrix(0,np,1)
      
      for(j in 1:np){
        A1a<- mvtelliptical(-Inf, y1[j], mu1[j],sigma2, "Normal", n=1e6)
        aux1MomW[j]<-c(A1a$EY)
      }
      
      
      
      E01[cc==1]<- aux1MomW
      
    }
    
    la.eq<-  glmnet(x, E01, lambda=lambda,family="gaussian", intercept = F, alpha=alpha)
    beta<-la.eq$beta[,1] 
    mu<-x%*%beta
    sigma2<-sum((E01-mu)^2)/n		
    
    
    
    ################################################################################
    auxpdf<-dnorm(y[cc==0], mu[cc==0], sqrt(sigma2),log = TRUE)
    auxcdf<- pnorm((y[cc==1]- mu[cc==1])/sqrt (sigma2),log.p = TRUE)
    lk<-sum(auxpdf)+sum(auxcdf)
    logver<-lk
    ################################################################################
    
    criterio <- abs((lk/lkante-1))
    
    lkante <- lk
    
    if (cont==iter.max){
      criterio <- error/10
    }
    
    
  }
  #print(cont)
  ychap<-E01 
  teta_novo<-matrix(c(beta,sigma2),ncol=1)
  
  nf<-(p-length(which(beta ==0)))+1
  
  
  bic<-  (-2*lk + log(n)*nf)/n
  
  return(list(theta=teta_novo, iter=cont,logver=logver,BIC=bic,yhat=ychap))	
}


EM.NormalCens_Adaptive <- function(cc, x, y,alpha = 1, gamma = 1,error = 1e-5,iter.max = 100,nfolds = 5){
  ################################################################################
  # cc         : censoring indicator (length n); 1 => left-censored, 0 => uncensored
  # x          : design matrix (n x p)
  # y          : response vector (length n)
  # alpha      : 0=Ridge, 1=Lasso, or anything in (0,1) => elastic net
  # gamma      : exponent for adaptive weights => w_j = 1 / (|beta_j| + eps)^gamma
  # error      : convergence tolerance (change in log-likelihood)
  # iter.max   : maximum EM iterations
  # nfolds     : number of folds for cv.glmnet to select lambda
  #
  # Returns a list with:
  #   theta  : final parameters (beta, sigma^2)
  #   iter   : number of EM iterations
  #   logver : final log-likelihood
  #   BIC    : approximate BIC measure
  #   yhat   : final E[Y]
  ################################################################################
  
  p <- ncol(x)
  n <- nrow(x)
  
  # 1) Initial non-adaptive fit (to get beta_init for weights)
  init_fit <- glmnet(x, y, lambda=0.01, family="gaussian",
                     intercept=FALSE, alpha=alpha)
  beta_init <- as.numeric(init_fit$beta)
  
  # 2) Compute one-shot adaptive weights
  eps <- 1e-6
  w <- 1 / (abs(beta_init) + eps)^gamma
  
  # 3) Initialize EM
  mu <- x %*% beta_init
  sigma2 <- mean((y - mu)^2)
  cont <- 0
  criterio <- 1
  lkante <- 1
  
  while (criterio > error && cont < iter.max){
    cont <- cont + 1
    
    #-----------------------
    # E-step
    #-----------------------
    E01 <- y
    E02 <- y^2
    
    if (any(cc == 1)) {
      mu_cens <- mu[cc == 1]
      y_cens <- y[cc == 1]
      np <- length(mu_cens)
      aux <- matrix(0, np, 2)
      
      for(j in 1:np){
        # Must define or have mvtelliptical for truncated normal integrals:
        out_j <- mvtelliptical(-Inf, y_cens[j],
                               mu_cens[j], sigma2,
                               "Normal", n=1e5)
        aux[j, ] <- c(out_j$EY, out_j$EYY)
      }
      E01[cc == 1] <- aux[,1]
      E02[cc == 1] <- aux[,2]
    }
    
    #-----------------------
    # M-step
    #   - We do not update w anymore! (One-shot)
    #-----------------------
    # Fit model with penalty.factor=w and cross-validate for lambda
    mod.eta <- cv.glmnet(x, E01, family="gaussian", intercept=FALSE,
                         alpha=alpha, penalty.factor=w,
                         nfolds=nfolds)
    lambda.sel <- mod.eta$lambda.1se
    
    fit <- glmnet(x, E01, lambda=lambda.sel, family="gaussian",
                  intercept=FALSE, alpha=alpha, penalty.factor=w)
    beta <- as.numeric(fit$beta)
    beta[abs(beta) < 1e-8] <- 0
    
    mu <- x %*% beta
    sigma2 <- mean(E02 - 2*E01*mu + mu^2)
    
    #-----------------------
    # Log-likelihood
    #-----------------------
    auxpdf <- dnorm(y[cc==0], mean=mu[cc==0], sd=sqrt(sigma2), log=TRUE)
    auxcdf <- pnorm((y[cc==1] - mu[cc==1]) / sqrt(sigma2), log.p=TRUE)
    lk <- sum(auxpdf) + sum(auxcdf)
    
    criterio <- abs(lk/lkante - 1)
    lkante <- lk
  }
  
  # Post-processing
  yhat <- E01
  theta <- c(beta, sigma2)
  n_nonzero <- sum(abs(beta) > 1e-8)
  nf <- n_nonzero + 1  # +1 for sigma^2
  bic <- (-2*lk + log(n)*nf)/n
  
  list(
    theta  = theta,  # (beta_1,...,beta_p, sigma^2)
    iter   = cont,
    logver = lk,
    BIC    = bic,
    yhat   = yhat
  )
}

