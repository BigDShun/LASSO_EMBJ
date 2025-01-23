###############################################################################
# Title: Censored Regression with Bootstrap Variable Selection
# Author: Dashun Liu
# Description: This script demonstrates a application using
#   the PSID1976 data from the AER package. We use EM-based approaches for
#   left-censored outcomes, plot convergence, and then perform a bootstrap
#   procedure to assess variable selection stability.
###############################################################################

# -----------------------------
# 0. Load Libraries
# -----------------------------
library(glmnet)      # Regularized regression
library(relliptical) # Functions for elliptical distributions, including EM.NormalCens
library(caret)       # Utility functions (e.g., data splitting, pre-processing)
library(latex2exp)   # LaTeX expressions in plots
library(AER)         # Contains the PSID1976 dataset
library(ggplot2)     # Plotting
library(zoo)         # Rolling functions (rollmean, etc.)

# -----------------------------
# 1. Load and Preprocess the Data
# -----------------------------
data("PSID1976")       # From AER package
wage <- PSID1976

# Create binary indicators
wage$city       <- ifelse(wage$city      == 'yes', 1, 0)
wage$college    <- ifelse(wage$college   == 'yes', 1, 0)
wage$hcollege   <- ifelse(wage$hcollege  == 'yes', 1, 0)
wage$participation <- ifelse(wage$participation == 'yes', 1, 0)

# Define response and covariates
y <- wage$wage
X <- wage[, -c(1, 7, 8)]  # Remove 'wage', 'occupation', 'participation'

                                                      
X[-c(15,17,18)]  <- scale(X[-c(15,17,18)])
X=as.matrix(X)
# Create censoring indicator vector: 0 = uncensored, 1 = censored

cc <- c(rep(0, 428), rep(1, 325))
wage$cc <- cc
N <- length(cc)

# -----------------------------
# Load neccessary functions
# -----------------------------

BJ.NormalCens <- function(cc, x, y, lambda=0.1, alpha=1, error = 0.00001, iter.max = 100){
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
        A1a<- mvtelliptical(-Inf, y1[j], mu[j], sigma2, "Normal", n=1e6)
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
  print(cont)
  ychap<-E01 
  teta_novo<-matrix(c(beta,sigma2),ncol=1)
  
  nf<-(p-length(which(beta ==0)))+1
  
  
  bic<-  (-2*lk + log(n)*nf)/n
  
  return(list(theta=teta_novo, iter=cont,logver=logver,BIC=bic,yhat=ychap))	
}

EM.NormalCens <- function(cc, x, y, alpha=1, error = 0.00001, iter.max = 100){
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
  
  reg <-  glmnet(x, y, lambda=0.01,family="gaussian",intercept = F, alpha=alpha)
  beta<- reg$beta[,1]
  mu<- x%*%beta
  sigma2 <- sum((y-mu)^2)/n
  
  lkHistory <- numeric(iter.max)
  
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
    
    
    #lambda1<-lambda*sigma2
    #la.eq<-  glmnet(x, E01, weights=rep(1/sigma2,n), lambda=lambda,family="gaussian", intercept = F, alpha=alpha)
    #beta<-la.eq$beta[,1]
    
    mod.eta <- cv.glmnet(x, E01,  family="gaussian",intercept=FALSE,
                         nfolds=20,
                         alpha=alpha)
    
    lambda.atual.eta <- mod.eta$lambda.1se
    la.eq2.eta <- glmnet(x, E01,  lambda=lambda.atual.eta,family="gaussian",
                         intercept=FALSE, alpha=alpha)
    
    beta <- matrix(as.numeric(la.eq2.eta$beta),p,1)
    beta[which(abs(beta) < 0.00000001)] <- 0
    mu<-x%*%beta
    sigma2<-sum(E02-2*E01*mu+mu^2)/n
    
    
    
    ################################################################################
    auxpdf<-dnorm(y[cc==0], mu[cc==0], sqrt(sigma2),log = TRUE)
    auxcdf<- pnorm((y[cc==1]- mu[cc==1])/sqrt (sigma2),log.p = TRUE)
    lk<-sum(auxpdf)+sum(auxcdf)
    logver<-lk
    lkHistory[cont] <- lk
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
  
  lkHistory <- lkHistory[1:cont]
  
  return(list(theta=teta_novo, iter=cont,logver=logver,BIC=bic,yhat=ychap,lkHistory=lkHistory))
}



# -----------------------------
# 2. Fit EM Model and Plot Convergence
# -----------------------------
set.seed(123)  # For reproducibility

# EM fit with alpha = 1 ( Normal distribution )
estEM <- EM.NormalCens(cc = cc, x = X, y = y,
                       alpha = 1,
                       error = 1e-6,
                       iter.max = 1000)




# Extract log-likelihood history from the EM fit
ll_history <- estEM$lkHistory

# Convergence Plots
par(mfrow = c(1, 2))

# Full iteration history
plot(ll_history, type = 'l',
     xlab = 'Iteration', ylab = 'Log-Likelihood',
     main = 'Log-Likelihood (Full Iteration)')

df_ll <- data.frame(
  Iteration = seq_along(ll_history),
  LogLik    = ll_history
)
p11 <- ggplot(df_ll, aes(x = Iteration, y = LogLik)) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Log-Likelihood"
  ) +
  theme_minimal()



#### BIC plot for BJ

aa<-seq(0.0001,0.3,length.out=50)

BICc<-matrix(0,length(aa),1)
k<-0
est<-list()

for (j in aa){
  k<-k+1
  est[[k]]<-BJ.NormalCens(cc, X,y, lambda=aa[k], alpha=1, error = 0.000001, iter.max = 1000)
  BICc[k]<-est[[k]]$BIC
} 

p1=data.frame(cbind(aa,BICc))
mlam = aa[BICc==min(BICc)]

p2 <- ggplot(p1, aes(x = aa, y = BICc)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = mlam, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(
    x = expression(lambda),
    y = "BIC"
  ) +
  theme_minimal()



library(patchwork)
p11+p2


# -----------------------------
# 3. Bootstrap for Variable Selection
# -----------------------------



# 1. Split indices by censoring status
censored_data <- wage[1:428, ] # First 428 are censored
uncensored_data <- wage[429:753, ]  # Remaining are uncensored

# 2. Number of bootstrap iterations
B <- 300 # Adjust as needed

# 3. Set up containers
candidate_vars <- colnames(X)[-19]
selected_count <- setNames(numeric(length(candidate_vars)), candidate_vars)

# 4. Stratified bootstrap loop
set.seed(123)

aa=Sys.time()

for (i in seq_len(B)) {
  
  # Bootstrapping censored data
  censored_indices <- sample(1:nrow(censored_data), size = nrow(censored_data), replace = TRUE)
  censored_sample <- censored_data[censored_indices, ]
  
  # Bootstrapping uncensored data
  uncensored_indices <- sample(1:nrow(uncensored_data), size = nrow(uncensored_data), replace = TRUE)
  uncensored_sample <- uncensored_data[uncensored_indices, ]
  

  # c) Create bootstrap sample
  data_boot <- rbind(censored_sample, uncensored_sample)
  

  
  # Re-extract cc, y, X from the bootstrap sample
  cc_boot <- data_boot$cc
  y_boot  <- data_boot$wage
  X_boot  <- data_boot[, -c(1,7,8,22)]  
  
  X_boot[-c(15,17,18)]  <- scale(X_boot[-c(15,17,18)])
  X_boot=as.matrix(X_boot)
  #X_boot  <- scale(X_boot)
  est_boot <- EM.NormalCens(cc_boot, X_boot, y_boot,
                             alpha = 1,
                             error = 1e-5,
                             iter.max = 1000)
   estimates <- est_boot$theta[-19]
  
   

  
  # e) Which variables are "selected"?
  selected_vars <- candidate_vars[abs(estimates) > 0.1]

  # f) Update selection counts
  for (v in selected_vars) {
    selected_count[v] <- selected_count[v] + 1
  }

}



bb=Sys.time()
bb-aa
# 5. Compute selection proportions
selected_prop_EM <- selected_count/ B
selected_prop_EM*100


#################   boots   traping  for BJ    ######################

censored_data <- wage[1:428, ] # First 428 are censored
uncensored_data <- wage[429:753, ]  # Remaining are uncensored


set.seed(123)  # For reproducibility

# Number of bootstrap samples
R <- 3



BJboot=matrix(0,R,18)

for (i in 1:R) {
  # Bootstrapping censored data
  censored_indices <- sample(1:nrow(censored_data), size = nrow(censored_data), replace = TRUE)
  censored_sample <- censored_data[censored_indices, ]
  
  # Bootstrapping uncensored data
  uncensored_indices <- sample(1:nrow(uncensored_data), size = nrow(uncensored_data), replace = TRUE)
  uncensored_sample <- uncensored_data[uncensored_indices, ]
  
  # Combine censored and uncensored data
  boots <- rbind(censored_sample, uncensored_sample)
  y <- boots$wage
  X <- boots[,-c(1,7,8,22)]
  X[-c(15,17,18)]  <- scale(X[-c(15,17,18)])
  X=as.matrix(X)
  cc<- c(rep(0,428),rep(1,325))
  N<-length(cc)
  
  ###BJ
  aa<-seq(0.1,0.3,length.out=40)
  BICc<-matrix(0,length(aa),1)
  k<-0
  est<-list()
  
  for (j in aa){
    k<-k+1
    est[[k]]<-BJ.NormalCens(cc, X,y, lambda=aa[k], alpha=1, error = 0.000001, iter.max = 1000)
    BICc[k]<-est[[k]]$BIC
  } 
  
  mlam = aa[BICc==min(BICc)]
  #select variable
  estBJ=BJ.NormalCens(cc, X,y, lambda=mlam, alpha=1, error = 0.000001, iter.max = 1000)
  
  BJboot[i,]=(abs(estBJ$theta[-19])>0.00001)+0
  
  
  
  
}

apply(BJboot, 2, mean)*100

# -----------------------------
# 4. Prediction results through Cross validation
# -----------------------------

set.seed(123)

est   <- EM.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
X1    <- X[, est$theta[1:ncol(X)] != 0, drop=FALSE]  # keep selected columns
y1    <- est$yhat
dfr   <- data.frame(y1 = y1, X1)
EMres <- train(y1 ~ ., data = dfr, method = "lm", trControl = ctrl)


######BJ
aa<-seq(0.0001,0.3,length.out=50)

BICc<-matrix(0,length(aa),1)
k<-0
est<-list()

for (j in aa){
  k<-k+1
  est[[k]]<-BJ.NormalCens(cc, X,y, lambda=aa[k], alpha=1, error = 0.000001, iter.max = 1000)
  BICc[k]<-est[[k]]$BIC
} 

###minimum lambda

mlam = aa[BICc==min(BICc)]

#select variable

estBJ=BJ.NormalCens(cc, X,y, lambda=mlam, alpha=1, error = 0.000001, iter.max = 1000)

ctrl <- trainControl(method = "cv", number = 20)

x1=X[, estBJ$theta[1:ncol(X)] != 0, drop=FALSE]
y1<- estBJ$yhat
dfr<-data.frame(y1=y1,X=x1)

BJ.res <- train(y1 ~., data = dfr, method = "lm", trControl = ctrl)



c(BJ.res$results$RMSE, BJ.res$results$Rsquared,BJ.res$results$MAE)
c(EMres$results$RMSE, EMres$results$Rsquared, EMres$results$MAE)


###############################################################################
# END OF SCRIPT
###############################################################################