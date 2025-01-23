###############################################################################
# Title: Censored Regression and Variable Selection on BNA Data
# Author: [Dashun Liu]
# Description:
#   This script demonstrates how to:
#     1. Load and preprocess the BNA data.
#     2. Construct a censoring indicator (left-censored at Y=20).
#     3. Perform a bootstrap variable selection procedure.
#     4. Fit a final model using the selected variables (via CensReg.SMN).
# Requirements:
#   - BNA.txt: Input dataset in working directory.
#   - Packages: relliptical, caret, AER, ggplot2, zoo, etc.
###############################################################################

# -----------------------------
# 0. Load Libraries
# -----------------------------
library(relliptical)  # for EM.NormalCens, BJ.NormalCens
library(caret)        # for data splitting, modeling utilities
library(AER)          # not strictly necessary unless you use the included data
library(ggplot2)      # plotting
library(zoo)          # for rollmean
library(xtable)       # for LaTeX table output


# -----------------------------
# 1. Load and Preprocess Data
# -----------------------------

BNA <- read.csv("BNA.txt", sep = "")
# Response: Y

y <- log10(BNA$Y)


# Covariates: columns 6 to 60
X <- BNA[, 6:60]
X <- as.matrix(X)
X <- scale(X)

# Construct censoring indicator:
#   cc = 1 if Y == 20 (censored), else 0

N  <- nrow(BNA)
cc <- rep(0, N)
cc[which(BNA$Y == 20)] <- 1
BNA$cc <- cc


# -----------------------------
# 2. Bootstrap for Variable Selection
# -----------------------------

# 1. Split indices by censoring status
idx_cc0 <- which(BNA$cc == 0)  # e.g., uncensored
idx_cc1 <- which(BNA$cc == 1)  # e.g., censored

n_cc0 <- length(idx_cc0)  # number of uncensored rows
n_cc1 <- length(idx_cc1)  # number of censored rows

# 2. Number of bootstrap iterations
B <- 30  # Adjust as needed

# 3. Set up containers
candidate_vars <- colnames(X)
selected_count <- setNames(numeric(length(candidate_vars)), candidate_vars)

# 4. Stratified bootstrap loop
set.seed(123)

aa=Sys.time()

for (i in seq_len(B)) {
  
  # a) Sample within each group, maintaining original group sizes
  boot_idx0 <- sample(idx_cc0, size = n_cc0, replace = TRUE)  # uncensored
  boot_idx1 <- sample(idx_cc1, size = n_cc1, replace = TRUE)  # censored
  
  # b) Combine the two groups
  boot_idx  <- c(boot_idx0, boot_idx1)
  
  # c) Create bootstrap sample
  data_boot <- BNA[boot_idx, ]

  # Re-extract cc, y, X from the bootstrap sample
  cc_boot <- data_boot$cc
  y_boot  <- data_boot$Y
  X_boot  <- data_boot[, 6:60]  
  X_boot  <- scale(X_boot)
  

   est_boot <- EM.NormalCens(cc_boot, X_boot, y_boot,
                             alpha = 1,
                             error = 1e-5,
                             iter.max = 1000)
   estimates <- est_boot$theta[-56]
  

  
  # e) Which variables are "selected"?
  selected_vars <- candidate_vars[estimates != 0]
  
  # f) Update selection counts
  for (v in selected_vars) {
    selected_count[v] <- selected_count[v] + 1
  }
}

bb=Sys.time()
bb-aa
# 5. Compute selection proportions
selected_prop <- selected_count / B
selected_prop


# -----------------------------
# 3. Summarize/Save Results
# -----------------------------
BNA_res <- selected_prop
df_BNA  <- as.data.frame(BNA_res)
colnames(df_BNA) <- c("Selection_Prop")

# Example: Keep variables selected > 0.5
df_selected <- df_BNA[df_BNA$Selection_Prop > 0.5, , drop = FALSE]
cat("\nVariables with selection proportion > 0.5:\n")
print(df_selected)


#################   boots   traping  for BJ    ######################

censored_data <- BNA[which(cc==1), ] 
uncensored_data <- BNA[which(cc!=1), ]  


set.seed(123)  # For reproducibility

# Number of bootstrap samples
R <- 30



BJboot=matrix(0,R,55)

time1=Sys.time()
for (i in 1:R) {
  # Bootstrapping censored data
  censored_indices <- sample(1:nrow(censored_data), size = nrow(censored_data), replace = TRUE)
  censored_sample <- censored_data[censored_indices, ]
  
  # Bootstrapping uncensored data
  uncensored_indices <- sample(1:nrow(uncensored_data), size = nrow(uncensored_data), replace = TRUE)
  uncensored_sample <- uncensored_data[uncensored_indices, ]
  
  # Combine censored and uncensored data
  boots <- rbind(censored_sample, uncensored_sample)
  
  y=boots$Y
  y=log10(y)
  #y=scale(y)
  X=boots[,6:60]
  X=as.matrix(X)
  X=scale(X)
  cc=boots$cc

  ###BJ
  aa<-seq(0.1,0.4,length.out=30)
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
  
  BJboot[i,]=(abs(estBJ$theta[-56])>0.00001)+0
  
  
  
  
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


