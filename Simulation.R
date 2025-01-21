###############################################################################
# Title: Censored Normal Regression with Various Penalties and Simulations
# Author: [Dashun Liu]
# Description:
#   1) Defines EM-based algorithms for Normal Censored Regression with Lasso,
#      BJ (Buckley-James style) approach, and an Adaptive Lasso approach.
#   2) Demonstrates usage in multiple simulation studies:
#       - Generating data with user-defined censoring rates
#       - Studying performance (TPR, FPR, Precision, Accuracy)
#       - Comparisons of different penalty/regularization approaches
#   3) Illustrates how to apply the models to p>n scenarios.
###############################################################################

rm(list = ls(all = TRUE))

# -----------------------------
# 0. Load Required Libraries
# -----------------------------
library(glmnet)
library(relliptical)   # for mvtelliptical function
library(caret)         # for trainControl, train
library(latex2exp)
library(MASS)

# -----------------------------
# 1. Function Definitions
# -----------------------------

###############################################################################
# EM.NormalCens
###############################################################################
# Performs an EM-based estimation for a left-censored linear model
# assuming Gaussian errors. Uses glmnet with cross-validation in M-step.
#
# Arguments:
#   cc      : censoring indicator (n x 1), 1 => censored, 0 => uncensored
#   x       : design matrix (n x p)
#   y       : response (length n)
#   alpha   : 0 => Ridge, 1 => Lasso, or in (0,1) => elastic net
#   error   : tolerance for EM convergence
#   iter.max: maximum number of EM iterations
#
# Returns:
#   A list containing:
#     theta     = final (beta, sigma^2)
#     iter      = number of EM iterations
#     logver    = final log-likelihood
#     BIC       = approximate BIC measure
#     yhat      = final E[Y_i]
#     lkHistory = vector of log-likelihood values per iteration
###############################################################################
EM.NormalCens <- function(cc, x, y, alpha = 1, error = 1e-5, iter.max = 100){
  p <- ncol(x)
  n <- nrow(x)
  
  # Initial estimate from glmnet (small lambda)
  reg  <- glmnet(x, y, lambda = 0.01, family = "gaussian",
                 intercept = FALSE, alpha = alpha)
  beta <- reg$beta[, 1]
  mu   <- x %*% beta
  sigma2 <- sum((y - mu)^2) / n
  
  lkHistory <- numeric(iter.max)
  
  cont     <- 0
  criterio <- 1
  lkante   <- 1
  
  while (criterio > error && cont < iter.max){
    cont <- cont + 1
    
    # E-step
    E01 <- y
    E02 <- y^2
    
    # For censored observations, compute truncated normal means & second moments
    if (sum(cc) > 0){
      mu_cens <- mu[cc == 1]
      y_cens  <- y[cc == 1]
      np      <- length(mu_cens)
      aux1MomW <- matrix(0, np, 2)
      
      for (j in seq_len(np)){
        A1a <- mvtelliptical(
          lower = -Inf, upper = y_cens[j],
          mu = mu_cens[j], Sigma = sigma2,
          dist = "Normal", n = 1e6
        )
        aux1MomW[j, ] <- c(A1a$EY, A1a$EYY)
      }
      
      E01[cc == 1] <- aux1MomW[, 1]
      E02[cc == 1] <- aux1MomW[, 2]
    }
    
    # M-step: fit with cv.glmnet
    mod.eta <- cv.glmnet(x, E01, family = "gaussian",
                         intercept = FALSE, nfolds = 10, alpha = alpha)
    lambda.atual.eta <- mod.eta$lambda.1se
    
    la.eq2.eta <- glmnet(x, E01, lambda = lambda.atual.eta,
                         family = "gaussian", intercept = FALSE, alpha = alpha)
    
    beta <- as.numeric(la.eq2.eta$beta)
    # Zero out extremely small coefficients
    beta[which(abs(beta) < 1e-8)] <- 0
    
    mu     <- x %*% beta
    sigma2 <- sum(E02 - 2 * E01 * mu + mu^2) / n
    
    # Compute log-likelihood
    auxpdf <- dnorm(y[cc == 0], mu[cc == 0], sqrt(sigma2), log = TRUE)
    auxcdf <- pnorm((y[cc == 1] - mu[cc == 1]) / sqrt(sigma2), log.p = TRUE)
    lk     <- sum(auxpdf) + sum(auxcdf)
    
    lkHistory[cont] <- lk
    
    criterio <- abs(lk / lkante - 1)
    lkante   <- lk
  }
  
  # Trim lkHistory to actual number of iterations
  lkHistory <- lkHistory[1:cont]
  
  # final outputs
  yhat  <- E01
  nf    <- (p - length(which(beta == 0))) + 1  # +1 for sigma^2
  bic   <- (-2 * lk + log(n) * nf) / n
  theta <- c(beta, sigma2)
  
  list(
    theta     = theta,
    iter      = cont,
    logver    = lk,
    BIC       = bic,
    yhat      = yhat,
    lkHistory = lkHistory
  )
}

###############################################################################
# BJ.NormalCens
###############################################################################
# Buckley-James style approach with Normal assumption for left-censored data.
#   This is a simplified approach that replaces censored observations by E(Y).
###############################################################################
BJ.NormalCens <- function(cc, x, y, alpha = 1, error = 1e-4, iter.max = 100){
  p <- ncol(x)
  n <- nrow(x)
  
  # Initial values from glmnet
  reg   <- glmnet(x, y, lambda = 0.01, family = "gaussian",
                  intercept = FALSE, alpha = alpha)
  beta  <- reg$beta[,1]
  mu    <- x %*% beta
  sigma2 <- sum((y - mu)^2) / n
  
  cont     <- 0
  criterio <- 1
  lkante   <- 1
  
  while (criterio > error && cont < iter.max){
    cont <- cont + 1
    
    # E-step
    E01 <- y
    if (sum(cc) > 0){
      mu_cens <- mu[cc == 1]
      y_cens  <- y[cc == 1]
      np      <- length(mu_cens)
      aux1MomW <- numeric(np)
      
      for (j in seq_len(np)){
        A1a <- mvtelliptical(
          lower = -Inf, upper = y_cens[j],
          mu = mu_cens[j], Sigma = sigma2,
          dist = "Normal", n = 1e6
        )
        aux1MomW[j] <- A1a$EY
      }
      E01[cc == 1] <- aux1MomW
    }
    
    # M-step
    mod.eta <- cv.glmnet(x, E01, family = "gaussian",
                         intercept = FALSE, nfolds = 20, alpha = alpha)
    lambda.atual.eta <- mod.eta$lambda.min
    
    la.eq2.eta <- glmnet(x, E01, lambda = lambda.atual.eta,
                         family = "gaussian", intercept = FALSE, alpha = alpha)
    
    beta  <- as.numeric(la.eq2.eta$beta)
    beta[abs(beta) < 1e-8] <- 0
    mu    <- x %*% beta
    sigma2 <- sum((E01 - mu)^2) / n
    
    # Log-likelihood
    auxpdf <- dnorm(y[cc==0], mu[cc==0], sqrt(sigma2), log = TRUE)
    auxcdf <- pnorm((y[cc==1] - mu[cc==1]) / sqrt(sigma2), log.p = TRUE)
    lk     <- sum(auxpdf) + sum(auxcdf)
    
    criterio <- abs((lk / lkante) - 1)
    lkante   <- lk
  }
  
  yhat       <- E01
  nf         <- (p - length(which(beta == 0))) + 1
  bic        <- (-2 * lk + log(n) * nf) / n
  theta      <- c(beta, sigma2)
  
  list(
    theta  = theta,
    iter   = cont,
    logver = lk,
    BIC    = bic,
    yhat   = yhat
  )
}

###############################################################################
# EM.NormalCens_Adaptive
###############################################################################
# Adaptive Lasso with EM for left-censored data. We compute one set of adaptive
# weights based on an initial glmnet fit, then proceed with the EM updates.
###############################################################################
EM.NormalCens_Adaptive <- function(cc, x, y,
                                   alpha = 1, gamma = 1,
                                   error = 1e-5, iter.max = 100,
                                   nfolds = 5){
  p <- ncol(x)
  n <- nrow(x)
  
  # 1) Initial fit to get beta_init
  init_fit  <- glmnet(x, y, lambda = 0.01, family = "gaussian",
                      intercept = FALSE, alpha = alpha)
  beta_init <- as.numeric(init_fit$beta)
  
  # 2) Adaptive weights
  eps <- 1e-6
  w   <- 1 / (abs(beta_init) + eps)^gamma
  
  # 3) Initialize
  mu     <- x %*% beta_init
  sigma2 <- mean((y - mu)^2)
  
  cont     <- 0
  criterio <- 1
  lkante   <- 1
  
  while (criterio > error && cont < iter.max){
    cont <- cont + 1
    
    # E-step
    E01 <- y
    E02 <- y^2
    
    if (any(cc == 1)){
      mu_cens <- mu[cc == 1]
      y_cens  <- y[cc == 1]
      np      <- length(mu_cens)
      aux     <- matrix(0, np, 2)
      
      for (j in seq_len(np)){
        out_j <- mvtelliptical(
          lower = -Inf, upper = y_cens[j],
          mu = mu_cens[j], Sigma = sigma2,
          dist = "Normal", n = 1e5
        )
        aux[j, ] <- c(out_j$EY, out_j$EYY)
      }
      E01[cc == 1] <- aux[,1]
      E02[cc == 1] <- aux[,2]
    }
    
    # M-step: fit with penalty.factor = w
    mod.eta <- cv.glmnet(x, E01, family = "gaussian", intercept = FALSE,
                         alpha = alpha, penalty.factor = w,
                         nfolds = nfolds)
    lambda.sel <- mod.eta$lambda.1se
    
    fit <- glmnet(x, E01, lambda = lambda.sel, family = "gaussian",
                  intercept = FALSE, alpha = alpha, penalty.factor = w)
    beta <- as.numeric(fit$beta)
    beta[abs(beta) < 1e-8] <- 0
    
    mu     <- x %*% beta
    sigma2 <- mean(E02 - 2 * E01 * mu + mu^2)
    
    # Log-likelihood
    auxpdf <- dnorm(y[cc == 0], mean = mu[cc == 0], sd = sqrt(sigma2), log = TRUE)
    auxcdf <- pnorm((y[cc == 1] - mu[cc == 1]) / sqrt(sigma2), log.p = TRUE)
    lk     <- sum(auxpdf) + sum(auxcdf)
    
    criterio <- abs(lk / lkante - 1)
    lkante   <- lk
  }
  
  # Finalize
  yhat       <- E01
  theta      <- c(beta, sigma2)
  n_nonzero  <- sum(abs(beta) > 1e-8)
  nf         <- n_nonzero + 1
  bic        <- (-2 * lk + log(n) * nf) / n
  
  list(
    theta  = theta,
    iter   = cont,
    logver = lk,
    BIC    = bic,
    yhat   = yhat
  )
}










# -----------------------------
# 2. Simulation 1: p < n
# -----------------------------


N <- 1000
p <- 7
censoring_rate <- 0.5
iterations     <- 1000

X <- matrix(rnorm(N * p, 0, 1), ncol = p)
sigma2 <- 1.2
beta   <- c(0, 1, 1, 0, 0, 1, 0)

cutoff_sum <- 0

cat("Finding average cutoff for desired censoring rate...\n")
for (i in seq_len(iterations)){
  y_i <- X %*% beta + rnorm(N, sd = sqrt(sigma2))
  aa  <- sort(y_i, decreasing = FALSE)
  c_i <- aa[ceiling(censoring_rate * N)]
  cutoff_sum <- cutoff_sum + c_i
}
cutoff <- cutoff_sum / iterations

cat("Average cutoff computed =", cutoff, "\n")

# Example usage: apply cutoff to y, set cc=1 if y <= cutoff
set.seed(123)
y <- X %*% beta + rnorm(N, sd = sqrt(sigma2))

selected_indices <- sample(seq_len(N), size = ceiling(censoring_rate * N), replace = FALSE)
y[selected_indices][y[selected_indices] < cutoff] <- cutoff
y[y < cutoff] <- cutoff

cc <- as.numeric(y <= cutoff)

hist(y, main = "Histogram of y after applying cutoff")

# -----------------------------
# 3. Additional Simulations
# -----------------------------
# Example placeholders for subsequent simulations 
beta1 <- c(0, -1, 1, 0, 0, 1, 0)
beta2 <- c(0, -0.2, 0.3, 0, 0, -0.3, 0)
beta3 <- c(-0.3, 0, 0, 0.5, 0, -0.3, 0)

sigma2 <- 1.2
LassoEM <- list()
LassoBJ <- list()

bEM <- matrix(0, 5, 8)
bBJ <- matrix(0, 5, 8)

# Example set of N and censoring levels
set <- expand.grid(c(300, 600, 1000), c(0.1, 0.3, 0.5))

# Example placeholders for storing results
EMres <- cbind(set, matrix(data = rep(0, 63), nrow = 9))
BJres <- cbind(set, matrix(data = rep(0, 63), nrow = 9))

metric <- data.frame(
  N = c(300, 300, 600, 600, 1000, 1000, 300, 300, 600, 600, 1000, 1000, 300, 300, 600, 600, 1000, 1000),
  P = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  type = c("EM", "BJ", "EM", "BJ", "EM", "BJ", "EM", "BJ", "EM", "BJ",
           "EM", "BJ", "EM", "BJ", "EM", "BJ", "EM", "BJ"),
  TPR       = rep(0, 18),
  FPR       = rep(0, 18),
  Precision = rep(0, 18),
  Accuracy  = rep(0, 18)
)

mont <- 1
set.seed(11)


# Example loop over scenarios
for (j in seq_len(9)){
  N    <- set[j, 1]
  perc <- set[j, 2]
  bEM  <- matrix(0, mont, 8)
  bBJ  <- matrix(0, mont, 8)
  LassoEM <- list()
  LassoBJ <- list()
  
  X <- matrix(rnorm(N * 7, 0, 1), ncol = 7)  # as example
  for (k in 1:mont){
    cat("Scenario j=", j, " iteration k=", k, "\n")
    y  <- X %*% beta1 + rnorm(N, sd = sqrt(sigma2))  # normal case
    aa <- sort(y, decreasing = FALSE)
    
    cutof <- aa[ceiling(perc * N)]
    cc    <- as.numeric(y <= cutof)
    y[cc == 1] <- cutof
    
    estEM <- EM.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
    estBJ <- BJ.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
    
    LassoEM[[k]] <- estEM$theta
    LassoBJ[[k]] <- estBJ$theta
    
    bEM[k, ] <- as.numeric(abs(LassoEM[[k]]) > 1e-5)
    bBJ[k, ] <- as.numeric(abs(LassoBJ[[k]]) > 1e-5)
  }
  
  # Summaries, TPR, FPR, etc.
  EMres[j, 3:9] <- apply(bEM, 2, sum)[1:7]
  BJres[j, 3:9] <- apply(bBJ, 2, sum)[1:7]
  
  TP_EM <- sum(EMres[j, c(4,5,8)])
  FP_EM <- sum(EMres[j, c(3,6,7,9)])
  FN_EM <- 3*mont - TP_EM
  TN_EM <- 4*mont - FP_EM
  
  TP_BJ <- sum(BJres[j, c(4,5,8)])
  FP_BJ <- sum(BJres[j, c(3,6,7,9)])
  FN_BJ <- 3*mont - TP_BJ
  TN_BJ <- 4*mont - FP_BJ
  
  metric[2*j - 1, 4:7] <- c(
    TP_EM/(TP_EM + FN_EM),
    FP_EM/(FP_EM + TN_EM),
    TP_EM/(TP_EM + FP_EM),
    (TP_EM + TN_EM)/(TP_EM + TN_EM + FN_EM + FP_EM)
  )
  metric[2*j, 4:7] <- c(
    TP_BJ/(TP_BJ + FN_BJ),
    FP_BJ/(FP_BJ + TN_BJ),
    TP_BJ/(TP_BJ + FP_BJ),
    (TP_BJ + TN_BJ)/(TP_BJ + TN_BJ + FN_BJ + FP_BJ)
  )
}


cat("\nFinal metric table:\n")
print(metric)


###############################################################################
# Simulation II:  p>n scenario, compare EM, BJ, Naive, True using caret train()
###############################################################################

Mont <- 1
LassoEM    <- matrix(0, Mont, 3)
LassoBJ    <- matrix(0, Mont, 3)
LassoNaive <- matrix(0, Mont, 3)
LassoTrue  <- matrix(0, Mont, 3)

ctrl <- trainControl(method = "cv", number = 5)


for (k in 1:Mont){
  cat("p>n scenario, iteration k=", k, "\n")
  
  y <- X %*% beta + rnorm(N, sd = sqrt(sigma2))
  yo <- y
  
  if (perc > 0){
    y[y < cutoff] <- cutoff
    cc <- as.numeric(y <= cutoff)
  } else {
    cc <- rep(0, N)
  }
  
  # 1) EM fit
  est   <- EM.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
  X1    <- X[, est$theta[1:p] > 0, drop=FALSE]  # keep selected columns
  y1    <- est$yhat
  dfr   <- data.frame(y1 = y1, X1)
  EMres <- train(y1 ~ ., data = dfr, method = "lm", trControl = ctrl)
  
  # 2) Naive approach (ignore censoring => cc=0)
  Naive.aux <- EM.NormalCens(rep(0, N), X, y, error = 1e-8, iter.max = 1000)
  X2        <- X[, Naive.aux$theta[1:p] > 0, drop=FALSE]
  dfr1      <- data.frame(y1 = y, X2)
  Naive.res <- train(y1 ~ ., data = dfr1, method = "lm", trControl = ctrl)
  
  # 3) True approach: Suppose we know the true relevant columns
  X3 <- X[, c(2,3,6), drop=FALSE]
  dfr2 <- data.frame(y1 = yo, X3)
  True.res <- train(y1 ~ ., data = dfr2, method = "lm", trControl = ctrl)
  
  # 4) BJ approach
  BJres <- BJ.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
  X4    <- X[, BJres$theta[1:p] > 0, drop=FALSE]
  y2    <- BJres$yhat
  dfr3  <- data.frame(y1 = y2, X4)
  BJ.reg <- train(y1 ~ ., data = dfr3, method = "lm", trControl = ctrl)
  
  # Store results (RMSE, R2, MAE)
  LassoEM[k, ]    <- c(EMres$results$RMSE, EMres$results$Rsquared, EMres$results$MAE)
  LassoBJ[k, ]    <- c(BJ.reg$results$RMSE, BJ.reg$results$Rsquared, BJ.reg$results$MAE)
  LassoNaive[k, ] <- c(Naive.res$results$RMSE, Naive.res$results$Rsquared, Naive.res$results$MAE)
  LassoTrue[k, ]  <- c(True.res$results$RMSE, True.res$results$Rsquared, True.res$results$MAE)
}


# Example: Saving results to text files
# write.table(LassoEM,    file="EM20.txt", sep=" ", row.names=FALSE, col.names=FALSE)
# write.table(LassoBJ,    file="BJ20.txt", sep=" ", row.names=FALSE, col.names=FALSE)
# write.table(LassoNaive, file="Naive20.txt", sep=" ", row.names=FALSE, col.names=FALSE)
# write.table(LassoTrue,  file="True20.txt", sep=" ", row.names=FALSE, col.names=FALSE)

# Example: Reading results and plotting
# naive <- read.table("Naive20 0.5.txt")
# EM    <- read.table("EM20 0.5.txt")
# BJ    <- read.table("BJ20 0.5.txt")
# True  <- read.table("True20 0.5.txt")
# ... create boxplots, etc.

###############################################################################
# Simulation III (p>n, correlation, different alpha, adaptive, etc.)
###############################################################################

alpha.grid <- seq(0, 1, by = 0.25)  # try alpha in 0,0.25,0.5,0.75,1
best.alpha <- NA
best.bic   <- Inf
best.res   <- NULL

for (a in alpha.grid){
  cat("Testing alpha =", a, "\n")
  resA <- EM.NormalCens(cc = cc, x = X, y = y,
                        alpha = a, error = 1e-5, iter.max = 100)
  if (resA$BIC < best.bic){
    best.bic   <- resA$BIC
    best.alpha <- a
    best.res   <- resA
  }
}
cat("Best alpha =", best.alpha, "with BIC=", best.bic, "\n")

coef.best <- best.res$theta  # (beta, sigma^2) of best alpha

###############################################################################
# Another p>n scenario: correlated data
###############################################################################
N <- 100
p <- 300
perc <- 0.5
beta <- c(rep(1, p/10), rep(0, (9*p/10)))  # 10% nonzero

mont <- 2
bEM   <- matrix(0, mont, p+1)
bNet  <- matrix(0, mont, p+1)
bAdap <- matrix(0, mont, p+1)

LassoEM   <- list()
Elastic   <- list()
Adaptive  <- list()

# Two blocks of correlated columns
rho1 <- 0.8
cov_matrix1 <- matrix(rho1, nrow = 15, ncol = 15)
diag(cov_matrix1) <- 1

rho2 <- 0.2
cov_matrix2 <- matrix(rho2, nrow = 15, ncol = 15)
diag(cov_matrix2) <- 1

# Generate data
set.seed(123)
X1 <- mvrnorm(n = N, mu = rep(0, 15), Sigma = cov_matrix1)
X2 <- mvrnorm(n = N, mu = rep(0, 15), Sigma = cov_matrix2)
X3 <- matrix(rnorm(N*(p-30), 0, 1), ncol = p-30)
X  <- as.matrix(cbind(X1, X2, X3))

for (k in 1:mont){
  cat("p>n correlated scenario, iteration k =", k, "\n")
  
  y <- X %*% beta + rnorm(N, sd = sqrt(sigma2))
  
  # Induce censoring
  aa    <- sort(y, decreasing = FALSE)
  cutof <- aa[ceiling(perc*N)]
  cc    <- as.numeric(y <= cutof)
  y[cc == 1] <- cutof
  
  estEM    <- EM.NormalCens(cc, X, y, error = 1e-5, iter.max = 1000)
  estAdapt <- EM.NormalCens_Adaptive(cc, X, y, alpha = 1, nfolds = 5)
  estNET   <- EM.NormalCens(cc, X, y, alpha = 0.75, iter.max = 1000)
  
  LassoEM[[k]]  <- estEM$theta
  Elastic[[k]]  <- estNET$theta
  Adaptive[[k]] <- estAdapt$theta
  
  bEM[k, ]   <- as.numeric(abs(LassoEM[[k]]) > 1e-4)
  bNet[k, ]  <- as.numeric(abs(Elastic[[k]]) > 1e-4)
  bAdap[k, ] <- as.numeric(abs(Adaptive[[k]]) > 1e-4)
}

# Summaries
EMres  <- apply(bEM, 2, sum)[1:p]
Netres <- apply(bNet, 2, sum)[1:p]
adpres <- apply(bAdap, 2, sum)[1:p]

# True positives: first p/10
TP_EM  <- sum(EMres[1:(p/10)])
FP_EM  <- sum(EMres[((p/10)+1):p])
FN_EM  <- (p/10)*mont - TP_EM
TN_EM  <- (9*p/10)*mont - FP_EM

TP_net <- sum(Netres[1:(p/10)])
FP_net <- sum(Netres[((p/10)+1):p])
FN_net <- (p/10)*mont - TP_net
TN_net <- (9*p/10)*mont - FP_net

TP_adp <- sum(adpres[1:(p/10)])
FP_adp <- sum(adpres[((p/10)+1):p])
FN_adp <- (p/10)*mont - TP_adp
TN_adp <- (9*p/10)*mont - FP_adp

EM_metric  <- c(
  TPR       = TP_EM/(TP_EM + FN_EM),
  FPR       = FP_EM/(FP_EM + TN_EM),
  precision = TP_EM/(TP_EM + FP_EM),
  accuracy  = (TP_EM + TN_EM)/(TP_EM + FP_EM + FN_EM + TN_EM)
)
Net_metric <- c(
  TPR       = TP_net/(TP_net + FN_net),
  FPR       = FP_net/(FP_net + TN_net),
  precision = TP_net/(TP_net + FP_net),
  accuracy  = (TP_net + TN_net)/(TP_net + FP_net + FN_net + TN_net)
)
adp_metric <- c(
  TPR       = TP_adp/(TP_adp + FN_adp),
  FPR       = FP_adp/(FP_adp + TN_adp),
  precision = TP_adp/(TP_adp + FP_adp),
  accuracy  = (TP_adp + TN_adp)/(TP_adp + FP_adp + FN_adp + TN_adp)
)

ribing <- rbind(EM_metric, Net_metric, adp_metric)
colnames(ribing) <- c("TPR", "FPR", "Precision", "Accuracy")

cat("\nFinal results for correlated p>n scenario:\n")
print(ribing)


###############################################################################
# End of Script
###############################################################################