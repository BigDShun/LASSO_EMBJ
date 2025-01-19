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
y <- scale(y)

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
set.seed(123)  # for reproducibility

B <- 300  # Number of bootstrap resamples

candidate_vars <- colnames(X)  # variable names in the model
selected_count <- setNames(numeric(length(candidate_vars)), candidate_vars)



for (i in seq_len(B)) {
  # (a) Resample indices (with replacement)
  boot_idx <- sample(seq_len(N), size = N, replace = TRUE)
  
  # (b) Create bootstrap sample
  data_boot <- BNA[boot_idx, ]
  
  # Preprocess the bootstrap sample
  cc_boot <- data_boot$cc
  y_boot  <- log10(data_boot$Y)  # log10 transform
  y_boot  <- scale(y_boot)       # standardize
  X_boot  <- scale(as.matrix(data_boot[, 6:60]))
  
  # (c) Run selection algorithm

     est_boot <- EM.NormalCens(cc_boot, X_boot, y_boot,
                               alpha = 1,      # Normal distribution
                               error = 1e-6,
                               iter.max = 1000)
     
     
     estimates <- est_boot$theta[-56]

  selected_vars <- candidate_vars[abs(estimates) > 1e-5]
  

  for (v in selected_vars) {
    selected_count[v] <- selected_count[v] + 1
  }
}



# (f) Compute selection proportions
selected_prop <- selected_count / B
cat("\nSelection proportions:\n")
print(selected_prop)

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


# -----------------------------
# 4. Fit/Check Model Convergence 
# -----------------------------

estEM2 <- EM.NormalCens(cc, X, y,
                        alpha    = 1,      # Normal
                        error    = 1e-6,
                        iter.max = 1000)




# Plot log-likelihood history
temp <- estEM2$lkHistory[1:100]
plot(temp, type = 'l',
     ylab = 'Log-likelihood', xlab = 'Iteration',
     main = 'Log-likelihood (first 100 iterations)')
lines(lowess(1:length(temp), temp), col = "red", lwd = 2)

plot(estEM2$lkHistory, type = 'l',
     ylab = 'Log-likelihood', xlab = 'Iteration',
     main = 'Log-likelihood (full iteration)')
lines(lowess(seq_along(estEM2$lkHistory), estEM2$lkHistory), col = "red", lwd = 2)

# Rolling mean
roll_ll <- zoo::rollmean(estEM2$lkHistory, k = 10, fill = NA)
plot(roll_ll, type = "l",
     main = "Rolling Mean of Log-Likelihood (window=10)",
     xlab = "Iteration", ylab = "Rolling LL")
lines(lowess(seq_along(estEM2$lkHistory), estEM2$lkHistory), col = "red", lwd = 2)

# -----------------------------
# 5. Final Model Fit 
# -----------------------------

chosen_cols <- c(4, 8, 41, 42, 49, 50, 53)
x1 <- cbind(1, X[, chosen_cols])  # add intercept

# Fit final censored model with CensReg.SMN
final_model <- CensReg.SMN(cc = cc,
                           x  = x1,
                           y  = as.vector(y),  # ensure y is a vector
                           cens = "left",
                           dist = "Normal",
                           show.envelope = TRUE)

summary(final_model)

