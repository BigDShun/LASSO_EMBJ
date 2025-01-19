###############################################################################
# Title: Censored Regression with Bootstrap Variable Selection
# Author: [Your Name or Institution]
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
library(rslurm)      # Parallelization via Slurm (optional)
library(parallel)    # Base R parallelization
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
X <- scale(X)             # Standardize covariates

# Create censoring indicator vector: 0 = uncensored, 1 = censored
# (Based on the user's code, 428 zeroes, 325 ones)
cc <- c(rep(0, 428), rep(1, 325))
wage$cc <- cc
N <- length(cc)

# -----------------------------
# 2. Fit EM Model and Plot Convergence
# -----------------------------
set.seed(123)  # For reproducibility

# EM fit with alpha = 1 ( Normal distribution ), error tolerance, and max iterations
estEM <- EM.NormalCens(cc = cc, x = X, y = y,
                       alpha = 1,
                       error = 1e-6,
                       iter.max = 1000)

# Extract log-likelihood history from the EM fit
ll_history <- estEM$lkHistory

# Convergence Plots
par(mfrow = c(1, 2))

# Plot of the first 100 iterations
plot(ll_history[1:100], type = 'l',
     xlab = 'Iteration', ylab = 'Log-Likelihood',
     main = 'Log-Likelihood (First 100 iterations)')
lines(lowess(1:100, ll_history[1:100]), col = "red", lwd = 2)

# Full iteration history
plot(ll_history, type = 'l',
     xlab = 'Iteration', ylab = 'Log-Likelihood',
     main = 'Log-Likelihood (Full Iteration)')
lines(lowess(seq_along(ll_history), ll_history), col = "red", lwd = 2)

# Optional: Rolling mean plot

par(mfrow = c(1, 1))
roll_ll <- rollmean(ll_history, k = 10, fill = NA)
plot(roll_ll, type = "l",
     main = "Rolling Mean of Log-Likelihood (window=10)",
     xlab = "Iteration", ylab = "Rolling Mean of LL")
lines(lowess(seq_along(ll_history), ll_history), col = "red", lwd = 2)



# -----------------------------
# 3. Bootstrap for Variable Selection
# -----------------------------
set.seed(123)  # For reproducibility

B <- 30 # Number of bootstrap replications

# Candidate variable names (columns of X)
candidate_vars <- colnames(X)

# Track how many times each variable is "selected"
selected_count <- setNames(numeric(length(candidate_vars)), candidate_vars)

for (i in seq_len(B)) {
  # 1. Bootstrap indices
  boot_idx <- sample(seq_len(N), size = N, replace = TRUE)
  
  # 2. Bootstrap sample
  data_boot <- wage[boot_idx, ]
  cc_boot   <- data_boot$cc
  y_boot    <- data_boot$wage
  
  # Drop wage, occupation, participation, cc columns
  X_boot <- data_boot[, -c(1, 7, 8, 22)]
  X_boot <- scale(X_boot)
  
  # 3. Fit selection algorithm with your chosen function (BJ, EM, or minBIC).
  #    Example: Using the BJ.NormalCens approach with alpha = 1.
  est_boot <- EM.NormalCens(cc_boot, X_boot, y_boot,
                            alpha = 1,
                            error = 1e-6,
                            iter.max = 1000)
  
  # Remove intercept from the selection criterion if needed
  estimates <- est_boot$theta[-19]  # Adjust index as appropriate
  
  # 4. Identify variables considered "selected"
  #    e.g., those with absolute value > 0.00001
  selected_vars <- candidate_vars[abs(estimates) > 1e-5]
  
  # 5. Update selection counts
  for (v in selected_vars) {
    selected_count[v] <- selected_count[v] + 1
  }
}

# Selection proportions
selected_prop <- selected_count / B
print("Selection proportions across the bootstrap:")
selected_prop

# -----------------------------
# 4. Fit Final Model on Chosen Covariates
# -----------------------------
# Example: Suppose we select columns 1, 2, 5, 11, 16 from X
# (Adjust based on your own selection logic or the results from 'selected_prop')
chosen_cols <- c(1, 2, 5, 11, 16)
x1 <- cbind(1, X[, chosen_cols])  # Add intercept column

# Fit final censored regression model:
#   - cc: censoring indicator (0=uncensored, 1=censored)
#   - x1: design matrix with chosen covariates
#   - y:  outcome vector
#   - nu=3, cens="left", dist="Normal"
#   - show.envelope="TRUE" for diagnostic plots
final_model <- CensReg.SMN(cc, x1, y,
                           nu = 3,
                           cens = "left",
                           dist = "Normal",
                           show.envelope = TRUE)

# Print final model results
summary(final_model)

###############################################################################
# END OF SCRIPT
###############################################################################