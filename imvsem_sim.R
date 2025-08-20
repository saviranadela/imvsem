suppressPackageStartupMessages({
  library(mirt)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

set.seed(123)

## -------------------------------
## 0. settings
## -------------------------------
N <- 10000    # examinees
J <- 50       # items
nbins <- 10   # number of theta percentile bins

## item parameter generators (2PL)
gen_items <- function(J){
  a <- rlnorm(J, meanlog=-0.2, sdlog=0.3)  # discrimination
  b <- rnorm(J, 0, 1)                      # difficulty
  list(a=a, b=b)
}

## ICC (2PL)
P_fun <- function(theta, a, b){
  plogis(a * (theta - b))
}

## -------------------------------
## 1. truth
## -------------------------------
theta_true <- rnorm(N)
items <- gen_items(J)
a_true <- items$a
b_true <- items$b

## 2. probability matrix P(theta, a, b)
Pmat <- outer(theta_true, seq_len(J),
              \(th, j) P_fun(th, a_true[j], b_true[j]))

## 3. responses for estimation
X_train <- matrix(rbinom(N*J, 1, as.vector(Pmat)), N, J)
colnames(X_train) <- paste0("i", seq_len(J))

## 4. fit 2PL model, get SEM curve
mod <- mirt(X_train, 1, itemtype="2PL", SE=TRUE, technical=list(NCYCLES=500))

coef_items <- coef(mod, IRTpars=TRUE, simplify=TRUE)$items
a_hat <- coef_items[,"a"]
b_hat <- coef_items[,"b"]

theta_grid <- seq(-3, 3, length.out=121)
Pg <- sapply(seq_len(J), function(j){
  plogis(a_hat[j] * (theta_grid - b_hat[j]))
})
Ij_grid <- sapply(seq_len(J), function(j){
  a_hat[j]^2 * Pg[,j]*(1 - Pg[,j])
})
Itest_grid <- rowSums(Ij_grid)
SEM_grid <- 1 / sqrt(Itest_grid)

## 5. regenerate responses for evaluation
X_test <- matrix(rbinom(N*J, 1, as.vector(Pmat)), N, J)
colnames(X_test) <- paste0("i", seq_len(J))

## predicted probs from fitted model
theta_for_pred <- theta_true   # use true theta for predictions (cleaner)
P_hat <- outer(theta_for_pred, seq_len(J),
               \(th, j) plogis(a_hat[j] * (th - b_hat[j])))

## 6. bin by true theta percentiles
cuts <- quantile(theta_true, probs = seq(0, 1, length.out = nbins + 1))
bin_id <- cut(theta_true, cuts, include.lowest = TRUE, labels = FALSE)
bin_centers <- sapply(seq_len(nbins), function(k){
  mean(theta_true[bin_id == k], na.rm=TRUE)
})

eps <- 1e-12
clamp <- function(p) pmin(pmax(p, eps), 1 - eps)

## -------------------------------
## 7. IMV calculation
## -------------------------------

calc_imv_binbaseline <- function(bin_k){
  idx <- which(bin_id == bin_k)
  Xk  <- X_test[idx, , drop=FALSE]
  Phk <- clamp(P_hat[idx, , drop=FALSE])
  
  ## model log-loss
  L_model <- -mean(Xk*log(Phk) + (1 - Xk)*log(1 - Phk))
  
  ## baseline: per-item prevalence inside bin
  pbar_item <- clamp(colMeans(Xk))
  L_base <- -mean( sweep(Xk, 2, log(pbar_item), `*`) + 
                     sweep(1 - Xk, 2, log(1 - pbar_item), `*`) )
  
  L_base - L_model
}

calc_imv_globalbaseline <- function(bin_k){
  idx <- which(bin_id == bin_k)
  Xk  <- X_test[idx, , drop=FALSE]
  Phk <- clamp(P_hat[idx, , drop=FALSE])
  
  ## model log-loss
  L_model <- -mean(Xk*log(Phk) + (1 - Xk)*log(1 - Phk))
  
  ## baseline: per-item global prevalence
  pbar_item <- clamp(colMeans(X_test))
  L_base <- -mean( sweep(Xk, 2, log(pbar_item), `*`) + 
                     sweep(1 - Xk, 2, log(1 - pbar_item), `*`) )
  
  L_base - L_model
}

imv_binbaseline    <- sapply(seq_len(nbins), calc_imv_binbaseline)
imv_globalbaseline <- sapply(seq_len(nbins), calc_imv_globalbaseline)

## -------------------------------
## 8. summarize SEM for bins
## -------------------------------
grid_bin <- cut(theta_grid, cuts, include.lowest=TRUE, labels=FALSE)
SEM_bin <- sapply(seq_len(nbins), function(k) mean(SEM_grid[grid_bin==k], na.rm=TRUE))
InvSEM_bin <- 1 / SEM_bin

df_bin <- tibble(
  bin = seq_len(nbins),
  theta_center = bin_centers,
  IMV_bin = imv_binbaseline,
  IMV_global = imv_globalbaseline,
  SEM = SEM_bin,
  InvSEM = InvSEM_bin
)

## -------------------------------
## 9. visualization
## -------------------------------

# overlay comparison
df_plot <- df_bin %>%
  mutate(
    IMV_bin_s = rescale(IMV_bin),
    IMV_global_s = rescale(IMV_global),
    InvSEM_s = rescale(InvSEM)
  ) %>%
  pivot_longer(cols = c(IMV_bin_s, IMV_global_s, InvSEM_s),
               names_to = "metric", values_to = "value")

ggplot(df_plot, aes(theta_center, value, color=metric)) +
  geom_line() +
  geom_point() +
  labs(
    title = "IMV vs 1/SEM across ability (2PL)",
    subtitle = "Comparing bin-baseline vs global-baseline IMV",
    x = expression(theta~bin~center),
    y = "scaled value"
  ) +
  theme_minimal(base_size = 12)

# direct relationship plot for global-baseline IMV
ggplot(df_bin, aes(SEM, IMV_global)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "IMV (global baseline) vs SEM per bin",
    x = "SEM in bin",
    y = "IMV (log-loss improvement)"
  ) +
  theme_minimal(base_size = 12)
