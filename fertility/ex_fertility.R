Fertility <- readRDS("~/GitHub/mmrr-suppl/fertility/fertility_stat2data.Rds")
library(mmrr)

Xtemp <- cbind(1, Fertility$Age, Fertility$LowAFC, Fertility$MeanAFC, Fertility$FSH)
p <- dim(Xtemp)[2] - 1
r <- 4

# Create design matrix where each response has its own coefficients
X <- matrix(0, nrow=4*333, ncol=4*5)
counter <- 1
for(kk in 1:333){
  for(jj in 1:r){
     X[counter, ((jj-1)*(p+1) + 1):(jj*(p+1))] <- Xtemp[kk,]
     counter <- counter + 1
  }
}
Y <- cbind(sqrt(Fertility$E2), log(Fertility$TotalGn), Fertility$Embryo, Fertility$Oocytes)
type <- c(1, 1, 3, 3)
psi <- c(1e-2, 1e-2, 1e-1, 1e-1)

# Fit models with diagonal covariance
M <- matrix(0, nrow=r, ncol=r)
diag(M) <- NA
ptm <- proc.time()
fitNULL_diag <- try(mmrr(Y = Y, X = X, type = type, M = M,
                                relative = FALSE,
                                quiet = c(F, T, T, T),
                                maxit = c(50, 100, 500, 100),
                                tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                                psi = psi,
                                uni_fit = TRUE))

fitNULL_trust <- try(mmrr(Y = Y, X = X, type = type, M = M,
                                relative = FALSE,
                                quiet = c(F, T, T, T),
                                maxit = c(50, 100, 500, 100),
                                tol = c(1e-5, 1e-12, 1e-12, 1e-12),
                                psi = psi,
                                pgd = FALSE,
                                Beta = fitNULL_diag$Beta,
                                Sigma = fitNULL_diag$Sigma,
                                W = fitNULL_diag$W))

    
fitNULL_pgd <- mmrr(Y = Y, X = X, type = type, M = M,
                          relative = FALSE,
                          quiet = c(F, T, T, T),
                          maxit = c(50, 100, 500, 100),
                          tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                          psi = psi,
                          pgd = TRUE,
                          Beta = fitNULL_trust$Beta,
                          Sigma = fitNULL_trust$Sigma,
                          W = fitNULL_trust$W)

# Fit models with general covariance
M <- matrix(NA, nrow=r, ncol=r)
fit_diag <- try(mmrr(Y = Y, X = X, type = type, M = M,
                                relative = FALSE,
                                quiet = c(F, T, T, T),
                                maxit = c(50, 100, 500, 100),              
                                tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                                psi = psi,
                                uni_fit = TRUE))

fit_trust <- try(mmrr(Y = Y, X = X, type = type, M = M,
                                relative = FALSE,
                                quiet = c(F, T, T, T),
                                maxit = c(50, 100, 500, 100),
                                tol = c(1e-5, 1e-12, 1e-12, 1e-12),
                                psi = psi,
                                pgd = FALSE,
                                Beta = fit_diag$Beta,
                                Sigma = fit_diag$Sigma,
                                W = fit_diag$W))
    
fit_pgd <- mmrr(Y = Y, X = X, type = type, M = M,
                          relative = FALSE,
                          quiet = c(F, T, T, T),
                          maxit = c(50, 100, 500, 100),
                          tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                          psi = psi,
                          pgd = TRUE,
                          Beta = fit_trust$Beta,
                          Sigma = fit_trust$Sigma,
                          W = fit_trust$W)

# Approximate likelihood ratio test of diagonal covariance matrix, i.e.,
# independent responses.
lrt_approx(fit_null = fitNULL_pgd, fit_full = fit_pgd)

proc.time() - ptm


# Get latent correlation matrix, and correlation matrix of responses given
# predictors evaluated at sample average of predictors.

# Latent covariance matrix:
cov2cor(fit_pgd$Sigma)

# Conditional correlation matrix at average predictors:
X_bar <- colMeans(Xtemp)
X_bar_big <- kronecker(diag(r), t(X_bar))
cov2cor(cov_mmrr(X = X_bar_big, Beta = fit_pgd$Beta, Sigma = fit_pgd$Sigma,
                 psi = psi, type = type, num_nodes = 10))


# Test whether small antra follicle count is important predictor
fit_null_safc <- mmrr(Y = Y, X = X[, -c(seq(3, 20, by = 5))],
                                 type = type, M = M,
                                 relative = FALSE,
                                 quiet = c(F, T, T, T),
                                 maxit = c(50, 100, 500, 100),
                                 tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                                 psi = psi,
                                 pgd = TRUE,
                                 Beta = fit_pgd$Beta[-seq(3, 20, by = 5)],
                                 Sigma = fit_pgd$Sigma,
                                 W = fit_pgd$W)
lrt_approx(fit_null = fit_null_safc, fit_full = fit_pgd)

# Get 95 % confidence region for cov(W_3, W_4) = Sigma_{34} by inverting
# the approximate likelihood ratio test statistic
M <- matrix(NA, r, r)
sig_seq <- seq(0.15, 0.25, length.out = 20)
test_stat <- rep(0, length(sig_seq))
M_null <- M
for(jj in 1:length(sig_seq)){
  M_null[3, 4] <- M_null[4, 3] <- sig_seq[jj]
  fit_null_sig34 <- mmrr(Y = Y, X = X,
                         type = type,
                         M = M_null,
                         relative = FALSE,
                         quiet = c(F, T, T, T),
                         maxit = c(50, 100, 500, 100),
                         tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                         psi = psi,
                         pgd = TRUE,
                         Beta = fit_pgd$Beta,
                         Sigma = fit_pgd$Sigma,
                         W = fit_pgd$W)
  test_stat[jj] <- lrt_approx(fit_null = fit_null_sig34, fit_full = fit_pgd)$stat 
}

# Values in the confidence interval
sig_seq[test_stat <= qchisq(0.95, 1)]

# Approximate corresponding correlations
sig_seq[test_stat <= qchisq(0.95, 1)] / sqrt(prod(diag(fit_pgd$Sigma)[3:4]))


