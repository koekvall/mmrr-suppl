# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("CNTools")
# install.packages("~/GitHub/mmrr-suppl/somatmut/TCGA2STAT_1.2.tar.gz",
#  repos = NULL, type="source")
# library(TCGA2STAT)
# mut.brca <- getTCGA(disease="BRCA", data.type="Mutation", type="somatic", clinical=TRUE)
# gene.brca <- getTCGA(disease="BRCA", data.type="RNASeq2", type="RPKM", clinical=TRUE)

###############################################################################
# Load data
###############################################################################
library(TCGA2STAT)
dat_list <- readRDS("~/GitHub/mmrr-suppl/somatmut/sumatmut_tcga.Rds")
mut.brca <- dat_list[[1]]
gene.brca <- dat_list[[2]]
rm(dat_list)
###############################################################################


###############################################################################
# Construct responses and predictors
###############################################################################
m2 <- OMICSBind(dat1 = mut.brca$merged.dat, dat2 = gene.brca$merged.dat)
mut <- m2$X
keeps <- which(colMeans(as.matrix(mut[,-c(1:3)])) > .05)
mutations <- names(keeps)
geneX <- m2$Y[match(m2$X$bcr, m2$Y$bcr),match(mutations, colnames(m2$Y))]
mutX <- mut[,match(mutations,colnames(mut))]
X <- mut.brca$clinical[match(m2$Y$bcr,rownames(mut.brca$clinical)), c(2, 11)]
X <- as.numeric(as.character(X[,1]))
Y <- cbind(geneX, mutX)
colnames(Y) <- c(paste(mutations, "-GX", sep=""), paste(mutations, "-SM", sep=""))
rms <- unique(which(is.na(X)), which(rowSums(is.na(Y)) > 0))
Xtemp <- X[-rms]
Y <- Y[-rms,]
rms <- which(is.na(rowSums(Y)))
Xtemp <- Xtemp[-rms]
Y <-  Y[-rms,]
n <- dim(Y)[1]
r <- dim(Y)[2]
X <- matrix(0, nrow=n*r, ncol=2*r)
counter <- 1
for(j in 1:n){
  for(k in 1:r){
      X[counter,(1:2) + (2*(k-1)) ] <- c(1, Xtemp[j])
      counter <- counter + 1
  }
}
Y <- as.matrix(Y)
Y[,1:(r/2)] <- log(Y[,1:(r/2)] + 1)
###############################################################################


###############################################################################
# Do data example
###############################################################################
library(mmmr)
library(pheatmap)
library(viridis)
type <- c(rep(1, r/2), rep(2, r/2))
psi <- c(rep(1e-2, r/2), rep(1, r/2))
M <- matrix(NA, r, r)
M[1:(r/2), (r/2 + 1):r] <- 0
M[(r/2 + 1):r, 1:(r/2)] <- 0
diag(M) <- NA
diag(M)[(r/2 + 1):r] <- 1
M # Null hypothesis says binary and continuous are independent

fitNULL_diag <- try(mmrr(Y = Y, X = X, type = type, M = M,
                                 relative = FALSE,
                                 quiet = c(F, T, T, T),
                                 maxit = c(50, 100, 500, 100),
                                 tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                                 psi = psi,
                                 uni_fit = TRUE))

fitNULL_pgd <- mmrr(Y = Y, X = X, type = type, M = M,
                            relative = FALSE,
                            quiet = c(F, F, F, F),
                            maxit = c(50, 100, 500, 100),
                            tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                            psi = psi, 
                            pgd = TRUE,
                            Beta = fitNULL_diag$Beta,
                            Sigma = fitNULL_diag$Sigma,
                            W = fitNULL_diag$W)

# Change restrictions to fit unconstrained model, but Sigma[jj, jj] = 1 for
# Bernoulli responses for identifiability.
M <- matrix(NA, r, r)
diag(M)[(r/2 + 1):r] <- 1


fit_pgd <- mmrr(Y = Y, X = X, type = type, M = M,
                            relative = FALSE,
                            quiet = c(F, F, F, F),
                            maxit = c(50, 100, 500, 100),
                            tol = c(1e-5, 1e-7, 1e-10, 1e-8),
                            psi = psi, 
                            pgd = TRUE,
                            Beta = fitNULL_pgd$Beta,
                            Sigma = fitNULL_pgd$Sigma,
                            W = fitNULL_pgd$W)

out_list <- lrt_approx(fit_null = fitNULL_pgd, fit_full = fit_pgd)

# saveRDS(list("fitNULL_diag" = fitNULL_diag, "fitNULL_pgd" = fitNULL_pgd,
#             "fit_pgd" = fit_pgd, "lrt_list" = out_list),
#        "~/GitHub/mmrr-suppl/somatmut/fitted_list.Rds")

CorEst <- cov2cor(fit_pgd$Sigma)
colnames(Y) <- c(paste(mutations, "-GEx", sep=""),
                 paste(mutations, "-SM", sep=""))

rownames(CorEst) <- colnames(CorEst) <- colnames(Y)
pheatmap(CorEst, treeheight_row = 0, treeheight_col = 0,
         color = magma(100, begin=.1, end=.9))
###############################################################################
