# TO DO: Explain here what this script does.
library(varbvs)
library(lattice)
library(latticeExtra)

# SCRIPT PARAMETERS
# -----------------
# These script parameters specify the candidate prior log-odds
# settings and which trait to analyze.
logodds <- seq(-5,-3,0.25)
trait   <- "testis"

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPE AND PHENOTYPE DATA
# --------------------------------
cat("LOADING DATA.\n")
load("cfw.RData")
y <- pheno[,trait]
if (trait == "edl" | trait == "soleus") {
  family <- "gaussian"
  Z      <- pheno[,c("batch16","tibia")]
} else if (trait == "testis") {
  family <- "gaussian"
  Z      <- pheno[,"sacwt"]
} else if (trait == "abnormal.bmd") {
  family <- "binary"
  Z      <- pheno[,"batch16"]
}
Z <- as.matrix(Z)

# Only analyze samples for which the phenotype and all the covariates
# are observed.
i <- which(apply(cbind(Z,y),1,function (x) sum(is.na(x)) == 0))
y <- y[i]
Z <- Z[i,]
X <- geno[i,]
rm(pheno,geno)

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
cat("FITTING MODEL TO DATA.\n")
fit <- varbvs(X,Z,y,family,logodds = logodds)

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- c(varbvsindep(fit,X,Z,y)$alpha %*% w)

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
varbvsprint(fit)

# Show three "genome-wide scans".
# TO DO.
