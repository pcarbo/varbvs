# TO DO: Explain here what this script does.
library(varbvs)
library(lattice)
library(latticeExtra)

# SCRIPT PARAMETERS
# -----------------
# These script parameters specify the candidate prior log-odds
# settings and which trait to analyze.
logodds <- seq(-5,-3,0.25)
trait   <- "soleus"

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
  sa     <- 0.05
} else if (trait == "testis") {
  family <- "gaussian"
  Z      <- pheno[,"sacwt"]
  sa     <- 0.05
} else if (trait == "abnormal.bmd") {
  family <- "binomial"
  Z      <- pheno[,"batch16"]
  sa     <- 1
}
Z <- as.matrix(Z)

# Only analyze samples for which the phenotype and all the covariates
# are observed.
i <- which(apply(cbind(Z,y),1,function (x) sum(is.na(x)) == 0))
y <- y[i]
Z <- as.matrix(Z[i,])
X <- geno[i,]
rm(pheno,geno)

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
cat("FITTING MODEL TO DATA.\n")
fit <- varbvs(X,Z,y,family,sa = sa,logodds = logodds)

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- c(varbvsindep(fit,X,Z,y)$alpha %*% w)

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
varbvsprint(fit)

# Show three genome-wide scans: (1) one using the posterior inclusion
# probabilities (PIPs) computed in the joint analysis of all
# variables; (2) one using the PIPs that ignore correlations between
# the variables.
trellis.device(height = 5,width = 10)
trellis.par.set(axis.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.7))
i <- which(fit$alpha %*% w > 0.5)
r <- gwscan.gemma[[trait]]
r[is.na(r)] <- 0
print(varbvsplot(fit,groups = map$chr,vars = i,gap = 1500,cex = 0.6,
                 ylab = "posterior prob."),
      split = c(1,1,1,3),more = TRUE)
print(varbvsplot(fit,groups = map$chr,score = log10(pip + 0.001),vars = i,
                 cex = 0.6,gap = 1500,ylab = "log10 posterior prob."),
      split = c(1,2,1,3),more = TRUE)
print(varbvsplot(fit,groups = map$chr,score = r,vars = i,cex = 0.6,
                 gap = 1500,score.line = 5.71,ylab = "-log10 p-value"),
     split = c(1,3,1,3),more = TRUE)
