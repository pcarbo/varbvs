# Map QTLs for phenotypes measured in CFW (Carworth Farms White)
# outbred mice. Phenotypes include muscle weights (EDL and soleus
# muscle) and testis weight measured at sacrifice. Running this script
# with trait = "testis" reproduces the results and figures shown in
# the second example of Carbonetto et al (2016).
library(varbvs)
library(lattice)

# SCRIPT PARAMETERS
# -----------------
# These script parameters specify the candidate prior log-odds
# settings, the prior variance of the coefficients, and which trait to
# analyze.
#
# Set trait to "edl", "soleus" or "testis".
trait   <- "testis"
logodds <- seq(-5,-3,0.25)
sa      <- 0.05

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPE AND PHENOTYPE DATA
# --------------------------------
cat("LOADING DATA.\n")
load("cfw.RData")
y <- pheno[,trait]
if (trait == "edl" | trait == "soleus") {
  Z <- pheno[,c("batch16","tibia")]
} else if (trait == "testis") {
  Z <- pheno[,"sacwt"]
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
runtime <- system.time(fit <- varbvs(X,Z,y,sa = sa,logodds = logodds))
cat(sprintf("Modeling fitting took %0.2f minutes.\n",runtime["elapsed"]/60))

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- c(varbvsindep(fit,X,Z,y)$alpha %*% w)

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
print(summary(fit))

# Show three genome-wide scans: (1) one using the posterior inclusion
# probabilities (PIPs) computed in the joint analysis of all
# variables; (2) one using the PIPs that ignore correlations between
# the variables; and (3) and one using the p-values computed using GEMMA.
trellis.device(height = 5,width = 10)
trellis.par.set(axis.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.7))
i <- which(fit$alpha %*% w > 0.5)
r <- gwscan.gemma[[trait]]
r[is.na(r)] <- 0
print(plot(fit,groups = map$chr,vars = i,gap = 1500,cex = 0.6,
           ylab = "posterior prob.",vars.xyplot.args = list(cex = 0.6)),
      split = c(1,1,1,3),more = TRUE)
print(plot(fit,groups = map$chr,score = log10(pip + 0.001),vars = i,
           cex = 0.6,gap = 1500,ylab = "log10 posterior prob.",
           vars.xyplot.args = list(cex = 0.6)),
      split = c(1,2,1,3),more = TRUE)
print(plot(fit,groups = map$chr,score = r,vars = i,cex = 0.6,gap = 1500,
           draw.threshold = 5.71,ylab = "-log10 p-value",
           vars.xyplot.args = list(cex = 0.6)),
     split = c(1,3,1,3),more = FALSE)
