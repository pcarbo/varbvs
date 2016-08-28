# Map QTLs for phenotypes measured in CFW (Carworth Farms White)
# outbred mice. Phenotypes include muscle weights---EDL and soleus
# muscle---and testis weight measured at sacrifice. Running this
# script with trait = "testis" reproduces the results and figures
# shown in the second example of Carbonetto et al (2016).
library(varbvs)
library(lattice)

# SCRIPT PARAMETERS
# -----------------
# These script parameters specify the candidate prior log-odds
# settings, the prior variance of the coefficients, and which trait to
# analyze. Set trait to "edl", "soleus" or "testis".
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

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
print(summary(fit))

# Show three genome-wide scans: (1) one using the posterior inclusion
# probabilities (PIPs) computed in the BVS analysis of all SNPs; (2)
# one using the p-values computed using GEMMA; and (3) one using the
# PIPs computed from the BVSR model in GEMMA.
trellis.device(height = 5,width = 10)
trellis.par.set(axis.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.7))
w <- c(normalizelogweights(fit$logw))
j <- which(fit$alpha %*% w > 0.5)
r <- gwscan.gemma[[trait]]
r[is.na(r)] <- 0
chromosomes   <- levels(gwscan.bvsr$chr)
xticks        <- rep(0,length(chromosomes))
names(xticks) <- chromosomes
pos           <- 0
for (i in chromosomes) {
  n         <- sum(gwscan.bvsr$chr == i)
  xticks[i] <- pos + n/2
  pos       <- pos + n + 25
}
print(plot(fit,groups = map$chr,vars = j,gap = 1500,cex = 0.6,
           ylab = "varbvs prob.",vars.xyplot.args = list(cex = 0.6)),
      split = c(1,1,1,3),more = TRUE)
print(plot(fit,groups = map$chr,score = r,vars = j,cex = 0.6,gap = 1500,
           draw.threshold = 5.71,ylab = "-log10 p-value",
           vars.xyplot.args = list(cex = 0.6)),
     split = c(1,2,1,3),more = TRUE)
print(xyplot(p1 ~ plot.x,gwscan.bvsr,pch = 20,col = "midnightblue",
             scales = list(x = list(at = xticks,labels = chromosomes),
                           y = list(at = c(0,0.5,1))),
             ylab = "BVSR prob."),
      split = c(1,3,1,3),more = FALSE)
