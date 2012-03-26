# In this small example, we explore the posterior distribution of the
# coefficients in a linear regression, with spike and slab priors. In this
# idealized case, the variables are independent.

# SCRIPT PARAMETERS.
p  <- 1e3  # The number of variables (SNPs).
n  <- 500  # The number of samples.
na <- 20   # Number of variables that effect the outcome ("causal" SNPs).

# The level of residual noise is set so that the proportion of variance
# explained should be a little less than 30% on average (according to my
# rough calculations).
se <- 9

# These two parameters specify the beta prior on the proportion of variables
# (SNPs) that are included in the linear model of Y.
a <- 0.02
b <- 1

# This parameter specifies the prior on the variance of the regression
# coefficients ('sa'). For more information on this prior, see
# 'varsimbvs'.
ca <- 0.02

# Candidate values of the variance of the residual (sigma), the prior
# variance of the regression coefficients ('sa'), and logarithm of the
# prior inclusion probability ('log10q').
sigma  <- seq(8,13)
sa     <- seq(0.025,0.4,0.025)
log10q <- seq(-2.5,-1,0.25)

# Set the random number generator seed.
set.seed(1)

# CREATE THE DATA.
# Note that X and y are centered.
cat("Creating data.\n")
snps <- create.snps(p,na)
data <- create.data(snps$maf,snps$beta,se,n)

# Calculate the proportion of variance explained. Here, SZ is the sample
# genetic variance.
sz <- var(c(data$X %*% snps$beta))
cat(sprintf("Proportion of variance explained is %0.3f.\n",sz/(sz + se)))

# Generate the importance samples at regular intervals in a predefined
# grid.
grid                  <- grid3d(sigma,sa,log10q)
names(grid)           <- c("sigma","sa","log10q")
dimnames(grid$sigma)  <- list(sigma=NULL,sa=NULL,log10q=NULL)
dimnames(grid$sa)     <- dimnames(grid$sa)
dimnames(grid$log10q) <- dimnames(grid$sa)

# COMPUTE VARIATIONAL ESTIMATES.
cat("Computing variational estimates.\n")
result <- varsimbvs(data$X,data$y,grid$sigma,grid$sa,grid$log10q,a,b,ca)
w      <- result$w
alpha  <- result$alpha
mu     <- result$mu
remove(result)

# SHOW POSTERIOR MEAN OF HYPERPARAMETERS.
cat("Approximate posterior means of hyperparameters:\n")
cat("Approximate posterior means of hyperparameters:\n")
cat(sprintf("log10(sigma) %6.3f\n",sum(w*log10(grid$sigma))))
cat(sprintf("log10(sa)    %6.3f\n",sum(w*log10(grid$sa))))
cat(sprintf("log10(q)     %6.3f\n",sum(w*grid$log10q)))

# DISPLAY POSTERIOR DISTRIBUTIONS OF HYPERPARAMETERS.
# Plot the variational estimate of the residual variance (sigma).
barplot(apply(w,1,sum),space=0,border=NA,names.arg=sigma,ylim=c(0,0.8),
        xlab="sigma",ylab="posterior",col="yellowgreen",axes=FALSE,
        cex.names=0.9,cex.axis=0.9)
axis(2,at=c(0,0.4,0.8),line=NA,las=1)
## subplot(3,3,1);
## bar(x,N,1,'LineStyle','none','FaceColor',rgb('yellowgreen'));

## Plot the variational estimate of the prior variance (sa).
barplot(apply(w,2,sum),space=0,border=NA,names.arg=sa,ylim=c(0,0.8),
        xlab="sa",ylab="posterior",col="yellowgreen",axes=FALSE)
axis(2,at=c(0,0.4,0.8),line=NA,las=1)
