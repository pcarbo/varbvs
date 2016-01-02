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
# This part requires the 'grid' and 'ggplot2' packages.
if (!suppressWarnings(require(grid,quietly = TRUE) &&
                      require(ggplot2,quietly = TRUE))) {
  cat(paste("Not showing plots because packages 'grid' and 'ggplot2'",
            "are not available.\n"))
} else {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2,2)))

  # Plot the variational estimate of the residual variance (sigma).
  N        <- data.frame(sigma,apply(w,c("sigma"),sum))
  names(N) <- c("sigma","w")
  p        <- ggplot(data = N,aes(x = sigma,y = w)) +
              geom_bar(fill = I("yellowgreen"),stat = "identity") +
              xlab("sigma") + ylab("posterior") + ylim(0,0.8) +
              opts(panel.background = theme_blank())
            
  print(p,vp = viewport(layout.pos.row = 1,layout.pos.col = 1))

  # Plot the variational estimate of the prior variance (sa).
  N  <- apply(w,c("sa"),sum)
  S1 <- seq(1,length(N),2)
  S2 <- seq(2,length(N),2)
  x  <- (sa[S1]+sa[S2])/2
  N  <- N[S1] + N[S2]
  N  <- data.frame((sa[S1]+sa[S2])/2,N)
  names(N) <- c("sa","w")
  p        <- ggplot(data = N,aes(x =sa,y = w)) +
              geom_bar(fill = I("yellowgreen"),stat = "identity") +
              xlab("sa") + ylab("posterior") + ylim(0,0.8) +
              opts(panel.background = theme_blank())
  print(p,vp = viewport(layout.pos.row = 1,layout.pos.col = 2))

  # Plot the variational estimate of the prior inclusion probability (q).
  N        <- data.frame(log10q,apply(w,c("log10q"),sum))
  names(N) <- c("log10q","w")
  p        <- ggplot(data = N,aes(x = log10q,y = w)) +
              geom_bar(fill = I("yellowgreen"),stat = "identity") +
              xlab("log10q") + ylab("posterior") + ylim(0,0.8) +
              opts(panel.background = theme_blank())
  print(p,vp = viewport(layout.pos.row = 2,layout.pos.col = 1))

  # Plot the ROC curve from the variational method.
  S   <- order(-alpha)
  e   <- as.double(snps$beta != 0)
  tp  <- c(0,cumsum(e[S]),sum(e[S]))
  fp  <- c(0,cumsum(1-e[S]),sum(1-e[S]))
  p   <- ggplot(data = data.frame(tp,fp),aes(x=fp,y=tp)) +
         geom_line(color = I("dodgerblue")) + xlab("no. false +ves") +
         ylab("no. true +ves") + xlim(c(0,20)) + ylim(c(0,20)) +
         opts(panel.background = theme_blank())
  print(p,vp = viewport(layout.pos.row = 2,layout.pos.col = 2))
}
