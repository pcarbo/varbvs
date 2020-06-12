library(MASS)
pkgbuild::compile_dll()
devtools::load_all()

# Script settings.
n  <- 200
p  <- 400
sd <- c(0,   1,    2)
w  <- c(0.9, 0.05, 0.05)
s  <- 0.1
s0 <- 10^seq(-4,0,length.out = 12)
s0[1] <- 0

# Simulate data.
set.seed(1)
S       <- matrix(s,p,p)
diag(S) <- 1
X       <- mvrnorm(n,rep(0,p),S)
X       <- scale(X,center = TRUE,scale = TRUE)
X       <- X/sqrt(n-1)
k       <- sample(length(w),p,replace = TRUE,prob = w)
beta    <- sd[k] * rnorm(p)
y       <- drop(X %*% beta + rnorm(n)/4)

set.seed(1)
fit1 <- varbvsmix(X,NULL,y,4*s0)
plot(beta,rowSums(with(fit1,alpha * mu)),pch = 20)
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")
