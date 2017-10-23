# My single-SNP calculations.
load("../datafiles/Muscle_Skeletal.ACTN3.1Mb.RData")
storage.mode(X) <- "double"
n <- nrow(X)
p <- ncol(X)

Z   <- cbind(1,Z)
out <- varbvs:::remove.covariate.effects(X,Z,y)
X   <- out$X
y   <- out$y
rm(out)

i <- 2654
j <- 2656

markers <- setdiff(1:p,c(i,j))
y0 <- y
y  <- c(y - X[,markers] %*% with(fit,alpha[markers,1] * mu[markers,1]))

sigma <- fit$sigma[1]
sa    <- fit$sa[1]
s     <- sa*sigma/(sa*varbvs:::diagsq(X)[c(i,j)] + 1)
mu    <- s*c(y %*% X[,c(i,j)])/sigma
logBF <- log(s/(sa*sigma)) + mu^2/(2*s)

