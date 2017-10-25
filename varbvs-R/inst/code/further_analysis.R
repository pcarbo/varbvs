sa    <- 1.035   # sa*sigma is prior variance.
sigma <- 0.1381  # Residual variance.

i <- "chr11_66560624_C_T_b38"
j <- "chr11_66561248_T_C_b38"

res <- summary(lm(y ~ X[,i] + Z))$coef[2,]
mu1 <- res[1]
s1  <- res[2]^2
BF1 <- sqrt(s1/(sa*sigma)) * exp(mu1^2/(2*s1))

res <- summary(lm(y ~ X[,j] + Z))$coef[2,]
mu2 <- res[1]
s2  <- res[2]^2
BF2 <- sqrt(s2/(sa*sigma)) * exp(mu2^2/(2*s2))

# My calculations.
out <- varbvsproxybf(X,Z,y,fit,i)

# Results for SNP i.
print(rbind(c(mu1,out$mu[i,1]),
            sqrt(c(s1,out$s[i,1])),
            c(BF1,out$BF[i,1])))

# Results for SNP j.
print(rbind(c(mu2,out$mu[j,1]),
            sqrt(c(s2,out$s[j,1])),
            c(BF2,out$BF[j,1])))

              
