# A more detailed look at the FMO2 locus.
library(glmnet)
library(varbvs)

set.seed(1)

singlesnp <- function (X, Z, y) {
 p <- ncol(X)
 if (!is.null(Z))
   y <- resid(lm(y ~ Z))
 out <- matrix(0,p,2)
 for(i in 1:p) {
   g       <- summary(lm(y ~ X[,i]))
   out[i,] <- coef(g)[2,1:2]
 }
 return(list(betahat = out[,1], sebetahat = out[,2],residuals = y))
}

autoselect.mixsd = function(betahat,sebetahat,mult = sqrt(2)){
  # To avoid exact measure causing (usually by mistake)
        sebetahat = sebetahat[sebetahat!=0] 
        # so that the minimum is small compared with measurement precision
        sigmaamin = min(sebetahat)/10 
        if (all(betahat^2 <= sebetahat^2)) {
            # to deal with the occassional odd case where this could happen; 8 is arbitrary
            sigmaamax = 8*sigmaamin 
        } else {
            # this computes a rough largest value you'd want to use, 
            # based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
            sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) 
        }
        if(mult==0){
            return(c(0,sigmaamax/2))
        } else {
            npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
            return(mult^((-npoint):0) * sigmaamax)
        }
    }

load("../datafiles/Thyroid.FMO2.1Mb.RData")
X <- X[,3501:4500]

storage.mode(X) <- "double"
n <- nrow(X)
p <- ncol(X)

out.singlesnp <- singlesnp(X,Z,y)
mixsd <- with(out.singlesnp,autoselect.mixsd(betahat,sebetahat))

set.seed(1)
y0            <- resid(lm(y ~ Z))
out.cv.glmnet <- cv.glmnet(X,y0,nfolds = 20)
lambda.opt <- out.cv.glmnet$lambda.1se

fit.glmnet <- glmnet(X,y0,lambda = lambda.opt,intercept = TRUE,
                     standardize = FALSE)
mu0 <- out.singlesnp$betahat
alpha0 <- rep(0,p)
alpha0[c(115,341)] <- 1
fit.varbvsmix <- varbvs::varbvsmix(X,Z,y,sa = c(0,mixsd^2),verbose = FALSE)
fit.varbvs <- varbvs(X,Z,y,verbose = FALSE)

markers1 <- c("chr1_171172098_C_T_b38",
             "chr1_171199984_T_C_b38",
             "chr1_171122735_A_G_b38",
             "chr1_171133158_A_G_b38")

markers2 <- c("chr1_171172098_C_T_b38",
              "chr1_171122735_A_G_b38")

y0 <- resid(lm(y ~ Z))
summary(lm(y0 ~ X[,markers]))$r.squared


markers3 <- c("chr1_171168633_C_A_b38",
              "chr1_171147265_C_A_b38",
              "chr1_171164750_C_A_b38",
              "chr1_171178589_C_T_b38")
summary(lm(y0 ~ X[,markers]))$r.squared
