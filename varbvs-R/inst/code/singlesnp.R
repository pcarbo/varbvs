# My single-SNP calculations.
load("../datafiles/Muscle_Skeletal.ACTN3.1Mb.RData")
storage.mode(X) <- "double"

out <- remove.covariate.effects(X,Z,y)
X   <- out$X
y   <- out$y
rm(out)

