library(ggplot)
pos <- as.numeric(sapply(strsplit(colnames(X),"_"),function (x) x[2])/1e6

ggplot(data.frame(x = pos,y = singlesnp$alpha[,1]),
       aes(x = x,y = y)) +
  geom_point(color = "black") +
  scale_x_continuous(breaks = seq(169,173,0.2)) +
  theme_cowplot(font_size = 10) +
  labs(x = "base-pair position on chromosome 1 (Mb)",y = "single-SNP PIP")

# markers <- which(pos > 171 & pos < 171.4)
# markers <- which(pos > 171.15 & pos < 171.34)
markers <- 1:p
i       <- sort(sample(markers,140))
R       <- cor(X[,i])
                  
library(ggcorrplot)
rownames(R) <- NULL
colnames(R) <- NULL
ggcorrplot(R^2,type = "upper",outline.col = "white")
    
