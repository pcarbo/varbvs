# Summarize the variable selection results in a single plot.
varbvsplot <- function (fit, score = NULL, vars = NULL, groups = NULL,
                        gap = 0, col = "midnightblue",col.vars = "magenta",
                        xlab = "", ylab = "",...) {
  
  # PROCESS OPTIONS
  # ---------------
  # Calculate the posterior inclusion probabilities (PIPs) if a
  # "score" isn't provided as one of the inputs.
  if (is.null(score)) {
    w <- c(normalizelogweights(fit$logw))
    y <- fit$alpha %*% w
  } else
    y <- score
  p <- length(y)

  # Determine the grouping of the variables. By default, all the
  # variables are assigned to a single group.
  if (is.null(groups))
    groups <- rep(1,p)
  group.labels <- unique(groups)

  # GENERATE GENOME-WIDE SCAN PLOT
  # ------------------------------
  # Determine the positions of the variables along the horizontal axis.
  x      <- rep(0,p)
  pos    <- 0
  xticks <- NULL
  for (i in group.labels) {
    j      <- which(groups == i)
    m      <- length(j)
    x[j]   <- pos + 1:m
    xticks <- c(xticks,pos+m/2)
    pos    <- pos + m + gap
  }

  # Plot the posterior probabilities, highlighting and labeling the
  # selected variables.
  labels <- rownames(fit$alpha)
  return(xyplot(y ~ x,data.frame(x = x,y = y),pch = 20,col = col,
                scales = list(x = list(at = xticks,labels = group.labels)),
                xlab = xlab,ylab = ylab,...) +
         as.layer(xyplot(y ~ x,data.frame(x = x,y = y)[vars,],pch = 20,
                         col = col.vars,
                         panel = function (x,y,...) {
                           panel.xyplot(x,y,...);
                           ltext(x = x,y = y,labels = labels[vars],
                                 pos = 4,cex = 0.5)
                         })))
}
