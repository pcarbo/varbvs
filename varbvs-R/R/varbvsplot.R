# Summarize the variable selection results in a single plot.
varbvsplot <- function (fit, score = NULL, groups = NULL, gap = 0,
                        col = "midnightblue", vars = NULL, var.labels = NULL,
                        var.col = "magenta", pch = 20, xlab = "", ylab = "",
                        ltext.args = "col=\"black\",pos=4,cex=0.5",...) {
  
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

  # Determine the selected variable labels. By default, use the labels
  # stored in the varbvs data structure ("fit").
  if (is.null(var.labels))
    var.labels <- rownames(fit$alpha)[vars]
  
  # GENERATE GENOME-WIDE SCAN PLOT
  # ------------------------------
  # Determine the positions of the variables and, if necessary, group
  # labels along the horizontal axis.
  if (length(group.labels) == 1) {
    x            <- 1:p
    xticks       <- NULL
    group.labels <- NULL
  } else {
    x      <- rep(0,p)
    pos    <- 0
    xticks <- NULL
    for (i in group.labels) {
      j      <- which(groups == i)
      n      <- length(j)
      x[j]   <- pos + 1:n
      xticks <- c(xticks,pos+n/2)
      pos    <- pos + n + gap
    }
  }
  
  # Plot the posterior probabilities, or "scores", highlighting and
  # labeling selected variables.
  return(xyplot(y ~ x,data.frame(x = x,y = y),pch = pch,col = col,
                scales = list(x = list(at = xticks,labels = group.labels)),
                xlab = xlab,ylab = ylab,...) +
         as.layer(xyplot(y ~ x,data.frame(x = x,y = y)[vars,],pch = pch,
           col = var.col,
           panel = function (x,y,...) {
             panel.xyplot(x,y,...);
             eval(parse(text = paste("ltext(x=x,y=y,labels=var.labels,",
                          ltext.args,")")))
           })))
}
