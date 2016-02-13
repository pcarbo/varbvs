# Summarize the variable selection results in a single plot.
plot.varbvs <- function (fit, score, groups, vars, var.labels, gap = 0,
                         col = "midnightblue", var.col = "magenta", pch = 20,
                         ltext.args = "col=\"black\",pos=4,cex=0.5",
                         score.line, xlab = "", ylab = "", ...) {

  # CHECK INPUTS
  # ------------
  # Check that the first input is an instance of class "varbvs".
  if (!is(fit,"varbvs"))
    stop("Input fit must be an instance of class \"varbvs\".")
  
  # PROCESS OPTIONS
  # ---------------
  # Calculate the posterior inclusion probabilities (PIPs) if a
  # "score" isn't provided as one of the inputs.
  if (missing(score)) {
    w <- c(normalizelogweights(fit$logw))
    y <- fit$alpha %*% w
  } else
    y <- score
  p <- length(y)

  # Determine the grouping of the variables. By default, all variables
  # are assigned to one group.
  if (missing(groups))
    groups <- rep(1,p)
  group.labels <- unique(groups)

  # Determine the selected variable labels. By default, use the labels
  # stored in the varbvs data structure ("fit").
  if (missing(var.labels))
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
                panel = function (x,y,...) {
                  panel.xyplot(x,y,...);
                  if (!missing(score.line))
                    panel.abline(a = score.line,b = 0,lty = "dotted",
                                 col = "orangered")
                },
                xlab = xlab,ylab = ylab,...) +
         as.layer(xyplot(y ~ x,data.frame(x = x,y = y)[vars,],pch = pch,
           col = var.col,
           panel = function (x,y,...) {
             panel.xyplot(x,y,...);
             eval(parse(text = paste("ltext(x=x,y=y,labels=var.labels,",
                          ltext.args,")")))
           })))
}
