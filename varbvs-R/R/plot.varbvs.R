# Summarize the variable selection results in a single plot.
plot.varbvs <-
    function (x, score, groups, vars = NULL, var.labels, draw.threshold = NA,
              gap = 0,col = "midnightblue", pch = 20, scales = NULL,
              xlab = "", ylab = "",
              abline.args = list(lty = "dotted",col = "orangered"),
              vars.xyplot.args = list(pch = 20,col = "magenta"),
              vars.ltext.args = list(col = "black",pos = 4,cex = 0.5), ...) {

  # Plotting defaults.
  abline.defaults      <- list(lty = "dotted",col = "orangered")
  vars.xyplot.defaults <- list(pch = 20,col = "magenta")
  vars.ltext.defaults  <- list(col = "black",pos = 4,cex = 0.5)
  
  # CHECK INPUTS
  # ------------
  # Check that the first input is an instance of class "varbvs".
  if (!is(x,"varbvs"))
    stop("Input argument x must be an instance of class \"varbvs\".")
  
  # PROCESS OPTIONS
  # ---------------
  # Calculate the posterior inclusion probabilities (PIPs) if a
  # "score" isn't provided as one of the inputs.
  if (missing(score)) {
    w <- c(normalizelogweights(x$logw))
    y <- x$alpha %*% w
  } else
    y <- score
  p <- length(y)

  # Determine the grouping of the variables. By default, all variables
  # are assigned to one group.
  if (missing(groups))
    groups <- rep(1,p)
  group.labels <- unique(groups)

  # Determine the selected variable labels. By default, use the labels
  # stored in the varbvs data structure.
  if (is.character(vars))
    vars <- match(vars,rownames(x$alpha))
  if (missing(var.labels))
    var.labels <- rownames(x$alpha)[vars]
  if (is.null(var.labels))
    var.labels <- rep("",length(vars))
  
  # Get the defaults for abline.args if not provided.
  abline.args <- c(abline.args,abline.defaults)
  abline.args <- abline.args[!duplicated(names(abline.args))]

  # Get the defaults for vars.xyplot.args if not provided.
  vars.xyplot.args <- c(vars.xyplot.args,vars.xyplot.defaults)
  vars.xyplot.args <- vars.xyplot.args[!duplicated(names(vars.xyplot.args))]

  # Get the defaults for vars.ltext.args if not provided.
  vars.ltext.args <- c(vars.ltext.args,vars.ltext.defaults)
  vars.ltext.args <- vars.ltext.args[!duplicated(names(vars.ltext.args))]
  
  # GENERATE GENOME-WIDE SCAN
  # -------------------------
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

  # CREATE XYPLOT
  # -------------
  # Plot the posterior probabilities, or "scores", highlighting and
  # labeling any selected variables.
  out <- xyplot(y ~ x,data.frame(x=x,y=y),pch = pch,col = col,
                scales = c(scales,
                  list(x = list(at = xticks,labels = group.labels))),
                xlab = xlab,ylab = ylab,
                panel = function (x,y,...) {
                  panel.xyplot(x,y,...);
                  if (!is.na(draw.threshold))
                    do.call("panel.abline",
                            c(list(a = draw.threshold,b = 0),abline.args))
                },,...)
  if (!is.null(vars))
    out <- out +
      as.layer(do.call("xyplot",
                       c(list(x = (y ~ x),data = data.frame(x=x,y=y)[vars,],
                              panel = function (x,y,...) {
                                panel.xyplot(x,y,...);
                                do.call("ltext",
                                        c(list(x=x,y=y,labels=var.labels),
                                          vars.ltext.args))
                              }),vars.xyplot.args)))
  return(out)
}
