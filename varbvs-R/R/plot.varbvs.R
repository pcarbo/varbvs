# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2018, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# Summarize the variable selection results in a single plot.
plot.varbvs <-
    function (x, score, groups, vars = NULL, var.labels, draw.threshold = NA,
              gap = 0,col = "midnightblue", pch = 20, scales = NULL,
              xlab = "variable", ylab = "posterior probability", 
              main = "fitted varbvs model: variable selection results",
              abline.args = list(lty = "dotted",col = "orangered"),
              vars.xyplot.args = list(pch = 20,col = "magenta"),
              vars.ltext.args = list(col = "black",pos = 4,cex = 0.5),
              par.settings = list(par.main.text = list(font = 1,cex = 0.8),
                                  layout.heights = list(axis.top = 0,
                                                        axis.bottom = 0)),
              ...) {

  # CHECK INPUTS
  # ------------
  # Check that the first input is an instance of class "varbvs".
  if (!is(x,"varbvs"))
    stop("Input argument x must be an instance of class \"varbvs\".")
  
  # PROCESS OPTIONS
  # ---------------
  # Get the posterior inclusion probabilities (PIPs) if a
  # "score" isn't provided as one of the inputs.
  if (missing(score)) {
    y <- x$pip
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
                xlab = xlab,ylab = ylab,main = main,
                par.settings = par.settings,
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
