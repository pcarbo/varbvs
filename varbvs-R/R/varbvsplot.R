# Summarize the variable selection results in a single plot.
varbvsplot <- function (fit, score = NULL, vars = NULL, groups = NULL,
                        gaps = 0) {

  # PROCESS OPTIONS
  # ---------------
  # Calculate the posterior inclusion probabilities (PIPs) if a
  # "score" isn't provided as one of the inputs.
  if (is.null(score))
    y <- score
  else {
    w <- c(normalizelogweights(fit$logw))
    y <- fit$alpha %*% w
  }
  p <- numel(y);

  # Determine the grouping of the variables. By default, all the
  # variables are assigned to a single group.
  if (is.null(groups))
    groups <- rep(1,p)
  group_labels <- unique(groups)

  # GENERATE GENOME-WIDE SCAN PLOT
  # ------------------------------
  # Determine the positions of the variables along the horizontal axis.
  x      <- rep(0,p)
  pos    <- 0
  xticks <- NULL
  for (i in group_labels) {
    j      <- which(groups == i)
    m      <- length(j)
    x[j]   <- pos + 1:m
    xticks <- c(xticks,pos+m/2)
    pos    <- pos + m + gap
  }

  # Plot the posterior probabilities, highlighting and labeling the
  # selected variables.
  return(xyplot(y ~ x,data.frame(x = x,y = y),pch = 20))
}
