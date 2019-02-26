# Output the roots x to the quadratic equation a*x^2 + b*x + c = 0.
roots2 <- function (a, b, c) {
  y <- -(b + sign(b)*sqrt(b^2 - 4*a*c))/2
  return(c(y/a,c/y))
}
