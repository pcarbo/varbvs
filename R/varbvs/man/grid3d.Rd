\name{grid3d}
\alias{grid3d}
\title{Generates possible combinations of 3 variables}

\description{It is equivalent to \code{ndgrid(x,y,z)} in MATLAB.}

\usage{grid3d(x, y, z)}

\arguments{
  \item{x,y,z}{vectors specifying the domains for variables X, Y and Z}
}

\value{Returns a list containing three components, \code{X}, \code{Y}
and \code{Z}. Each is an array with one element for every possible
combination of X, Y and Z, so that \code{c(X[i], Y[i], Z[i])} is a
possible assignment to variables X, Y and Z.}

\author{Peter Carbonetto}

\examples{

## Generate all possible combinations of three binary variables.
grid <- grid3d(0:1,0:1,0:1)

}
