# Large-scale Bayesian variable selection for R

[![CRAN status badge](http://www.r-pkg.org/badges/version/varbvs)](https://cran.r-project.org/package=varbvs)
[![Build Status](https://travis-ci.org/pcarbo/varbvs.svg?branch=master)](https://travis-ci.org/pcarbo/varbvs)
[![codecov](https://codecov.io/gh/pcarbo/varbvs/branch/master/graph/badge.svg)](https://codecov.io/gh/pcarbo/varbvs)

### Overview

We introduce *varbvs*, an R package for analysis of large-scale data
sets using Bayesian variable selection methods. To facilitate
application of Bayesian variable selection to a range of problems, the
*varbvs* interface hides most of the complexities of modeling and
optimization, while also providing many options for adaptation to
range of applications. The *varbvs* software has been used to
implement Bayesian variable selection for large problems with over a
million variables and thousands of samples, including analysis of
massive genome-wide data sets.

The R package been tested by
[Travis CI](https://travis-ci.org/pcarbo/varbvs.svg?branch=master),
and the tests' code coverage has been assessed by
[Codecov](https://codecov.io/gh/pcarbo/varbvs).

### Citing varbvs

If you find that this software is useful for your research project,
please cite our paper:

Carbonetto, P., and Stephens, M. (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. *Bayesian Analysis* **7**,
73-108.

### License

Copyright (c) 2012-2016, Peter Carbonetto.

The *varbvs* source code repository by
[Peter Carbonetto](http://github.com/pcarbo) is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *varbvs*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](../LICENSE) for
the full text of the license.

### Installing the package

To install the official release of the *varbvs* package available from
CRAN ([link](http://www.r-pkg.org/pkg/varbvs)), in R simply run:

```R
install.packages("varbvs")
```

Alternatively, you can to install the most up-to-date development
version. The easiest way to accomplish this is using the
[devtools](http://www.r-pkg.org/pkg/devtools) package:

```R
install.packages("devtools")
library(devtools)
install_github("pcarbo/varbvs",subdir = "varbvs-R")
```

Without devtools, it is a little more complicated, but not
hard. Begin by downloading the github repository for this project. The
simplest way to do this is to
[download the repository as a ZIP archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
you have extracted the files from the compressed archive, you will see
that the main directory has two subdirectories, one containing the
MATLAB code, and the other containing the files for the R package.

This subdirectory has all the necessary files to build and install a
package for R. To install this package, follow the
[standard instructions](https://cran.r-project.org/doc/manuals/R-admin.html)
for installing an R package from source. On a Unix or Unix-like
platform (e.g., Mac OS X), the following steps should install the R
package:

```sh
mv varbvs-R varbvs
R CMD build varbvs
R CMD INSTALL varbvs_2.0.9.tar.gz
```

### Using the package

Once you have installed the package, load the package in R by entering

```R
library(varbvs)
```

To get an overview of the package, enter

```R
help(package = "varbvs")
```

The key function in this package is function <code>varbvs</code>.
Here is an example in which we fit the variable selection model to the
Leukemia data:

```R
library(varbvs)
data(leukemia)
fit <- varbvs(leukemia$x,NULL,leukemia$y,family = "binomial",
              logodds = seq(-3.5,-1,0.1),sa = 1)
print(summary(fit))
```

To get more information about this function, type

```R
help(varbvs)
```

### Working examples

We have provided several R scripts in the [vignettes](vignettes) and
[tests](tests) folders to illustrate application of *varbvs* to small
and large data sets:

+ Script [demo.qtl.R](tests/testthat/demo.qtl.R) demonstrates how to
use the varbvs function for mapping a quantitative trait (*i.e.*, a
continuously valued outcome) in a small, simulated data set. Script
[demo.cc.R](tests/testthat/demo.cc.R) demonstrates mapping of a binary
valued outcome in a simulated data set.

+ The [leukemia.R](vignettes/leukemia.R) vignette demonstrates
application of both *glmnet* and *varbvs* to the Leukemia data. The
main aim of this script is to illustrate some of the different
properties of *varbvs* (Bayesian variable selection) and *glmnet*
(penalized sparse regression).

+ Like `demo.qtl.R`, the [cfw.R](vignettes/cfw.R) vignette
demonstrates *varbvs* for mapping genetic factors contributing to a
quantitative trait, but here it is applied to an actual data set
generated from an outbred mouse study.

+ Finally, the [cd.R](vignettes/cd.R) and [cytokine.R](cytokine.R)
illustrate how the *varbvs* package can be applied to a very large
data set to map genetic loci and test biological hypotheses about
genetic factors contributing to human disease risk. Although we cannot
share the data needed to run these scripts due to data privacy
restrictions, we have included these scripts because it is helpful to
be able to follow the steps given in these R scripts.

## How to build static HTML documentation

These are the R commands to build the website (make sure you are
connected to Internet while running these commands, and the working
directory is set to `varbvs-R`):

```R
library(pkgdown)
build_site(mathjax = FALSE)
```

Before building the website, move the `RcppExports.R` file temporarily
outside the `R` directory, otherwise `build_site` will report an
error. Once the pkgdown site is successfully built, move the
`RcppExports.R` file back to its original location, overwriting the
newly generated file.

### Credits

The *varbvs* software package was developed by:<br>
[Peter Carbonetto](http://pcarbo.github.io)<br>
Dept. of Human Genetics, University of Chicago<br>
2012-2017

Xiang Zhou, Xiang Zhu and Matthew Stephens have also contributed to
the development of this software.
