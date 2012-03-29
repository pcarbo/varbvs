#Variational inference for Bayesian variable selection implemented in MATLAB and R

###Introduction

This software package includes
[MATLAB](http://www.mathworks.com/products/matlab/) and
[R](http://www.r-project.org) implementations of the variational
inference procedure for Bayesian variable selection, as described in
the [Bayesian Analysis](http://ba.stat.cmu.edu/) paper [Scalable
variational inference for Bayesian variable selection in regression,
and its accuracy in genetic association
studies](http://ba.stat.cmu.edu/journal/2012/vol07/issue01/carbonetto.pdf)
(*Bayesian Analysis* 7, March 2012, pages 73-108). This software has
been used to implement Bayesian variable selection for large problems
with over a million variables and thousands of samples.

The MATLAB implementation has been tested in version 7.10 (R2010a) of
MATLAB for 64-bit Linux. The R implementation has been tested in
version 2.14.2 of R for 64-bit Linux.

###License

Variational inference for Bayesian variable selection by [Peter
Carbonetto](http://www.cs.ubc.ca/spider/pcarbo) is free software: you
can redistribute it under the terms of the [GNU General Public
License](http://www.gnu.org/licenses/gpl.html). All the files in this
project are part of Variational inference for Bayesian variable
selection. This project is distributed in the hope that it will be
useful, but **without any warranty**; without even the implied
warranty of **merchantability or fitness for a particular
purpose**. See file LICENSE for the full text of the license.

Copyright 2012 Peter Carbonetto.

This project includes several MATLAB functions created by James
P. LeSage as part of the [Econometrics
Toolbox](http://www.spatial-econometrics.com/) (March 2010 revision),
which is also distributed under the GNU General Public License. None
of code from this toolbox has been modified. Also included (without
any changes) is the MATLAB function rgb.m written by [Kristján
Jónasson](http://www.hi.is/~jonasson).

###Quick start for MATLAB

Start by downloading the github repository for this project. The
simplest way to do this is to [download the repository as a ZIP
archive](http://github.com/pcarbo/Variational-inference-for-Bayesian-variable-selection/zipball/master). Once
you have extracted the files from the compressed archive, you will see
that the main directory contains two subdirectories, one for the
MATLAB functions, and one for the R code.

Next you will need to compile the C code into MATLAB executable
("MEX") files. To do this, you will need to have a [C compiler
supported by
MATLAB](http://www.mathworks.com/support/compilers/current_release/),
and you will need to configure MATLAB to build MEX files. See [this
webpage](http://www.mathworks.com/support/tech-notes/1600/1605.html)
for details. When you follow this step, it is important that you
configure MATLAB so that it uses the version of the C compiler that is
compatible with your version of MATLAB. Otherwise, you will encounter
errors when building the MEX files, or MATLAB may crash when
attempting to run the examples. If you run into these problems, you
may have to run the mex command with the -v flag to check what
compiler is being used, and you may have to edit the MEX configuration
file manually.

To build the necessary MEX files, run the **install.m** script in
MATLAB. 

Note that the beginning of this script sets some compiler and linker
flags. These flags tell the GCC compiler to use the ISO C99 standard,
and to optimize the code as much as possible. However, these flags may
not be relevant to your setup, especially if you are not using
[GCC](gcc.gnu.org). To avoid errors during installation, if you are
using a compiler other than GCC, it may be best to set variables
**cflags** and **ldflags** to empty strings before running the
**install.m** script.

If you modify the installation procedure to fit your compiler setup,
it is important that you define macro MATLAB_MEX_FILE, akin to the
\#define MATLAB_MEX_FILE directive in C. In GCC, this is accomplished
by including flag -DMATLAB_MEX_FILE when issuing the commands to build
MEX files.

Once you have built the MEX files, start by running the script
**example1.m**. This script demonstrates how the variational inference
algorithm is used to compute posterior probabilities for a small
linear regression example in which only a small subset of the
variables (single nucleotide polymorphisms, or SNPs) has affects the
outcome (a simulated quantitative trait). In this small example, the
variational estimates of the posterior probabilities are compared with
estimates obtained by MCMC simulation. Notice that it takes a
considerable amount of time to simulate the Markov chain.

###Quick start for R

Start by downloading the github repository for this project. The
simplest way to do this is to [download the repository as a ZIP
archive](http://github.com/pcarbo/Variational-inference-for-Bayesian-variable-selection/zipball/master). Once
you have extracted the files from the compressed archive, you will see
that the main directory has two subdirectories, one containing the
MATLAB code, and the other containing R files.

The subdirectory **R/varbvs** has all the necessary files to build and
install a package for R. To install this package, follow the [standard
instructions](http://cran.r-project.org/doc/manuals/R-admin.html) for
installing an R package from source. On a Unix or Unix-like platform
(such as Mac OS X), the installation command is
**R CMD INSTALL homedir/varbvs/R/varbvs**.

Once you have installed the package, you can load the functions in R
by running **library(varbvs)**. To get an overview of the package, run
**help(varbvs)** and **help(package="varbvs")** in R.

Once you have installed and loaded the **varbvs** package, start by
running the demonstration script with **demo(example1)**. This script
demonstrates how the variational inference algorithm is used to
compute posterior probabilities for a small linear regression example
in which only a small subset of the variables (single nucleotide
polymorphisms, or SNPs) has affects the outcome (a simulated
quantitative trait). 

Note that the demonstration script **example1.R** requires packages
[ggplot2](http://had.co.nz/ggplot2) and **grid** to show the plots
depicting the posterior distributions of the hyperparameters computed
using the variational inference method. The grid package is normally
included in the base distribution of R, and ggplot2 is available on
[CRAN](http://cran.r-project.org), and can be installed using the
**install.packages** function (for more recent versions of R).

###Overview of R and MATLAB functions 

The MATLAB subdirectory and R package contain several functions. Here
are the most interesting ones:

+ **varbvs** in MATLAB or **varbvsoptimize** in R returns
variational estimates of the posterior statistics for the linear
regression model with spike and slab priors, given choices for the
hyperparameters. It computes the posterior statistics by running the
coordinate ascent updates until they converge at a local minimum of
the Kullback-Leibler divergence objective (which corresponds to a
local maximum of the variational lower bound to the marginal
log-likelihood). This function implements the "inner loop" in the
*Bayesian Analysis* paper.

+ **varbvsbin** in MATLAB or **varbvsbinoptimize** in R is analogous
to **varbvs** or **varbvsoptimize**, except that it is meant for
logistic regression instead of linear regression. This is useful for
modeling a binary-valued outcome such as disease status in a
case-control study.

+ **varsimbvs** (MATLAB and R) demonstrates how to run the full
  variational inference procedure for Bayesian variable selection in
  linear regression, in which we fit both the regression coefficients
  and the hyperparameters to the data. It runs both the "inner" and
  "outer" loops of the inference algorithm, where the inner loop
  executes the coordinate ascent updates for a given value of the
  hyperparameters, and the outer loop runs importance sampling to
  estimate the posterior of the hyperparameters. This is precisely the
  variational inference procedure used in the two simulation studies
  presented in the *Bayesian Analysis* paper. This function assumes
  specific choices for priors on the hyperparameters, as described in
  the paper and in the comments at the top of the file.

###Who

Variational inference for Bayesian variable selection was developed by:<br>
[Peter Carbonetto]((http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
March 2012
