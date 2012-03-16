*Note: the implementation for R is still under development.*

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
MATLAB for 64-bit Linux. An implementation for R is forthcoming.

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
simplest way to do this is to download the repository as a ZIP
archive. Once you have extracted the files from the compressed
archive, you will see that the main directory contains two
subdirectories, one for the MATLAB functions, and one for the R
code.

Next you will need to compile the C++ code into MATLAB executable
("MEX") files. To do this, you will need to have a [C++ compiler
supported by
MATLAB](http://www.mathworks.com/support/compilers/current_release/),
and you will need to configure MATLAB to build MEX files. See [this
webpage](http://www.mathworks.com/support/tech-notes/1600/1605.html)
for details. When you follow this step, it is important that you
configure MATLAB so that it uses the version of the C++ compiler that
is compatible with your version of MATLAB. Otherwise, you will
encounter errors when building the MEX files, or MATLAB may crash when
attempting to run the examples. If you run into these problems, you may
have to run the mex command with the -v flag to check what compiler is
being used, and you may have to edit the MEX configuration file
manually.

To build the necessary MEX files, run the **install.m** script in
MATLAB. 

*Make note about including flag -DMATLAB_MEX_FILE when building MEX
 files.*

Start by running the script **example1.m**. This script demonstrates
how the variational inference algorithm is used to compute posterior
probabilities for a small linear regression example in which only a
small subset of the variables (single nucleotide polymorphisms, or
SNPs) has affects the outcome (a simulated quantitative trait). In
this small example, the variational estimates of the posterior
probabilities are compared with estimates obtained by MCMC
simulation. Notice that it takes a considerable amount of time to
simulate the Markov chain.

###Quick start for R

*Instructions go here.*

###Overview of MATLAB functions 

The MATLAB subdirectory contains over 40 functions. Here are the most
important ones:

+ **varbvs.m** returns variational estimates of the posterior statistics
for the linear regression model with spike and slab priors, given
choices for the hyperparameters. It computes the posterior statistics
by running the coordinate ascent updates until they converge at a
local minimum of the Kullback-Leibler divergence objective (which
corresponds to a local maximum of the variational lower bound to the
marginal log-likelihood). This function implements the "inner loop" in
the *Bayesian Analysis* paper.

+ **varbvsbin.m** is the same as **varbvs.m**, except that it is meant
  for a logistic regression model instead of a linear regression
  model. This is useful for modeling a binary-valued outcome, such as
  disease status in a case-control study.

+ **varsimbvs.m** demonstrates how to run the full variational
  inference procedure for Bayesian variable selection in linear
  regression, in which we fit both the regression coefficients and the
  hyperparameters to the data. It runs both the "inner" and "outer"
  loops of the inference algorithm, where the inner loop executes the
  coordinate ascent updates for a given value of the hyperparameters,
  and the outer loop runs importance sampling to estimate the
  posterior of the hyperparameters. This is precisely the variational
  inference procedure used in the two simulation studies presented in
  the *Bayesian Analysis* paper. This function assumes specific
  choices for priors on the hyperparameters, as described in the
  paper and in the comments at the top of the file.

###Overview of R functions

*Details about main functions go here.*

###Who

Variational inference for Bayesian variable selection was developed by:<br>
[Peter Carbonetto]((http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
March 2012
