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
of MATLAB code from this toolbox has been modified. Also included
(without modification) is the MATLAB function rgb.m written by
Kristján Jónasson.

###Quick start for MATLAB

Start by downloading the github repository for this project. The
simplest way is to download the repository as a ZIP archive. Once you
have extracted the files from the compressed archive, you will see
that the project directory contains two subdirectories for the MATLAB
and R code.

Next you will need to compile the C++ code into MATLAB executable
("MEX") files. To do this, you will need to have a [C++ compiler
supported by
MATLAB](http://www.mathworks.com/support/compilers/current_release/),
and you will need to configure MATLAB to build MEX files. See [this
webpage](http://www.mathworks.com/support/tech-notes/1600/1605.html)
for details. When you follow this step, it is important that you
configure MATLAB so that it uses the version of the C++ compiler that
is compatible with MATLAB. Otherwise, you will encounter errors, or
MATLAB may crash when attempting to run the code. If you run into
these problems, you may have to run the mex command with the -v flag
to check what compiler is being used, and you may have to edit the MEX
configuration file manually.

To build the necessary MEX files, run the **install.m** script in
MATLAB.

###Quick start for R

*Instructions go here.*

###Overview of MATLAB functions 

*Details about main functions go here.*

###Overview of R functions

*Details about main functions go here.*

###Who

Variational inference for Bayesian variable selection was developed by:<br>
[Peter Carbonetto]((http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
March 2012
