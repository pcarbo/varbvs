#Large-scale Bayesian variable selection for R and MATLAB
 
###Introduction

We introduce *varbvs*, a suite of functions writen in R and MATLAB for
analysis of large-scale data sets using Bayesian variable selection
methods. To facilitate application of Bayesian variable selection to a
range of problems, the *varbvs* interface hides most of the
complexities of modeling and optimization, while also providing many
options for adaptation to range of applications. The *varbvs* software
has been used to implement Bayesian variable selection for large
problems with over a million variables and thousands of samples,
including analysis of massive genome-wide data sets.

The MATLAB interface has been tested in version 8.6.0 (2015b). The R
package has been tested in version R versions 3.3.1 and 3.3.2.

If you find that this software is useful for your research project,
please cite our paper:

Carbonetto, P., and Stephens, M. (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. *Bayesian Analysis* **7**,
73--108.

###License

Copyright (c) 2012-2016, Peter Carbonetto.

The *varbvs* source code repository by
[Peter Carbonetto](http://github.com/pcarbo) is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *varbvs*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](LICENSE) for
the full text of the license.

###Quick start for R

Start by downloading the github repository for this project. The
simplest way to do this is to [download the repository as a ZIP
archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
you have extracted the files from the compressed archive, you will see
that the main directory has two subdirectories, one containing the
MATLAB code, and the other containing R files.

The subdirectory **R/varbvs** has all the necessary files to build and
install a package for R. To install this package, follow the [standard
instructions](http://cran.r-project.org/doc/manuals/R-admin.html) for
installing an R package from source. On a Unix or Unix-like platform
(such as Mac OS X), the installation command is
**R CMD INSTALL --clean homedir/varbvs/R/varbvs**.

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

**Note:** The demonstration script
[example1.R](R/varbvs/demo/example1.R) requires packages
[ggplot2](http://had.co.nz/ggplot2) and **grid** to show the plots
depicting the posterior distributions of the hyperparameters computed
using the variational inference method. The grid package is normally
included in the base distribution of R, and ggplot2 is available on
[CRAN](http://cran.r-project.org), and can be installed using the
**install.packages** function (for more recent versions of R).

**Another note:** By default, R compiles the C code with the debugging
information using the **-g** flag. The computations will run much
faster without this debugging information. For example, on a Linux
machine you can override the default package compilation by creating a
text file **.R/Makevars** in your home directory that contains a
single line: **CFLAGS=-fpic -O3**. 

###Quick start for MATLAB

Start by downloading the github repository for this project. The
simplest way to do this is to [download the repository as a ZIP
archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
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

To build the necessary MEX files, run the
[install.m](MATLAB/install.m) script in MATLAB.

Note that the beginning of this script sets some compiler and linker
flags. These flags tell the GCC compiler to use the ISO C99 standard,
and to optimize the code as much as possible. However, these flags may
not be relevant to your setup, especially if you are not using
[gcc](http://gcc.gnu.org). To avoid errors during installation, if you
are using a compiler other than GCC, it may be best to set variables
**cflags** and **ldflags** to empty strings before running the
[install.m](MATLAB/install.m) script.

If you modify the installation procedure to fit your compiler setup,
it is important that you define macro MATLAB_MEX_FILE, akin to the
\#define MATLAB_MEX_FILE directive in C. In GCC, this is accomplished
by including flag -DMATLAB_MEX_FILE when issuing the commands to build
MEX files.

Once you have built the MEX files, start by running the script
[example1.m](MATLAB/example1.m). This script demonstrates how the
variational inference algorithm is used to compute posterior
probabilities for a small linear regression example in which only a
small subset of the variables (single nucleotide polymorphisms, or
SNPs) has affects the outcome (a simulated quantitative trait). In
this small example, the variational estimates of the posterior
probabilities are compared with estimates obtained by MCMC
simulation. Notice that it takes a considerable amount of time to
simulate the Markov chain.

###Credits

The *varbvs* software package was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Department of Human Genetics, University of Chicago<br>
and AncestryDNA, San Francisco, CA<br>
2012-2016

Other important contributors to this software include Xiang Zhu and
Matthew Stephens.
