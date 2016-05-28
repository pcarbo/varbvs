#Large-scale Bayesian variable selection for R and MATLAB
 
###Overview

We introduce *varbvs*, a suite of functions writen in R and MATLAB for
analysis of large-scale data sets using Bayesian variable selection
methods. To facilitate application of Bayesian variable selection to a
range of problems, the *varbvs* interface hides most of the
complexities of modeling and optimization, while also providing many
options for adaptation to range of applications. The *varbvs* software
has been used to implement Bayesian variable selection for large
problems with over a million variables and thousands of samples,
including analysis of massive genome-wide data sets.

The MATLAB interface has been tested extensively in MATLAB
version 8.6.0 (2015b). The R package has also been tested extensively
in R versions 3.2.1 and 3.2.2.

###Citing varbvs

If you find that this software is useful for your research project,
please cite our paper:

Carbonetto, P., and Stephens, M. (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. *Bayesian Analysis* **7**,
73-108.

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

To install the official release of the *varbvs* package available from
CRAN ([link](http://cran.r-project.org/web/packages/varbvs)), simply
run <code>install.packages("varbvs")</code> in R.

Alternatively, to install the most up-to-date development version,
begin by downloading the github repository for this project. The
simplest way to do this is to
[download the repository as a ZIP archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
you have extracted the files from the compressed archive, you will see
that the main directory has two subdirectories, one containing the
MATLAB code, and the other containing the files for the R package.

Subdirectory [varbvs-R](varbvs-R) has all the necessary files to build
and install a package for R. To install this package, follow the
[standard instructions](http://cran.r-project.org/doc/manuals/R-admin.html)
for installing an R package from source. On a Unix or Unix-like
platform (e.g., Mac OS X), the following steps should install the R
package:

    mv varbvs-R varbvs
	R CMD build varbvs
	R CMD INSTALL varbvs_2.0.0.tar.gz

Once you have installed the package, load the package in R by typing
<code>library(varbvs)</code>. To get an overview of the package, type
<code>help(package="varbvs")</code>. The most important function you
will use is function "varbvs". Type <code>help(varbvs)</code> to get
more information about this.

We have provided several R scripts in the
[vignettes](varbvs-R/vignettes) folder to illustrate application of
*varbvs* to small and large data sets:

+ Script [demo.qtl.R](varbvs-R/vignettes/demo.qtl.R) demonstrates how to
use the varbvs function for mapping a quantitative trait (*i.e.*, a
continuously valued outcome) in a small, simulated data set. Script
[demo.cc.R](varbvs-R/vignettes/demo.cc.R) demonstrates mapping of a
binary valued outcome in a simulated data set.

+ Script [demo.leukemia.R](varbvs-R/vignettes/demo.leukemia.R)
demonstrates application of both *glmnet* and *varbvs* to the Leukemia
data. The main aim of this script is to illustrate some of the
different properties of *varbvs* (Bayesian variable selection) and
*glmnet* (penalized sparse regression). This script also reproduces
the results and graphs presented in the first example of Carbonetto *et
al* (2016).

+ Like demo.qtl.R, script [demo.cfw.R](varbvs-R/vignettes/demo.cfw.R)
also demonstrates *varbvs* for mapping genetic factors contributing to
a quantitative trait, but here it is applied to an actual data set
generated from an outbred mouse study. Running this script with
<code>trait = "testis"</code> reproduces the results and figures given
in the second example of Carbonetto *et al* (2016).

+ Finally, scripts [demo.cd.R](varbvs-R/vignettes/demo.cd.R) and
[demo.cytokine.R](varbvs-R/vignettes/demo.cytokine.R) show how the
*varbvs* package can be applied to a very large data set to map
genetic loci and test biological hypotheses about genetic factors
contributing to human disease risk. Although we cannot share the data
needed to run these scripts due to data privacy restrictions, we have
included these scripts because it is helpful to be able to follow the
steps given in these R scripts. These scripts reproduce some of the
results and figures presented in Carbonetto *et al* (2016).

###Quick start for MATLAB

Begin by downloading the github repository for this project. The
simplest way to do this is to
[download the repository as a ZIP archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
you have extracted the files from the compressed archive, you will see
that the main directory contains two subdirectories, one containing
the MATLAB code, and the other containing the files for the R package.

Next, you will need to compile the C code into MATLAB executable
("MEX") files. To build the necessary MEX files, run the
[install.m](varbvs-MATLAB/install.m) script in MATLAB. For the MEX
files to be built successfully, you will need to have a
[C compiler supported by MATLAB](http://www.mathworks.com/support/compilers/current_release/),
and you will need to configure MATLAB to build MEX files. See
[this webpage](http://www.mathworks.com/support/tech-notes/1600/1605.html)
for more details.

When you take this step, it is important that you configure MATLAB so
that it uses the version of the C compiler that is compatible with
your version of MATLAB. Otherwise, you will encounter errors when
building the MEX files, or MATLAB may crash when attempting to run the
examples. If you run into these problems, you may have to run the mex
command with the -v flag to check what compiler is being used, and you
may have to edit the MEX configuration file manually.

Also note that the beginning of this script sets some compiler and
linker flags. These flags tell the GCC compiler to use the ISO C99
standard, and to optimize the code as much as possible. However, these
flags may not be relevant to your setup, especially if you are not
using [gcc](http://gcc.gnu.org). Do not remove flag -DMATLAB_MEX_FILE;
this is important for correctly compiling the C code for MATLAB.

We have provided a few MATLAB scripts to illustrate application of
*varbvs* to small and large data sets:

+ Script [demo_qtl.m](varbvs-MATLAB/demo_qtl.m) demonstrates how to
use the varbvs function for mapping a quantitative trait (*i.e.*, a
continuously valued outcome) in a small, simulated data
set. Additionally, script [demo_cc.m](varbvs-MATLAB/demo_cc.m)
demonstrates mapping of a binary valued outcome in a simulated data
set.

+ Scripts [demo_cd.m](varbvs-MATLAB/demo_cd.m),
[demo_celiac.m](varbvs-MATLAB/demo_celiac.m) and
[demo_cytokine.m](varbvs-MATLAB/demo_cytokine.m) show how the *varbvs*
package can be applied to very large data sets to map genetic loci and
test biological hypotheses about genetic factors contributing to human
diseases. Although we cannot share the data needed to run these
scripts due to data privacy restrictions, we have included these
scripts anyhow since it is helpful to be able to follow the steps
given in these MATLAB scripts. These scripts reproduce some of the
results and figures presented in Carbonetto *et al* (2016).

###Credits

The *varbvs* software package was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics, University of Chicago<br>
and AncestryDNA, San Francisco, California<br>
2012-2016

Xiang Zhou, Xiang Zhu and Matthew Stephens have also contributed to
the development of this software.
