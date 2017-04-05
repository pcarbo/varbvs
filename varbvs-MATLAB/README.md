# Large-scale Bayesian variable selection for MATLAB

### Overview

We introduce *varbvs*, a suite of MATLAB routines for analysis of
large-scale data sets using Bayesian variable selection methods. To
facilitate application of Bayesian variable selection to a range of
problems, the *varbvs* interface hides most of the complexities of
modeling and optimization, while also providing many options for
adaptation to range of applications. The *varbvs* software has been
used to implement Bayesian variable selection for large problems with
over a million variables and thousands of samples, including analysis
of massive genome-wide data sets.

The MATLAB interface has been tested extensively in MATLAB
version 8.6.0 (2015b). 

### Citing varbvs

If you find that this software is useful for your research project,
please cite our paper:

Carbonetto, P., and Stephens, M. (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. *Bayesian Analysis* **7**,
73-108.

### License

Copyright (c) 2012-2017, Peter Carbonetto.

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

Begin by downloading the github repository for this project. The
simplest way to do this is to
[download the repository as a ZIP archive](http://github.com/pcarbo/varbvs/archive/master.zip). Once
you have extracted the files from the compressed archive, you will see
that the main directory contains two subdirectories, one containing
the MATLAB code, and the other containing the files for the R package.

Next, you will need to compile the C code into MATLAB executable
("MEX") files. To build the necessary MEX files, run the
[install.m](install.m) script in MATLAB. For the MEX
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

### Using the package

The main function of this package is the <code>varbvs</code> function,
which fits the variable selection model to data. In the
[demo_qtl.m](demo_qtl.m) script, for example, the varbvs function call
is simply

```MATLAB
fit = varbvs(X,Z,y,labels,[],struct('logodds',(-3:0.1:-1)));
```

In this example, we have used varbvs to map a quantitative trait
(*i.e.*, a continuously valued outcome) in a small, simulated data
set. Additionally, script [demo_cc.m](demo_cc.m) demonstrates mapping
of a binary valued outcome in a simulated data set.

We have provided a few other MATLAB scripts to demonstrate the
application of *varbvs* to very large data sets:
[demo_cd.m](demo_cd.m), [demo_celiac.m](demo_celiac.m) and
[demo_cytokine.m](demo_cytokine.m). Although we cannot share the data
needed to run these scripts due to data privacy restrictions, we have
included these scripts anyhow since it is helpful to be able to follow
the steps given in these MATLAB scripts. These scripts reproduce some
of the results and figures presented in Carbonetto *et al* (2016).

### Credits

The *varbvs* software package was developed by:<br>
[Peter Carbonetto](http://pcarbo.github.io)<br>
Dept. of Human Genetics, University of Chicago<br>
2012-2017

Xiang Zhou, Xiang Zhu and Matthew Stephens have also contributed to
the development of this software.
