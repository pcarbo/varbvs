# Large-scale Bayesian variable selection for R and MATLAB

[![CRAN status badge](http://www.r-pkg.org/badges/version/varbvs)](http://cran.r-project.org/package=varbvs)
[![Travis Build Status](https://travis-ci.org/pcarbo/varbvs.svg?branch=master)](https://travis-ci.org/pcarbo/varbvs)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/xqnqhsk5na3855nq?svg=true)](https://ci.appveyor.com/project/pcarbo/varbvs)
[![codecov](https://codecov.io/gh/pcarbo/varbvs/branch/master/graph/badge.svg)](https://codecov.io/gh/pcarbo/varbvs)

### Overview

We introduce *varbvs*, a suite of functions writen in
[R](http://www.r-project.org) and
[MATLAB](http://www.mathworks.com/products/matlab) for analysis of
large-scale data sets using Bayesian variable selection methods. To
facilitate application of Bayesian variable selection to a range of
problems, the *varbvs* interface hides most of the complexities of
modeling and optimization, while also providing many options for
adaptation to range of applications. The *varbvs* software has been
used to implement Bayesian variable selection for large problems with
over a million variables and thousands of samples, including analysis
of massive genome-wide data sets.

For more details on the R package, see here the
[README](varbvs-R/docs/README.md) and the
[online documentation](http://pcarbo.github.io/varbvs/docs/index.html).
For example, only a few lines of R code are needed to
fit a variable selection model to the Leukemia data:

```R
library(varbvs)
data(leukemia)
fit <- varbvs(leukemia$x,NULL,leukemia$y,family = "binomial",
              logodds = seq(-3.5,-1,0.1),sa = 1)
print(summary(fit))
```

For more details on the MATLAB interface, see the
[README](varbvs-MATLAB/README.md) in the
[varbvs-MATLAB](varbvs-MATLAB) subdirectory. The MATLAB package also
provides a simple interface for fitting variable selection models;
this is the call to varbvs from the
[demo_qtl.m](varbvs-MATLAB/demo_qtl.m) example:

```MATLAB
fit = varbvs(X,Z,y,labels,[],struct('logodds',-3:0.1:-1));
```

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
fitness for a particular purpose**. See file [LICENSE](LICENSE) for
the full text of the license.

### Credits

The *varbvs* software package was developed by:<br>
[Peter Carbonetto](http://pcarbo.github.io)<br>
Dept. of Human Genetics, University of Chicago<br>
2012-2017

Xiang Zhou, Xiang Zhu, David Gerard and Matthew Stephens have also
contributed to the development of this software.
