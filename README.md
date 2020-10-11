
# R/`natmed2`

[![Travis-CI Build
Status](https://travis-ci.org/benkeser/natmed2.svg?branch=master)](https://travis-ci.org/benkeser/natmed2)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/benkeser/natmed2?branch=master&svg=true)](https://ci.appveyor.com/project/benkeser/natmed2)
[![Coverage
Status](https://img.shields.io/codecov/c/github/benkeser/natmed2/master.svg)](https://codecov.io/github/benkeser/natmed2?branch=master)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Multiply robust estimators of natural mediation parameters for designs
> where the measurement of the mediator is subject to two-phase
> sampling.

**Author:** [David
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)

-----

## Description

`natmed2` is an R package that computes multiply robust estimators of
natural mediation effects in settings where the mediator is measured
subject to two-phase sampling.

-----

## Installation

<!-- 
Install the current stable release from
[CRAN](https://cran.r-project.org/) via


```r
install.packages("natmed2")
```
 -->

A developmental release may be installed from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("benkeser/natmed2")
```

-----

## Usage

Below is an illustration of the the main function on simulated data.

``` r
n <- 500
W1 <- rbinom(n, 1, 0.5)
W2 <- rnorm(n, 0, 1)
A <- rbinom(n, 1, 0.5)
S <- W1 / 4 - W2 / 3 + A + rnorm(n)
Y <- rbinom(n, 1, plogis(-2 + A + W1 / 2 - S / 2))

# add censoring
C <- rbinom(n, 1, plogis(2 + W1 / 2 - W2 / 3))
# arbitrary fill in
Y[C == 0] <- -999

R <- rep(0, n)
# case-cohort sampling
R <- rbinom(n, 1, 0.25)
R[Y == 1] <- 1

library(natmed2)
fit <- natmed2(
  W = data.frame(W1 = W1, W2 = W2), 
  A = A, R = R, S = S, C = C, Y = Y
)
#> Warning in eval(family$initialize): non-integer #successes in a binomial
#> glm!
#> Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
#> ifelse(type == : prediction from a rank-deficient fit may be misleading

#> Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
#> ifelse(type == : prediction from a rank-deficient fit may be misleading
#> Warning in eval(family$initialize): non-integer #successes in a binomial
#> glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial
#> glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial
#> glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial
#> glm!

fit
#> $risk
#>               one_step         cil       ciu      cil_cv    ciu_cv
#> E[Y(1,S(1))] 0.2249777  0.17014223 0.2798132  0.17014223 0.2798132
#> E[Y(0,S(0))] 0.1494254  0.10314521 0.1957055  0.10314521 0.1957055
#> E[Y(1,S(0))] 0.1579031 -0.01509897 0.3309051 -0.01509897 0.3309051
#> E[Y(0,S(1))] 0.3840393  0.17889991 0.5891786  0.17889991 0.5891786
#> 
#> $eff
#>                               effect one_step_est       cil      ciu
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]     1.505619 1.0161249 2.230916
#> Direct   E[Y(1,S(0))] / E[Y(0,S(0))]     1.056735 0.3390251 3.293827
#> Indirect E[Y(1,S(1))] / E[Y(1,S(0))]     1.424784 0.5160847 3.933479
#>              cil.1    ciu.1
#> Total    1.0161249 2.230916
#> Direct   0.3390251 3.293827
#> Indirect 0.5160847 3.933479
#> 
#> $eff2
#>                               effect one_step_est       cil      ciu
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]    1.5056194 1.0161249 2.230916
#> Direct   E[Y(1,S(1))] / E[Y(0,S(1))]    0.5858196 0.3257179 1.053625
#> Indirect E[Y(0,S(1))] / E[Y(0,S(0))]    2.5701076 1.4251442 4.634936
#>              cil.1    ciu.1
#> Total    1.0161249 2.230916
#> Direct   0.3257179 1.053625
#> Indirect 1.4251442 4.634936
#> 
#> $cov
#>              eif_psi11    eif_psi00    eif_psi10    eif_psi01
#> eif_psi11 7.827292e-04 3.160710e-06 1.056723e-03 2.155223e-06
#> eif_psi00 3.160710e-06 5.575418e-04 1.182496e-05 2.505094e-04
#> eif_psi10 1.056723e-03 1.182496e-05 7.790949e-03 6.748625e-05
#> eif_psi01 2.155223e-06 2.505094e-04 6.748625e-05 1.095433e-02
#> 
#> $cov_cv
#>              eif_psi11    eif_psi00    eif_psi10    eif_psi01
#> eif_psi11 7.827292e-04 3.160710e-06 1.056723e-03 2.155223e-06
#> eif_psi00 3.160710e-06 5.575418e-04 1.182496e-05 2.505094e-04
#> eif_psi10 1.056723e-03 1.182496e-05 7.790949e-03 6.748625e-05
#> eif_psi01 2.155223e-06 2.505094e-04 6.748625e-05 1.095433e-02
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/benkeser/natmed2/issues).

-----

## License

© 2016-2019 [David C.
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)
    
    Copyright (c) 2016-2019 David C. Benkeser
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
