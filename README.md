
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

------------------------------------------------------------------------

## Description

`natmed2` is an R package that computes multiply robust estimators of
natural mediation effects in settings where the mediator is measured
subject to two-phase sampling.

------------------------------------------------------------------------

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

------------------------------------------------------------------------

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
#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!
#> Warning in predict.lm(object, newdata, se.fit, scale = 1, type = if (type == :
#> prediction from a rank-deficient fit may be misleading

#> Warning in predict.lm(object, newdata, se.fit, scale = 1, type = if (type == :
#> prediction from a rank-deficient fit may be misleading
#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!

fit
#> $risk
#>                one_step        cil       ciu     cil_cv    ciu_cv
#> E[Y(1,S(1))] 0.24561975 0.18962056 0.3016189 0.18962056 0.3016189
#> E[Y(0,S(0))] 0.15876201 0.10968346 0.2078405 0.10968346 0.2078405
#> E[Y(1,S(0))] 0.29097314 0.13876900 0.4431773 0.13876900 0.4431773
#> E[Y(0,S(1))] 0.09114202 0.03818788 0.1440962 0.03818788 0.1440962
#> 
#> $eff
#>                               effect one_step_est       cil      ciu    cil_cv
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]     1.547094 1.0555146 2.267614 1.0555146
#> Direct   E[Y(1,S(0))] / E[Y(0,S(0))]     1.832763 1.0022209 3.351577 1.0022209
#> Indirect E[Y(1,S(1))] / E[Y(1,S(0))]     0.844132 0.5237696 1.360443 0.5237696
#>            ciu_cv
#> Total    2.267614
#> Direct   3.351577
#> Indirect 1.360443
#> 
#> $eff2
#>                               effect one_step_est       cil       ciu    cil_cv
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]    1.5470940 1.0555146 2.2676141 1.0555146
#> Direct   E[Y(1,S(1))] / E[Y(0,S(1))]    2.6949124 1.4500698 5.0084159 1.4500698
#> Indirect E[Y(0,S(1))] / E[Y(0,S(0))]    0.5740795 0.3586164 0.9189966 0.3586164
#>             ciu_cv
#> Total    2.2676141
#> Direct   5.0084159
#> Indirect 0.9189966
#> 
#> $cov
#>              eif_psi11    eif_psi00     eif_psi10     eif_psi01
#> eif_psi11 8.163029e-04 6.860484e-06  9.099670e-04  1.589128e-05
#> eif_psi00 6.860484e-06 6.270052e-04  2.911345e-05  3.987886e-04
#> eif_psi10 9.099670e-04 2.911345e-05  6.030326e-03 -1.606312e-04
#> eif_psi01 1.589128e-05 3.987886e-04 -1.606312e-04  7.299408e-04
#> 
#> $cov_cv
#>              eif_psi11_cv eif_psi00_cv  eif_psi10_cv  eif_psi01_cv
#> eif_psi11_cv 8.163029e-04 6.860484e-06  9.099670e-04  1.589128e-05
#> eif_psi00_cv 6.860484e-06 6.270052e-04  2.911345e-05  3.987886e-04
#> eif_psi10_cv 9.099670e-04 2.911345e-05  6.030326e-03 -1.606312e-04
#> eif_psi01_cv 1.589128e-05 3.987886e-04 -1.606312e-04  7.299408e-04
#> 
#> $risk_lazy
#>                one_step       cil       ciu    cil_cv    ciu_cv
#> E[Y(1,S(1))] 0.24561975 0.1896206 0.3016189 0.1896206 0.3016189
#> E[Y(0,S(0))] 0.15876201 0.1096835 0.2078405 0.1096835 0.2078405
#> E[Y(1,S(0))] 0.28818635 0.1151032 0.4612695 0.1151032 0.4612695
#> E[Y(0,S(1))] 0.09338338 0.0242826 0.1624842 0.0242826 0.1624842
#> 
#> $eff_lazy
#>                               effect one_step_est       cil      ciu    cil_cv
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]    1.5470940 1.0555146 2.267614 1.0555146
#> Direct   E[Y(1,S(0))] / E[Y(0,S(0))]    1.8152098 0.9252997 3.560994 0.9252997
#> Indirect E[Y(1,S(1))] / E[Y(1,S(0))]    0.8522949 0.4879395 1.488722 0.4879395
#>            ciu_cv
#> Total    2.267614
#> Direct   3.560994
#> Indirect 1.488722
#> 
#> $eff2_lazy
#>                               effect one_step_est       cil      ciu    cil_cv
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]    1.5470940 1.0555146 2.267614 1.0555146
#> Direct   E[Y(1,S(1))] / E[Y(0,S(1))]    2.6302299 1.2157160 5.690564 1.2157160
#> Indirect E[Y(0,S(1))] / E[Y(0,S(0))]    0.5881973 0.3017323 1.146632 0.3017323
#>            ciu_cv
#> Total    2.267614
#> Direct   5.690564
#> Indirect 1.146632
#> 
#> $cov_lazy
#>                   eif_psi11    eif_psi00 eif_psi10_lazy eif_psi01_lazy
#> eif_psi11      8.163029e-04 6.860484e-06   9.361999e-04   1.180446e-05
#> eif_psi00      6.860484e-06 6.270052e-04   1.321039e-05   3.811615e-04
#> eif_psi10_lazy 9.361999e-04 1.321039e-05   7.798256e-03  -1.243435e-03
#> eif_psi01_lazy 1.180446e-05 3.811615e-04  -1.243435e-03   1.242950e-03
#> 
#> $cov_lazy_cv
#>                   eif_psi11_cv eif_psi00_cv eif_psi10_lazy_cv eif_psi01_lazy_cv
#> eif_psi11_cv      8.163029e-04 6.860484e-06      9.361999e-04      1.180446e-05
#> eif_psi00_cv      6.860484e-06 6.270052e-04      1.321039e-05      3.811615e-04
#> eif_psi10_lazy_cv 9.361999e-04 1.321039e-05      7.798256e-03     -1.243435e-03
#> eif_psi01_lazy_cv 1.180446e-05 3.811615e-04     -1.243435e-03      1.242950e-03
```

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/benkeser/natmed2/issues).

------------------------------------------------------------------------

## License

© 2021- [David C. Benkeser](https://www.davidbphd.com)

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)

    Copyright (c) 2021- David C. Benkeser

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
