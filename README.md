
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
#>                one_step         cil       ciu      cil_cv    ciu_cv
#> E[Y(1,S(1))] 0.22361271  0.17074847 0.2764769  0.17074847 0.2764769
#> E[Y(0,S(0))] 0.14736160  0.09832274 0.1964005  0.09832274 0.1964005
#> E[Y(1,S(0))] 0.07937785 -0.08523596 0.2439917 -0.08523596 0.2439917
#> E[Y(0,S(1))] 0.31651220  0.20578744 0.4272370  0.20578744 0.4272370
#> 
#> $eff
#>                               effect one_step_est        cil       ciu
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]    1.5174421 1.00916832  2.281711
#> Direct   E[Y(1,S(0))] / E[Y(0,S(0))]    0.5386603 0.06593618  4.400542
#> Indirect E[Y(1,S(1))] / E[Y(1,S(0))]    2.8170669 0.37722321 21.037586
#>              cil_cv    ciu_cv
#> Total    1.00916832  2.281711
#> Direct   0.06593618  4.400542
#> Indirect 0.37722321 21.037586
#> 
#> $eff2
#>                               effect one_step_est       cil      ciu
#> Total    E[Y(1,S(1))] / E[Y(0,S(0))]     1.517442 1.0091683 2.281711
#> Direct   E[Y(1,S(1))] / E[Y(0,S(1))]     0.706490 0.4630251 1.077972
#> Indirect E[Y(0,S(1))] / E[Y(0,S(0))]     2.147861 1.4429147 3.197213
#>             cil_cv   ciu_cv
#> Total    1.0091683 2.281711
#> Direct   0.4630251 1.077972
#> Indirect 1.4429147 3.197213
#> 
#> $cov
#>               eif_psi11     eif_psi00     eif_psi10     eif_psi01
#> eif_psi11  7.274645e-04  1.070057e-06  7.253333e-04 -2.400833e-06
#> eif_psi00  1.070057e-06  6.259919e-04 -4.285169e-07  4.545122e-04
#> eif_psi10  7.253333e-04 -4.285169e-07  7.053755e-03 -1.765371e-04
#> eif_psi01 -2.400833e-06  4.545122e-04 -1.765371e-04  3.191371e-03
#> 
#> $cov_cv
#>               eif_psi11     eif_psi00     eif_psi10     eif_psi01
#> eif_psi11  7.274645e-04  1.070057e-06  7.253333e-04 -2.400833e-06
#> eif_psi00  1.070057e-06  6.259919e-04 -4.285169e-07  4.545122e-04
#> eif_psi10  7.253333e-04 -4.285169e-07  7.053755e-03 -1.765371e-04
#> eif_psi01 -2.400833e-06  4.545122e-04 -1.765371e-04  3.191371e-03
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
