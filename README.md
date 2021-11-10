README
================

``` r
library(ls)
```

# 1\. Getting Started

  - To install, `devtools::install_github("eleafeit/latent_strat")`
  - This package is built from [Berman & Feit
    (2019)](https://arxiv.org/pdf/1911.08438.pdf)
  - Citation of the paper (APA): Berman, R. & Feit E.M. (2019, November
    19). Principal Stratification for Advertising Experiments. *The
    Wharton School Research Paper*, *Wharton Customer Analytics
    Initiative Research Paper*, 1-37.

# 2\. About the Package

This package contains the necessary functions of the latent
stratification model from Berman & Feit (2019).In advertising
experiments, noisy responses complicates the estimation of average
treatment effect (ATE) and the assessment of the return of interest
(ROI). Therefore, running normal t-tests will often yield estimates with
high variance. The latent stratification model divides the customers
into three strata: always buyer, who are positive under treatment and
control (A); influenced buyer, who are positive only under treatment
(B), and never buyer, who are zero under treatment and control (C). The
model assumes that defiers, who are positive only under control, does
not exist. It is able to improve the precision of the estimate by
lowering the variance of the estimate. At the same time, it is also able
to improve the accuracy by yielding an ATE closer to the true value.

# 3.1 Data Preparation

In order to demonstrate the package, we first simulate the data via the
`sim_latent_strat()` function. The number following the strata (A for
always buyer, B for influenced buyer, and C for never buyer) indicates
whether the subject receives treatment. In this case, if the number is
1, then it means that the subject has received treatment, and 0
otherwise. The function takes in the sample size (n), proportion of
always buyer (piA), proportion of influenced buyer (piB), mean of always
buyer who received treatment (muA1), mean of always buyers who are in
control, mean of the influenced buyer who received treatment, and the
variance of group A1, A0, and B1. We assume that the treatment
proportion is 50%. The `sim_latent_strat()` function generates a list
containing a data frame `data`, the ATE `ATE`, and a vector `par`.

``` r
set.seed(20030601)
sim <- sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
```

The input parameters are stored as the `par` vector. The data will be
simulated from these true parameters. For example, 20% of all
observations belongs to strata A, who are positive under treatment and
control. Observations in strata A who received treatment have a mean
outcome of 5, and those who are in control have a mean outcome of 4.5,
both with a variance of 0.3.

``` r
sim$par
```

    ##   piA   piB  muA1  muA0  muB1 sigma 
    ##   0.2   0.1   5.0   4.5   3.0   0.3

The true ATE can be calculated from these true parameters. In this
example, the true ATE is 0.4.

``` r
sim$ATE
```

    ## [1] 0.4

Let’s now examine the simulated data.

``` r
head(sim$data)
```

    ##          y z s   sB   sC
    ## 1 3.356879 1 B  1.1 -0.5
    ## 2 0.000000 1 C  0.1  0.5
    ## 3 0.000000 1 C  0.1  0.5
    ## 4 0.000000 1 C  0.1  0.5
    ## 5 3.167482 1 B  1.1 -0.5
    ## 6 4.749002 1 A -0.9 -1.5

The data frame contains the 10000 randomly generated value (y) for each
observation based on the strata (s; which can be A, always buyer; B,
influenced buyer; or C, never buyer) and treatment (z; 1 for treatment,
0 for control). Since the treatment proportion is 50%, the first 5000
observations will have z = 1. For example, in the first entry, the
subject belongs to strata B, who is an influenced buyer that only buys
under treatment. Therefore, the subject in the first entry received
treatment and has an outcome of 3.36. Following the same logic, the
second entry is a subject in the never buyer strata. Therefore, even
though he received treatment, his outcome is 0.

# 3.2 Estimating ATE Using Latent Stratification

To estimate ATE using latent stratification, we simply call the
`mle_ls()` function and input the data frame as a parameter.

``` r
test = mle_ls(sim$data)
```

    ## Optimization runs took 0.1 secs
    ## Greatest ll was achieved in run 1 at -7008.303 (worst was -7008.303 )
    ## Computing variance-covariance matrix took 0.1  secs

``` r
test$pars
```

    ##      par        est          se
    ## 1    ATE  0.4032336 0.013019384
    ## 2 expATE 13.6626614 0.401495726
    ## 3    piA  0.2000032 0.004000203
    ## 4    piB  0.1000800 0.004225939
    ## 5   muA1  4.9938814 0.009284889
    ## 6   muA0  4.4878030 0.009460486
    ## 7   muB1  3.0177472 0.013298883
    ## 8  sigma  0.2961601 0.004195862

The function returned the analysis result as a data frame. It includes
the estimates and standard errors of ATE, strata proportions, strata
outcome means, and variance. The latent stratification model provides an
accurate estimate of the ATE. This can be seen by ATE estimate of 0.403.
In addition, the estimate is also precise. The standard error of the
estimate under the latent stratification model is only 0.013.
