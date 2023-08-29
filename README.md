
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Build rd-file with devtools::build_readme() -->

# inlami <a href='https://github.com/emmaSkarstein/inlami'><img src='man/figures/inlami_transparent.png' align="right" height="131.5" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/emmaSkarstein/inlami/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/emmaSkarstein/inlami/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**NOTE** The package is still under development, and some important
functionality is still missing. Feel free to get in touch if you would
like to use the package, but would like to know more about it’s current
limitations.

------------------------------------------------------------------------

Fitting measurement error models and missing data imputation models in
INLA is not trivial, and requires several workarounds in order to fit
the model. At the same time, a Bayesian hierarchical framework is very
attractive for modeling these types of data, since the hierarchical
model structure alllows us to describe how the data were collected, and
any errors they may have, through the model specification. By supplying
informative priors, we may be able to adjust for biases due to these
errors, and propagate any uncertainty they cause. Until recently, it has
been complicated to implement these kinds of models in R-INLA. This
package provides a helpful interface that makes measurement error and
missing data modelling in R-INLA much more feasible.

## Installation

You can install the development version of inlami from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("emmaSkarstein/inlami")
```

## When should I use this package?

This package is designed for fitting models where you have one covariate
that has classical measurement error, Berkson measurement error, missing
observations, or any combination of these three. That could mean that
you only have missing data, and if so this package can do missing data
imputation.

The model itself must be of the class of models that is possible to fit
with R-INLA. That means that it can be used for most common regression
types, like linear regression and logistic regression, and you can
include as many error free covariates as needed. You can also include
random effects, the same way as you would normally include such effects
in R-INLA.

The package is unfortunately not able to adjust for measurement error or
missingness in multiple covariates, though this could be implemented, it
just hasn’t yet. Feel free to get in touch if that is functionality that
would be of interest for your research!

## How can I use this package?

The dataset `simple_data` is included in the package, and is a very
simple simulated data set, used to illustrate the package. In this data
set, we have a response variable $y$, an error free covariate $z$, and a
covariate that is observed with classical error, Berkson error and
missingness, called $x$. We wish to fit the model

$$
 y = \beta_0 + \beta_x x + \beta_z z + \varepsilon \ ,
$$

but adjusting for all the errors in $x$.

First load the package:

``` r
library(inlami)
```

Next, we need to specify the formula for the main model and the formula
for the imputation model. This is done in the standard way in *R*:

``` r
main_formula <- y ~ x + z
```

For the imputation model, we take advantage of any correlation between
$x$ and the error free covariate $z$, so the imputation model will be

$$
 x = \alpha_0 + \alpha_z z + \varepsilon_x \ .
$$

We write that as

``` r
imputation_formula <- x ~ z
```

When adjusting for measurement error, we are completely dependent on
having some information about the measurement error variaces
$\sigma_{u_c}^2$ (for the classical error) and $\sigma_{u_b}^2$ (for the
Berkson error), since the size of this variance will affect how the
estimates are biased. We can gain information about these variances in a
few different ways, if repeated measurements have been made then these
can be put directly into the model to estimate the error variance, or if
we have some expert knowledge about the error size, then that can be
used to specify an informative prior for the variance (or precision,
since in INLA the precision is used, rather than the variance).

In this case, since we have in fact simulated the data ourselves, we
know that the error variances are in both cases close to 1, so we
specify priors that have modes at 1.

In the `fit_inlami` function we also need to specify the likelihood for
the model of interest, which in this case is Gaussian.

``` r
simple_model <- fit_inlami(formula_moi = main_formula,
                           formula_imp = imputation_formula,
                           family_moi = "gaussian",
                           data = simple_data, 
                           error_type = c("classical", "berkson", "missing"),
                           prior.prec.y = c(10, 9),
                           prior.prec.u_b = c(10, 9),
                           prior.prec.u_c = c(10, 9),
                           prior.prec.r = c(10, 9),
                           initial.prec.y = 1,
                           initial.prec.u_b = 1,
                           initial.prec.u_c = 1,
                           initial.prec.r = 1)
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
```

Once we have fit the model, we can view the summary:

``` r
summary(simple_model)
#> Fixed effects for model of interest: 
#>            mean        sd 0.025quant 0.5quant 0.975quant     mode          kld
#> beta.0 1.060777 0.2160778  0.6607368 1.061048   1.480327 1.107720 3.394914e-07
#> beta.z 1.974736 0.3845800  1.3292609 1.994203   2.725384 2.074372 6.295227e-05
#> 
#> Coefficient for error prone variable: 
#>            mean        sd 0.025quant 0.5quant 0.975quant     mode
#> beta.x 1.934415 0.1872919   1.584007 1.927663   2.321449 1.901174
#> 
#> Fixed effects for imputation model: 
#>             mean         sd 0.025quant 0.5quant 0.975quant     mode
#> alpha.0 1.033112 0.05053598  0.9339588 1.033119   1.132223 1.033134
#> alpha.z 2.024636 0.05219671  1.9222915 2.024620   2.127071 2.024588
#>                  kld
#> alpha.0 2.108388e-12
#> alpha.z 9.127922e-12
#> 
#> Model hyperparameters (apart from beta.x): 
#>                                                 mean        sd 0.025quant
#> Precision for the Gaussian observations    1.1019393 0.3310944  0.5695484
#> Precision for the Gaussian observations[2] 1.0525910 0.3198897  0.5878718
#> Precision for the Gaussian observations[3] 0.9392355 0.1110450  0.7315933
#> Precision for the Gaussian observations[4] 0.9695980 0.1193234  0.7621798
#>                                             0.5quant 0.975quant      mode
#> Precision for the Gaussian observations    1.0632758   1.861922 0.9908255
#> Precision for the Gaussian observations[2] 0.9980265   1.834041 0.8917934
#> Precision for the Gaussian observations[3] 0.9357206   1.168491 0.9338017
#> Precision for the Gaussian observations[4] 0.9599274   1.231860 0.9371894
```

And we can use the default plot function to see a plot of the fixed
effects:

``` r
plot(simple_model)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />
