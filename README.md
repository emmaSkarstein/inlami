
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Build rd-file with devtools::build_readme() -->

# inlami <a href='https://github.com/emmaSkarstein/inlami'><img src='man/figures/inlami_transparent.png' style="float:right; height:200px;" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/emmaSkarstein/inlami/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/emmaSkarstein/inlami/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> :warning: **NOTE** The package is under development, and some
> important functionality is still missing. Feel free to get in touch if
> you would like to use the package, but would like to know more about
> it’s current limitations.

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
simple_model <- fit_inlami(data = simple_data, 
                         formula_moi = main_formula, 
                         formula_imp = imputation_formula, 
                         family_moi = "gaussian",
                         error_type = c("berkson", "classical"),
                         prior.prec.moi = c(10, 9),
                         prior.prec.berkson = c(10, 9), 
                         prior.prec.classical = c(10, 9),
                         prior.prec.imp = c(10, 9),
                         initial.prec.moi = 1,
                         initial.prec.berkson = 1,
                         initial.prec.classical = 1,
                         initial.prec.imp = 1)
```

Once we have fit the model, we can view the summary:

``` r
summary(simple_model)
#> Formula for model of interest: 
#> y ~ x + z
#> 
#> Formula for imputation model: 
#> x ~ z
#> 
#> Error types: 
#> [1] "berkson"   "classical"
#> 
#> Fixed effects for model of interest: 
#>            mean        sd 0.025quant 0.5quant 0.975quant mode          kld
#> beta.0 1.043220 0.2124403  0.6356073 1.043788   1.463736   NA 1.630568e-07
#> beta.z 1.939786 0.3776542  1.2496757 1.946568   2.705101   NA 8.641776e-06
#> 
#> Coefficient for error prone variable: 
#>            mean        sd 0.025quant 0.5quant 0.975quant mode
#> beta.x 1.940827 0.1769273   1.601524 1.937251   2.299834   NA
#> 
#> Fixed effects for imputation model: 
#>             mean         sd 0.025quant 0.5quant 0.975quant mode          kld
#> alpha.0 1.033071 0.05062282  0.9337409 1.033080   1.132348   NA 3.379187e-12
#> alpha.z 2.024728 0.05228673  1.9222192 2.024707   2.127353   NA 1.494428e-11
#> 
#> Model hyperparameters (apart from beta.x): 
#>                                                 mean        sd 0.025quant
#> Precision for the Gaussian observations    1.1348399 0.3469250  0.5811453
#> Precision for the Gaussian observations[2] 1.0492040 0.2794852  0.6201420
#> Precision for the Gaussian observations[3] 0.9400527 0.1102513  0.7362982
#> Precision for the Gaussian observations[4] 0.9599192 0.1204817  0.7501751
#>                                             0.5quant 0.975quant mode
#> Precision for the Gaussian observations    1.0930371   1.934834   NA
#> Precision for the Gaussian observations[2] 1.0085920   1.714990   NA
#> Precision for the Gaussian observations[3] 0.9355375   1.170696   NA
#> Precision for the Gaussian observations[4] 0.9502789   1.224720   NA
```

And we can use the default plot function to see a plot of the fixed
effects and estimated coefficient for the variable with error:

``` r
plot(simple_model)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />
