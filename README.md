

![logo](www/logo.png){width=30%}

# recalc

{recalc} is a forensic meta-science R package. It allows you to recalculate *p*-values and standardized effect sizes (e.g., Cohen's *d*) from reported summary statistics and test statistics, taking a) the rounding of the inputs and outputs and b) various analytic choices into account (e.g., Students t vs Welch's t; whether SEs may have been confused with SDs, etc).  



## Inspired by

As a recalculation method, recalc is conceptually related to other methods such as [statcheck](https://statcheck.io/).

recalc was inspired by:

- rpsychi R package (v0.8, no longer on CRAN due to lack of maintenance)
- Nick Brown's f_range() function (https://steamtraen.blogspot.com/2018/01/)



## Related projects

Mark Bolland's {reappraised} R package and [Shiny app](https://reappraised.shinyapps.io/check_p_vals_cont/) aims aims to do the same as this package, with different degress of analtyic flexibility. 

- different versions and implementtions of between-groups cohen's d: https://rpubs.com/metinbulus/welch
- different versions and implementations of within-subjects cohen's d: https://github.com/ianhussey/versions-of-cohens-d
- Extracting pre-post correlations from summary stats: https://matthewbjane.quarto.pub/pre-post-correlations/
- Lisa Debruine's {faux} package has some error detection functions https://github.com/debruine/faux
- Lisa Debruine's {within} repo assesses the plausibility of a within subject t-test by showing the range of possible between timepoint correlations given the reported results https://github.com/debruine/within



## Installation

Currently, recalc is only available on GitHub. To install it, run the following code in your R Console:

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("ianhussey/recalc")
```



## Web app

Shiny app available [here](https://errors.shinyapps.io/recalc_independent_t_test/).



## TODO

- Expand to include dependent t-test. The 'within' repo has some code for this, but the basis should be the code in this repo. Need to add choice of also using the reported t value, or p value, or pre-post correlation, or a reasonable estimate of this (eg .5 to .7, or perhaps .25 to .9) to be able to reproduce the p and d. produce estimates of both d_z and d_av. 
- Use independent t test code to expand to one-way and two-way ANOVA, following Nick's original goal.
- document the regression recalc
- add regression from correlation matrix recalc
- add mediation analysis recalc - eg Rutherford et al 2017
- **Integrate granularity and bounds tests more generally.** GRIM / GRIMMER-style granularity checks (does a reported mean / SD have a denominator that a stated N can produce?) and physical-bounds checks (is the reported mean within the scale's logical range; is the reported SD within the Bhatia-Davis bound implied by the scale endpoints and the reported mean / N) are currently only wired into `recalc_missing_subgroup()` (E2). They apply equally to E1's recalculated quantities, to `independent_t_test_summary()`'s reported group statistics, and to any other check whose inputs include means and SDs on a bounded scale. Factor `scale_min` / `scale_max` arguments, a shared Bhatia-Davis bound helper, and a granularity helper into a common layer used by every reported-descriptive-statistic check.
- Validate `recalc_sd_concentration()` (E3): the lower-tail Bartlett interpretation as a TIVA-style forensic signal is currently flagged as exploratory in its docstring and vignette. Run null-calibration sims at realistic (k, n_g) and under mild non-normality; characterize power against plausible fabricated-SD patterns; only then graduate it from "diagnostic" to "check" and consider inclusion in any `audit_*()` wrapper.
