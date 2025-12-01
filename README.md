

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
