# recalc



# Built upon / inspired by

rpsychi R package (v0.8, no longer on CRAN due to lack of maintenance)

Nick Brown's f_range() function (https://steamtraen.blogspot.com/2018/01/)



# TODO

- Add option to t-test function to check if the result recalcultes if we assume the SD is actually an SE.
- Add way to directly compare multple among t/df/p, cohens d, M/SD/N?  

# INTERNAL NOTES

## Ideas for scope

### Within scope

Recalculation of statistical results from reported results.

Calculating discrete values or intervals of possible values of a given statisics (e.g., test stat, p value, effect size) given different assumptions about the reported values (e.g., summary stats). 

- For example, given reported M and SD for both experimental and control conditons, what range of cohen's d values are consistent with this result, given (a) different rounding, (b) method of calculating cohen's d.

### Out of scope

Some of the more inferential tests about the flagging (possible) errors.



## Ideas for features

- Consider rounding in reported results.



## Resources

- different versions and implementations of within-subjects cohen's d: https://github.com/ianhussey/versions-of-cohens-d
- different versions and implementtions of between-groups cohen's d: https://rpubs.com/metinbulus/welch
- rpsychi package has some useful functions. code has been copied into the `/dev` directory
- nick brown has an ANOVA from summary stats function that accounts for rounding, code in `/dev`
- Extracting pre-post correlations from summary stats: https://matthewbjane.quarto.pub/pre-post-correlations/
- Lisa Debruine's {faux} package has some error detection functions https://github.com/debruine/faux
- Lisa Debruine's {within} repo assesses the plausibility of a within subject t-test by showing the range of possible between timepoint correlations given the reported results https://github.com/debruine/within

# README DRAFT

The goal of *recalc* is to facilitate error checking in research articles. To this end, it offers functions that recalculate statistical tests from reported summary statistics. In addition, it offers functions that can help to assess the plausibility of effect sizes, and significance starts (e.g., in correlation tables).

## Installation

Currently, recalc is only available on GitHub. To install it, run the following code in your R Console:

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("ianhussey/recalc")
```



