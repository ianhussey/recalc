\encoding{UTF-8}
\name{rpsychi-package}
\alias{rpsychi-package}
\alias{rpsychi}
\docType{package}
\title{
Statistics for psychiatric research
}
\description{
The rpsychi offers a number of functions for psychiatry, psychiatric nursing, clinical psychology. Functions are primarily for statistical significance testing using published work. For example, you can conduct a factorial analysis of variance (ANOVA), which requires only the mean, standard deviation, and sample size for each cell, rather than the individual data. This package covers fundamental statistical tests such as t-test, chi-square test, analysis of variance, and multiple regression analysis. With some exceptions, you can obtain effect size and its confidence interval. These functions help you to obtain effect size from published work, and then to conduct a priori power analysis or meta-analysis, even if a researcher do not report effect size in a published work.
}

\details{
\tabular{ll}{
Package: \tab rpsychi\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2011-11-12\cr
License: \tab GPL (version 2 or later)\cr
LazyLoad: \tab yes\cr
}
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 

\references{
Cohen B (2000) Calculating a factorial ANOVA from means and standard deviations. Understanding Statistics, 1, 191-203.

Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Cohen J (1988). Statistical power analysis for the behavioral sciences (2nd ed). Hillsdale, NJ: Erlbaum.

Cohen J, Cohen P, Aiken LS (2003) Applied multiple regression/correlation analysis for the behavioral sciences (3rd ed). Mahwah, NJ: Erlbaum.

Iwasaki M (2002) Foundations of incomplete data analysis (in Japanese) Tokyo: EconomistSha Publishing.

Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.

Smithson M (2001) Correct confidence intervals for various regression effect sizes and parameters: The importance of noncentral distributions in computing intervals, 61, 605-632.

Tabachnick BG, Fidell LS (2007) Experimental designs using ANOVA. Belmont, CA: Thomson. 

Toyoda H (1998) Introduction to structural equation modeling (in Japanese) Tokyo: Asakura Publishing.
}


\keyword{ package }


\examples{
##Kline (2004) Table 6.15
my.cont <- matrix(c(-5,-3,-1,1,3,5,
                   5,-1,-4,-4,-1,5), ncol=6, nrow=2, byrow=TRUE)
dep.oneway.second(m = c(11.77,21.39,27.5,31.02,32.58,34.2), 
                  sd = c(7.6,8.44,8.95,9.21,9.49,9.62), 
                  n = 137, 
                  corr=lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
                  .84,.8,.91,.88,.93)),
                  contr=my.cont)
}
