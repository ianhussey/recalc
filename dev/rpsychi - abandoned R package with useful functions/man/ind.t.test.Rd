\encoding{UTF-8}
\name{ind.t.test}
\alias{ind.t.test}

\title{A t-test with independent samples using individual data: Reporting effect size}

\description{
   \code{ind.t.test} conducts a t-test with independent samples using individual data. 
}
\usage{
ind.t.test(formula, data, correct=TRUE,
                sig.level = 0.05, digits = 3)
}

\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of one independent variable containing a factor with two levels}
  \item{data}{a data frame contains the variables in the \code{fomrmula}}
  \item{correct}{a logical indicating whether to compute an unbiased standardized mean difference (\eqn{delta}) or not (\code{correct} = TRUE)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}
\details{
This function conducts a t-test with independent samples using individual data. 
Statistical power is calculated using the following specifications: 

(a) small (\eqn{d = 0.20}), medium (\eqn{d = 0.50}), and large (\eqn{d = 0.80}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{formula} and \code{data}

(c) significance level specified by \code{sig.level}
 
}

\value{
  The returned object of \code{ind.t.test} contains the following components:
  \item{samp.stat}{returns the means, standard deviations, and sample sizes}
  \item{raw.difference}{returns a raw mean difference, its' confidence interval, and standard error}
  \item{standardized.difference}{returns a standardized mean difference (Hedges's \eqn{g}), its' approximate confidence interval for population standardized mean difference, and standard error}
  \item{power}{returns statistical power for detecting 
  small (\eqn{d = 0.20}), medium (\eqn{d = 0.50}), and large (\eqn{d = 0.80}) population effect sizes}
}

\references{
Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.
}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 

\seealso{\code{\link{ind.t.test.second}}, \code{\link{samplesize.d}}}



\examples{
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5)
                  )
ind.t.test(y~x, data=dat, correct=FALSE)
}

\keyword{design}
\keyword{htest}

