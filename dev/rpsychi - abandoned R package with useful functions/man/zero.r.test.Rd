\encoding{UTF-8}
\name{zero.r.test}
\alias{zero.r.test}
\title{A significance testing of a product moment correlation using individual data}
\description{
    \code{zero.r.test} conducts a significance testing of a product moment correlation using individual data
}
\usage{
zero.r.test(formula, data, sig.level = 0.05, digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of one independent variable containing a numeric variable}
  \item{data}{a data frame contains the variables in the \code{fomrmula}}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function conducts a significance testing of a product moment correlation using individual data.
Statistical power is calculated using the following specifications: 

(a) small (\eqn{r = 0.10}), medium (\eqn{r = 0.30}), and large (\eqn{r = 0.50}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{data}

(c) significance level specified by \code{sig.level}

}

\value{
  The returned object of \code{zero.r.test} contains the following components:
  \item{samp.stat}{returns the means and unbiased standard deviations}

  \item{correlation}{returns a product moment correlation, its' approximate confidence interval for population correlation, and standard error}
  \item{power}{returns statistical power for detecting 
small (\eqn{r = 0.10}), medium (\eqn{r = 0.30}), and large (\eqn{r = 0.50}) population effect sizes}
}


\references{
Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.
}




\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{zero.r.test.second}}, \code{\link{samplesize.r}}}

\examples{
dat <- data.frame(x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1),
                  y = c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8))          
zero.r.test(y~x, data=dat)
}
\keyword{design}
\keyword{htest}

