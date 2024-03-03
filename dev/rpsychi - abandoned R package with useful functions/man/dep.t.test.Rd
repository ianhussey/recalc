\encoding{UTF-8}
\name{dep.t.test}
\alias{dep.t.test}

\title{A t-test with dependent samples using individual data: Reporting effect size}

\description{
  \code{dep.t.test} conducts a t-test with dependent samples using individual data. 
}


\usage{
    dep.t.test(formula, data, block, 
                    sig.level=.05, digits=3)
}

\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of one independent variable containing a factor with two levels}
  \item{data}{a data frame contains the variables in the \code{fomrmula}}
  \item{block}{a character string specify the blocking variable}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function conducts a t-test with dependent samples using individual data. 
}

\value{
  The returned object of \code{dep.t.test.second} contains the following components:
  
  \item{samp.stat}{returns the means, standard deviations, sample size, and correlation}
  \item{raw.difference}{returns a raw mean difference, its' confidence interval, and standard error}
  \item{standardized.difference}{returns a standardized mean difference (Hedges's \eqn{g}) and its' approximate confidence interval for a population standardized mean difference}
}

\references{Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{dep.t.test.second}}}

\examples{ 
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=2)
                  )
dep.t.test(y~x, block="subj", data=dat)

}
\keyword{design}
\keyword{htest}
