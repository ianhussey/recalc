\encoding{UTF-8}
\name{dep.oneway}
\alias{dep.oneway}

\title{A one-way design with dependent samples using individual data: Reporting effect size}

\description{
    \code{dep.oneway} conducts a one-way design with dependent samples, namely one-way repeated-measures analysis of variance, using individual data. 
}


\usage{
dep.oneway(formula, data, block, 
            contr = NULL, sig.level = 0.05, digits = 3)
}

\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of one independent variable containing a factor with two or more levels}
  \item{data}{a data frame contains the variables in the \code{fomrmula} and \code{block}}
  \item{block}{a character string specify the blocking variable}
  \item{contr}{a matrix or vector contains the contrast weights}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}
\details{
This function conducts a one-way design with dependent samples, namely one-way repeated-measures analysis of variance, using individual data.
If you do not specify \code{contr}, all possible pairwise contrasts will be calculated.
}
\value{
  The returned object of \code{dep.oneway} contains the following components:
  \item{anova.table}{returns a ANOVA table containing sums of squares, degrees of freedom, mean squares, \eqn{F} values, and a \eqn{p} value}
  \item{omnibus.es}{returns a omnibus effect size which is a partial \eqn{\eta^2}}
  \item{raw.contrasts}{returns raw mean differences, their confidence intervals, and standard errors}
  \item{standardized.contrasts}{returns standardized mean differences for the contrasts (Hedges's \eqn{g}) and their approximate confidence intervals for population standardized mean differences}
}

\references{
Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association. 
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 



\seealso{\code{\link{dep.oneway.second}}}

\examples{
##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=3)
                  )
dep.oneway(formula = y~x, data=dat, block="subj")


##contrast 1: a - c, contrast 2 : 1/2(a + c) - b
my.cont <- matrix(c(1,0,-1,1/2,-1,1/2), ncol=3, nrow=2, byrow=TRUE)
dep.oneway(formula = y~x, data=dat, block="subj", contr=my.cont)
}
\keyword{design}
\keyword{htest}

