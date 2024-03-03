\encoding{UTF-8}
\name{ind.prop.second}
\alias{ind.prop.second}

\title{A Z test for the equality of two proportions using published workc}
\description{
    \code{ind.prop.second} conducts a Z test for the equality of two proportions using published work. 
}
\usage{
ind.prop.second(x, n, sig.level = 0.05, digits = 3, ref.ind=1)
}

\arguments{
  \item{x}{a numeric vector (length(\code{x}) = 2) contains the counts of successes}
  \item{n}{a numeric vector (length(\code{n}) = 2) contains the sample sizes}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
  \item{ref.ind}{the reference for the independent variable (default 1)}
}
\details{
This function conducts a Z test for the equality of two proportions using published work. 
Statistical power is calculated using the following specifications: 

(a) small (\eqn{h = 0.20}), medium (\eqn{h = 0.50}), and large (\eqn{h = 0.80}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{n}

(c) significance level specified by \code{sig.level}

}

\value{
  The returned object of \code{ind.prop.second} contains the following components:
  \item{samp.stat}{returns the proportions, sample sizes, and a effect size index (Cohen's \eqn{h})}
  \item{risk.difference}{returns a risk difference, its' confidence interval, and standard error}
  \item{risk.ratio}{returns a risk ratio, its' confidence interval, and a standard error of a log-transformed risk ratio}
  \item{odds.ratio}{returns a odds ratio, its' confidence interval, and a standard error of a log-transformed odds ratio}
  \item{power}{returns statistical power for detecting 
  small (\eqn{h = 0.20}), medium (\eqn{h = 0.50}), and large (\eqn{h = 0.80}) population effect sizes}
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

\seealso{\code{\link{ind.prop}}, \code{\link{samplesize.h}}}

\examples{
##Kline (2004) Chapter 5
x1 <- c("relapsed", "not relapsed")
y1 <- c("control", "treatment")

dat <- data.frame(y =         
factor(c(rep(x1, c(60, 40)), rep(x1, c(40, 60))), levels=x1),
x = factor(rep(y1, each=100), levels=y1)
)
tab <- xtabs(~x+y, data=dat)
tab
ind.prop.second(x=tab[,1], n = rowSums(tab))             #Risk for relapse is lower in treatment than control condition.
ind.prop.second(x=tab[,1], n = rowSums(tab), ref.ind=2)  #Risk for relapse is higher in control than treatment condition.
}

\keyword{design}
\keyword{htest}

