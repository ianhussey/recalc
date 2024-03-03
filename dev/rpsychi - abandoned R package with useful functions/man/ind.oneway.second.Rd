\encoding{UTF-8}
\name{ind.oneway.second}
\alias{ind.oneway.second}

\title{A one-way design with independent samples using published work: Reporting effect size}
\description{
    \code{ind.oneway.second} conducts a one-way design with independent samples, namely one-way randomized-group analysis of variance, using published work. 
}
\usage{
ind.oneway.second(m, sd, n, 
            unbiased = TRUE, contr = NULL, sig.level = 0.05, digits = 3)
}


\arguments{
  \item{m}{a numeric vector contains the means (length(\code{m}) >= 2)}
  \item{sd}{a numeric vector contains the sample/unbiased standard deviations (length(\code{sd}) >= 2)}
  \item{n}{a numeric contains the sample size (length(\code{n}) >= 2)}
  \item{unbiased}{\code{sd} contains unbiased standard deviations (\code{unbiased} = TRUE, default) or sample standard deviations (\code{unbiased} = FALSE)}
  \item{contr}{a matrix or vector contains the contrast weights}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function conducts a one-way design with independent samples, namely one-way randomized-group analysis of variance, using published work.
If you do not specify \code{contr}, all possible pairwise contrasts will be calculated. 
Statistical power is calculated using the following specifications: 

(a) small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), and large (\eqn{\eta^2 = 0.14}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{n}

(c) significance level specified by \code{sig.level}

}

\value{
  The returned object of \code{ind.oneway.second} contains the following components:
  \item{anova.table}{returns a ANOVA table containing sums of squares, degrees of freedom, mean squares, \eqn{F} values}
  \item{omnibus.es}{returns a omnibus effect size which is a \eqn{\eta^2}, and its' confidence interval}
  \item{raw.contrasts}{returns raw mean differences, their confidence intervals, and standard errors}
  \item{standardized.contrasts}{returns standardized mean differences for the contrasts (Hedges's \eqn{g}), their approximate confidence intervals for population standardized mean differences, and standard errors}
  \item{power}{returns statistical power for detecting small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), 
  and large (\eqn{\eta^2 = 0.14}) population effect sizes}
}

\references{
Cohen B (2000) Calculating a factorial ANOVA from means and standard deviations. Understanding Statistics, 1, 191-203.

Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{ind.oneway}}, \code{\link{samplesize.etasq}}}


\examples{
##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5)
                  )                 

##contrast 1: a - c, contrast 2: 1/2(a + c) - b
my.cont <- matrix(c(1,0,-1,1/2,-1,1/2), ncol=3, nrow=2, byrow=TRUE)


ind.oneway.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n= tapply(dat$y, dat$x, length)) 

ind.oneway.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n= tapply(dat$y, dat$x, length),
                  contr = my.cont)   

}

\keyword{design}
\keyword{htest}


