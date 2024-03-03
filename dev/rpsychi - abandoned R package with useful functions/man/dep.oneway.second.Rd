\encoding{UTF-8}
\name{dep.oneway.second} 
\alias{dep.oneway.second}
\title{A one-way design with dependent samples using published work: Reporting effect size} 
\description{ 
    \code{dep.oneway.second} conducts a one-way design with dependent samples, namely one-way repeated-measures analysis of variance, using published work. 
}

\usage{
dep.oneway.second(m, sd, n, corr, 
            unbiased = TRUE, contr = NULL, sig.level = 0.05, digits = 3)
}


\arguments{
  \item{m}{a numeric vector contains the means (length(\code{m}) >= 2)}
  \item{sd}{a numeric vector contains the sample/unbiased standard deviations (length(\code{sd}) >= 2)}
  \item{n}{a numeric contains the sample size (length(\code{n}) >= 2)}
  \item{corr}{a matrix or data frame contains the correlation matrix}
  \item{unbiased}{\code{sd} contains unbiased standard deviations (\code{unbiased} = TRUE, default) or sample standard deviations (\code{unbiased} = FALSE)}
  \item{contr}{a matrix or vector contains the contrast weights}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function conducts a one-way design with dependent samples, namely one-way repeated-measures analysis of variance, using published work.
If you do not specify \code{contr}, all possible pairwise contrasts will be calculated.
}

\value{
  The returned object of \code{dep.oneway.second} contains the following components:
  \item{anova.table}{returns a ANOVA table containing sums of squares, degrees of freedom, mean squares, \eqn{F} values, and a \eqn{p} value}
  \item{omnibus.es}{returns a omnibus effect size which is a partial \eqn{\eta^2}}
  \item{raw.contrasts}{returns raw mean differences, their confidence intervals, and standard errors}
  \item{standardized.contrasts}{returns standardized mean differences for the contrasts (Hedges's \eqn{g}) and their approximate confidence intervals for population standardized mean differences}
}

\references{Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 

\seealso{\code{\link{dep.oneway}}, \code{\link{lower2R}}}
 
\examples{ 
##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=3)
                  )
dep.oneway(formula = y~x, data=dat, block="subj")

datwide <- reshape(dat, direction="wide", idvar="subj", timevar="x")
tmp <- datwide[,-1]
dep.oneway.second(m = apply(tmp, 2, mean), apply(tmp, 2, sd), n = nrow(tmp), corr=cor(tmp))



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
\keyword{design}
\keyword{htest}
