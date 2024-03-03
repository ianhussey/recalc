\encoding{UTF-8}
\name{dep.t.test.second}
\alias{dep.t.test.second}

\title{A t-test with dependent samples using published work: Reporting effect size}

\description{
  \code{dep.t.test.second} conducts a t-test with dependent samples using published work. 
}


\usage{
    dep.t.test.second(m, sd, n, corr, 
            unbiased = TRUE, sig.level = 0.05, digits = 3)
}

\arguments{
  \item{m}{a numeric vector contains the means (length(\code{m}) = 2)}
  \item{sd}{a numeric vector contains the sample/unbiased standard deviations (length(\code{sd}) = 2)}
  \item{n}{a numeric contains the sample size (length(\code{n}) = 2)}
  \item{corr}{a numeric contains the correlation}
  \item{unbiased}{\code{sd} contains unbiased standard deviations (\code{unbiased} = TRUE, default) or sample standard deviations (\code{unbiased} = FALSE)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}
\details{
This function conducts a t-test with dependent samples using published work. 
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


\seealso{\code{\link{dep.t.test}}}

\examples{ 
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=2)
                  )
datwide <- reshape(dat, direction="wide", idvar="subj", timevar="x")

dep.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = nlevels(dat$subj),
                  corr = cor(datwide[,2:3])[1,2]
                  )

dep.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = 30,
                  corr = cor(datwide[,2:3])[1,2]
                  )
}
\keyword{design}
\keyword{htest}
