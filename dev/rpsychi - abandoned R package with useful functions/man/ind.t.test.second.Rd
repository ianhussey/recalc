\encoding{UTF-8}
\name{ind.t.test.second}
\alias{ind.t.test.second}

\title{A t-test with independent samples using published work: Reporting effect size}
\description{
   \code{ind.t.test.second} conducts a t-test with independent samples using published work. 
}
\usage{
ind.t.test.second(m, sd, n, 
        unbiased = TRUE, correct=TRUE, sig.level = 0.05, digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{a numeric vector contains the means (length(\code{m}) = 2)}
  \item{sd}{a numeric vector contains the sample/unbiased standard deviations (length(\code{sd}) = 2)}
  \item{n}{a numeric contains the sample size (length(\code{n}) = 2)}
  \item{unbiased}{\code{sd} contains unbiased standard deviations (\code{unbiased} = TRUE, default) or sample standard deviations (\code{unbiased} = FALSE)}
  \item{correct}{a logical indicating whether to compute an unbiased standardized mean difference (\eqn{delta}) or not (\code{correct} = TRUE)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}
\details{
 This function conducts a t-test with independent samples using published work. 
Statistical power is calculated using the following specifications: 

(a) small (\eqn{d = 0.20}), medium (\eqn{d = 0.50}), and large (\eqn{d = 0.80}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{n}

(c) significance level specified by \code{sig.level}
}

\value{
  The returned object of \code{ind.t.test.second} contains the following components:
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

\seealso{\code{\link{ind.t.test}}, \code{\link{samplesize.d}}}




\examples{
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5)
                  )
ind.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = tapply(dat$y, dat$x, length), correct=FALSE
                  )
ind.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = tapply(dat$y, dat$x, length), correct=TRUE
                  )     #approximate unbiased estimator of delta
}

\keyword{design}
\keyword{htest}

