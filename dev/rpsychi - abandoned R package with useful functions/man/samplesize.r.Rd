\encoding{UTF-8}
\name{samplesize.r}
\alias{samplesize.r}


\title{Sample size estimation in the significance testing of a product moment correlation}
\description{
 \code{samplesize.r} determines the sample size, when specified the effect size to be detected, desired statistical power, and significance level. 
 \code{samplesize.r} always specifies two-sided test.
}
\usage{
samplesize.r(delta, power = 0.8, sig.level = 0.05)
}


\arguments{
  \item{delta}{a numeric contains the effect size (product moment correlation) to be detected}
  \item{power}{a numeric contains the desired statistical power (default 0.80)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
}

\value{
  Return a numeric containing the appropriate sample size.
}


\references{
Cohen J (1988). Statistical power analysis for the behavioral sciences (2nd ed). Hillsdale, NJ: Erlbaum.
}


\author{Yasuyuki Okumura\cr
        Department of Social Psychiatry, \cr
        National Institute of Mental Health, \cr
        National Center of Neurology and Psychiatry \cr
        \email{yokumura@blue.zero.jp}
        }

\seealso{\code{\link{zero.r.test}}, \code{\link{zero.r.test.second}}}

\examples{
##Cohen (1988) ex.3.4
samplesize.r(delta=.3, power=.8, sig.level=.05)
}

\keyword{design}
\keyword{htest}

