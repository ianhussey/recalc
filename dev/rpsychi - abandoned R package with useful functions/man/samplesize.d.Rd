\encoding{UTF-8}
\name{samplesize.d}
\alias{samplesize.d}

\title{Sample size estimation in the t test for means}

\description{
 \code{samplesize.d} determines the sample size, when specified the effect size to be detected, desired statistical power, and significance level. 
 \code{samplesize.d} always specifies two-sided test.
}
\usage{
samplesize.d(delta, power = 0.8, sig.level = 0.05)
}


\arguments{
  \item{delta}{a numeric contains the effect size (Cohen's \eqn{d}) to be detected}
  \item{power}{a numeric contains the desired statistical power (default 0.80)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
}



\value{
    Return a numeric containing the appropriate sample size in each group.
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

\seealso{\code{\link{ind.t.test}}, \code{\link{ind.t.test.second}}}

\examples{
##Cohen (1988) ex.2.9
samplesize.d(delta=.20, power=.95, sig.level=.05)
}

\keyword{design}
\keyword{htest}
