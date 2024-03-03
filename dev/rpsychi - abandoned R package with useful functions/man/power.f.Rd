\encoding{UTF-8}
\name{power.f}
\alias{power.f}

\title{Statistical power for F tests on means in the analysis of variance}

\description{
    This is an internal function. \code{power.f} computes statistical power for F tests on means in the analysis of variance.
}
\usage{
power.f(u, n, delta, sig.level = 0.05)
}

\arguments{
  \item{u}{degree of freedom of the numerator of the F ratio}
  \item{n}{sample size}
  \item{delta}{Cohen's \eqn{f}}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
}

\value{
    Return a numeric containing the statistical power.
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




\examples{
##Cohen (1988) ex.8.1
power.f(u = 3, n = 20, delta=.28, sig.level=.05)
power.f(u = 3, n = 20, delta=.28, sig.level=.10)


##Cohen (1988) ex 8.2
power.f(u = 2, n = 200, delta=.23, sig.level=.01)
power.f(u = 2, n = 200, delta=.33, sig.level=.01)
}

\keyword{internal}
