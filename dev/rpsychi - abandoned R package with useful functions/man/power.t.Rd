\encoding{UTF-8}
\name{power.t}
\alias{power.t}

\title{Statistical power for the t test for means}
\description{
    This is an internal function. \code{power.t} computes statistical power for the t test for means.
}

\usage{
power.t(n1, n2, delta, sig.level = 0.05)
}

\arguments{
  \item{n1}{a numeric contains the sample size in group 1}
  \item{n2}{a numeric contains the sample size in group 2}
  \item{delta}{Cohen's \eqn{d}}
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
##Cohen (1988) ex.2.1
power.t(sig.level=.05, delta=0.5, n1=30, n2=30)
}

\keyword{internal}


