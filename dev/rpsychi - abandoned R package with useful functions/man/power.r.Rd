\encoding{UTF-8}
\name{power.r}
\alias{power.r}

\title{Statistical power for the significance testing of a product moment correlation}
\description{
    This is an internal function. \code{power.r} computes power for the significance of a product moment correlation.

}

\usage{
power.r(n, delta, sig.level = 0.05)
}


\arguments{
  \item{n}{sample size}
  \item{delta}{population correlation}
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
##Cohen (1988) ex.3.1
power.r(n=50, delta=0.3, sig.level=0.05)
}


\keyword{internal}

