\encoding{UTF-8}

\name{power.prop}
\alias{power.prop}

\title{Statistical power for differences between proportions}

\description{
    This is an internal function. \code{power.prop} computes statistical power for differences between proportions.
}
\usage{
power.prop(n, h = NULL, sig.level = 0.05, one.sided = FALSE, prop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a numeric vector (length(\code{n}) = 2) contains the sample size in each group}
  \item{h}{Cohen's \eqn{h}}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{one.sided}{nondirectional (\code{one.sided} = FALSE, default) or directional (\code{one.sided} = TRUE)}
  \item{prop}{...}
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
##Cohen (1988) ex.6.2
power.prop(h=.4, n=c(100,100), sig.level=0.05)
}

\keyword{internal}


