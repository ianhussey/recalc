\encoding{UTF-8}

\name{power.multi} 
 
\alias{power.multi} 


\title{Statistical power for F tests of variance proportions}


\description{
    This is an internal function. \code{power.multi} computes statistical power for F tests of variance proportions.
}
\usage{
power.multi(n, n.ind, delta, sig.level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{n.ind}{number of independent variables}
  \item{delta}{Cohen's \eqn{f^2}}
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
##Cohen (1988) ex 9.1
power.multi(n=95, n.ind=5, delta=.1111, sig.level=.05)

}

\keyword{internal}

