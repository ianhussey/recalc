\encoding{UTF-8}
\name{samplesize.rsq}
\alias{samplesize.rsq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Sample size estimation in the F tests of variance proportions}

\description{
 \code{samplesize.rsq} determines the sample size, when specified the effect size to be detected, desired statistical power, and significance level. 
 
}
\usage{
samplesize.rsq(delta, n.ind, power = 0.8, sig.level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{a numeric contains the effect size (Cohen's \eqn{R^2}) to be detected}
  \item{n.ind}{a numeric contains the number of independent variables}
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


\seealso{\code{\link{multreg}}, \code{\link{multreg.second}}}

\examples{
##Cohen (1988) ex.9.18
samplesize.rsq(delta=.16, n.ind=20, power = .90, sig.level=.05)
}

\keyword{design}
\keyword{htest}

