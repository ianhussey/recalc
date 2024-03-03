\encoding{UTF-8}
\name{samplesize.etasq}
\alias{samplesize.etasq}


\title{Sample size estimation in the analysis of variance}
\description{
 \code{samplesize.etasq} determines the sample size, when specified the number of group, effect size to be detected, desired statistical power, and significance level. 
}

\usage{
samplesize.etasq(k, delta, power=.80, sig.level=.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{k}{a numeric contains the number of group}
  \item{delta}{a numeric contains the effect size (\eqn{\eta^2}) to be detected}
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


\seealso{\code{\link{ind.oneway}}, \code{\link{ind.oneway.second}}}

\examples{
##Cohen (1988) ex.8.10
f <- .25
samplesize.etasq(k=4, delta= f^2/(1+f^2), power=.80, sig.level=.05)
}

\keyword{design}
\keyword{htest}

