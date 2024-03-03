\encoding{UTF-8}
\name{samplesize.h}
\alias{samplesize.h}

\title{Sample size estimation in the differences between proportions}

\description{
 \code{samplesize.h} determines the sample size, when specified the effect size to be detected, desired statistical power, and significance level. 
 \code{samplesize.h} always specifies two-sided test.
}
\usage{
samplesize.h(delta, power = 0.8, sig.level = 0.05)
}


\arguments{
  \item{delta}{a numeric contains the effect size (Cohen's \eqn{h}) to be detected}
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

\seealso{\code{\link{ind.prop}}, \code{\link{ind.prop.second}}}

\examples{
##Cohen (1988) ex.6.7
samplesize.h(delta=.20, power=.90, sig.level=.01)

##Cohen (1988) ex.6.8
samplesize.h(delta=.2828, power=.95, sig.level=.05)
}

\keyword{design}
\keyword{htest}
