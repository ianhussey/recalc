\encoding{UTF-8}
\name{power.f2}
\alias{power.f2}
\title{
Statistical power for F tests on means in the two-way analysis of variance}
\description{
    This is an internal function. \code{power.f2} computes statistical power for F tests on means in the two-way analysis of variance.
}
\usage{
power.f2(df1, df2, delta, sig.level = 0.05)
}

\arguments{
  \item{df1}{degree of freedom of the numerator of the F ratio}
  \item{df2}{degree of freedom of the denominator of the F ratio}
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
##Cohen (1988) ex 8.6
power.f2(df1=1,df2=96,sig.level=0.05,delta=.10)
power.f2(df1=2,df2=96,sig.level=0.05,delta=.25)
power.f2(df1=3,df2=96,sig.level=0.05,delta=.40)

power.f2(df1=1,df2=120,sig.level=0.05,delta=.10)
power.f2(df1=2,df2=120,sig.level=0.05,delta=.25)
power.f2(df1=3,df2=120,sig.level=0.05,delta=.40)
}

\keyword{internal}

