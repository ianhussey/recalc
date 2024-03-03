\encoding{UTF-8}
\name{FNONCT}
\alias{FNONCT}

\title{Compute the noncentrality parameter of an F distribution}



\description{
    This is an internal function. \code{FNONCT} compute the value of the noncentrality parameter of an F distribution. 
}
\usage{
FNONCT(x, df1, df2, prob, 
        interval = c(0, 10000), my.tol = 1e-06)
}


\arguments{
  \item{x}{a numeric contains the \eqn{F} value}
  \item{df1}{a numeric contains the numerator degrees of freedom parameter}
  \item{df2}{a numeric contains the denominator degrees of freedom parameter}
  \item{prob}{a numeric contains the probability}
  \item{interval}{a numeric vector contains the end-points of the interval to be searched}
  \item{my.tol}{the desired accuracy}
}


\value{
  Return a numeric contains the noncentrality parameter.
}

\references{Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp}
        } 



\examples{
##Kline (2004) Table 4.7
FNONCT(9.61, 1, 58, .975)
FNONCT(9.61, 1, 58, .025)
}

\keyword{internal}
