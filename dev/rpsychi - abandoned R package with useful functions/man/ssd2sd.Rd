\encoding{UTF-8}
\name{ssd2sd}
\alias{ssd2sd}


\title{Convert sample standard deviation into unbiased one}

\description{
    \code{ssd2sd} converts a sample standard deviation into unbiased one.
}
\usage{
ssd2sd(n, ssd)
}

\arguments{
  \item{n}{a numeric contains the sample size}
  \item{ssd}{a numeric contains the sample standard deviation}
}



\value{ 
  Return a numeric containing the unbiased standard deviation.
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 



\seealso{\code{\link{r2cov}}, \code{\link{svar}}, \code{\link{ssd}}, 
        \code{\link{svar}}, 
        \code{\link{lower2R}}
} 


\examples{
data(infert)
ssd2sd(nrow(infert), ssd(infert$age))
sd(infert$age)
}

\keyword{univar}

