\encoding{UTF-8}
\name{ssd}
\alias{ssd}

\title{Compute sample standard deviation}
\description{
    \code{ssd} computes a sample standard deviation. 
}

\usage{
ssd(x, na.rm = TRUE)
}


\arguments{
  \item{x}{a numeric vector}
  \item{na.rm}{missing values (NA) will be stripped before the computation proceeds (\code{na.rm} = TRUE, default)}
}

\value{
  Return a numeric containing a sample standard deviation.
}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{r2cov}}, 
        \code{\link{svar}}, \code{\link{svar}}, \code{\link{ssd2sd}}, 
        \code{\link{lower2R}}
}

\examples{
data(infert)
ssd(infert$age)     #sample standard deviation
sd(infert$age)      #unbiased standard deviation
}
\keyword{univar}

