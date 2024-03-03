\encoding{UTF-8}
\name{svar}
\alias{svar}


\title{Compute sample variance}
\description{
    \code{svar} computes a sample variance. 
}
\usage{
svar(x, na.rm = TRUE)
}

\arguments{
  \item{x}{a numeric vector}
  \item{na.rm}{missing values (NA) will be stripped before the computation proceeds (\code{na.rm} = TRUE, default)}
}

\value{
  Return a numeric containing a sample variance.
}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{r2cov}}, \code{\link{svar}}, \code{\link{ssd}}, 
        \code{\link{ssd2sd}}, 
        \code{\link{lower2R}}
} 


\examples{
data(infert)
svar(infert$age)    #sample variance
var(infert$age)     #unbiased variance
}
\keyword{univar}
