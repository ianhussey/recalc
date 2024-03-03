\encoding{UTF-8}
\name{r2cov} 

\alias{r2cov} 


\title{Convert correlation matrix into covariance matrix}


\description{ 
    \code{r2cov} converts correlation matrix and 
    sample/unbiased standard deviation into sample/unbiased 
    covariance matrix. 
} 

\usage{ r2cov(sd, R) }  

\arguments{ 
  \item{sd}{a numeric vector contains the sample/unbiased standard deviations}
  \item{R}{a matrix or data frame contains the correlation matrix}
} 

\details{This function converts correlation matrix and 
    sample/unbiased standard deviation into sample/unbiased 
    covariance matrix using the following equation: \eqn{S = D^{1/2} R 
    D^{1/2}}, where \eqn{S} is a sample/unbiased covariance matrix, 
    \eqn{R} is a correlation matrix, and \eqn{D^{1/2}} is a square 
    matrix with \code{sd} on the main diagonal and 0's elsewhere. 
    The length of \code{sd} should be equal to the number of rows 
    and columns in \code{R}.
} 

\value{ 
  Return a matrix containing the sample/unbiased covariance matrix.
}


\references{Toyoda H (1998) Introduction to structural equation modeling (in Japanese) Tokyo: Asakura Publishing.}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 

\seealso{\code{\link{svar}}, \code{\link{ssd}}, 
        \code{\link{svar}}, \code{\link{ssd2sd}}, 
        \code{\link{lower2R}}
} 

\examples{ 
##data(iris) 
x <- iris[,1:4] 
cov(x)
r2cov(apply(x, 2, sd), cor(x)) 


##Toyoda (1998) p.34 
r2cov(sd = sqrt(c(.862, 1.089, 0.606)), 
      R = lower2R(c(.505, -0.077, -.233)))
} 

\keyword{array}
