\encoding{UTF-8}
\name{formatted}
\alias{formatted}

\title{
  Convert a numeric vector into a character vector with the specified number of decimal place
}
\description{
  \code{formatted} convert a numeric vector into a character vector with the specified number of decimal place. 
} 
\usage{
formatted(x, digits = 2)
}

\arguments{
  \item{x}{a numeric vector}
  \item{digits}{the specified number of decimal places (default 2)}
}

\value{
  Return a character vector.
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\examples{
data(infert)
x <- svar(infert$age)    #sample variance
formatted(x)
formatted(x, digits=4)
}

\keyword{univar}
