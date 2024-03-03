\encoding{UTF-8}
\name{lower2R} 
\alias{lower2R}

\title{Convert a vector containing correlations into a correlation matrix} 

\description{ 
    
    \code{lower2R} converts a vector containing correlations into a correlation matrix.
}

\usage{
lower2R(x, varname=NULL)
}


\arguments{
  \item{x}{a numeric vector contains the correlations}
  \item{varname}{a character vector contains the row and column names}
}


\details{
This function converts a vector containing correlations into a correlation matrix. 
\code{x} contains the elements that consist of the lower triangle of the correlation matrix. 
The length of \code{varname} should be equal to the number of rows and columns in the correlation matrix.
}

\value{
  Return a matrix containing the correlation matrix.
}



\references{Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.}



\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{svar}}, \code{\link{ssd}}, 
        \code{\link{svar}}, \code{\link{ssd2sd}}, 
        \code{\link{r2cov}}
}



\examples{
lower2R(c(1:15))

##Kline (2004) Table 6.15
lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
        .84,.8,.91,.88,.93))

lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
        .84,.8,.91,.88,.93), 
        varname=paste("trial", 1:6, sep=""))
}

\keyword{array}

