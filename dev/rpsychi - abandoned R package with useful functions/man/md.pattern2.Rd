\encoding{UTF-8}
\name{md.pattern2}
\alias{md.pattern2}

\title{Display missing-data patterns}
\description{
  \code{md.pattern2} displays missing-data patterns
}
\usage{
md.pattern2(x)
}

\arguments{
  \item{x}{a data frame}
}

\details{
This function displays missing-data patterns. \code{md.pattern} in library(mice) is extended to sort variables.
}

\value{
 Return a matrix with \code{ncol(x) + 1} columns, in which each row corresponds
     to a missing data pattern (1 = observed, 0 = missing).  Rows and
     columns are sorted in increasing amounts of the sample sizes.
     The last column and row contain row and column counts, respectively.
     The row name contains the sample sizes in the specific missing data patterns.
}

\references{Iwasaki M (2002) Foundations of incomplete data analysis (in Japanese) Tokyo: EconomistSha Publishing.}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 



\examples{
##Iwasaki (2002)
dat <- data.frame(matrix(c(
71  ,  68,  72,  72,  90,  72,  77,  76,  84,  77,
1850,2000,2100,1700,  NA,2200,2150,  NA,  NA,  NA,
136 , 139, 147, 142,  NA, 150, 156,  NA, 152,  NA,
34  , 45 ,  50,  38,  NA,  41,  43,  52,  57,  48
), ncol=4))
md.pattern2(dat)


#sample sizes in the specific pattern
#^
#^               numbers of missing data in each pattern
#|                ^
#|                |
#    X2 X3 X4 X1 NA
#6    1  1  1  1  0
#2    0  0  1  1  2
#1    0  0  0  1  3
#1    0  1  1  1  1
#Sum  4  3  1  0  8 --> numbers of missing data in each variable

}

\keyword{array}
