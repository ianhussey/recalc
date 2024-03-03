\encoding{UTF-8}
\name{groupSummary}
\alias{groupSummary}

\title{Compute summary statistics by group}

\description{
    \code{groupSummary} summarize continuous variables by providing means and sample standard deviations 
    and categorical variables by giving numbers and percentages for each group. 
}

\usage{
groupSummary(data, group = NULL, digits = 3)
}


\arguments{
  \item{data}{a data frame contains the variables}
  \item{group}{a character string specify the grouping variable (default NULL)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function summarize continuous variables by providing means and sample standard deviations 
    and categorical variables by giving numbers and percentages for each group. 
    Missing values (NA) will be stripped before the computation proceeds.
    If you do not specify \code{group}, this function will summarize the variables for total group only.
}

\value{
  Return a data frame containing (1) means and/or numbers, (2) sample standard deviations and/or percentages, (3) numbers of available data for total and each group.
}

\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 

\seealso{\code{\link{ssd}}}

\examples{
data(infert)
infert$case <- factor(infert$case, labels=c("control", "case"))
infert$induced <- factor(infert$induced, labels=c("0","1","2 or more"))
infert$spontaneous <- factor(infert$spontaneous, labels=c("0","1","2 or more"))

#continuous and categorical variables
groupSummary(infert, group="case")

#continuous variables only
groupSummary(infert[, c(2,3,7,8, 5)], group="case")

#categorical variables only
groupSummary(infert[, c(1,4, 6, 5)],  group="case")    

#total sample
groupSummary(infert[, c(1,4, 6, 5)])

}

\keyword{univar}
