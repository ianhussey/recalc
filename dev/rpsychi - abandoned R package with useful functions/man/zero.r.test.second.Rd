\encoding{UTF-8}
\name{zero.r.test.second}
\alias{zero.r.test.second}

\title{A significance testing of a product moment correlation using published work}

\description{
    \code{zero.r.test.second} conducts a significance testing of a product moment correlation using published work.
}


\usage{
zero.r.test.second(r, n, sig.level = 0.05, digits = 3)
}

\arguments{
  \item{r}{a numeric contains the product moment correlation}
  \item{n}{a numeric contains the sample size}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}

\details{
This function conducts a significance testing of a product moment correlation using published work.
Statistical power is calculated using the following specifications: 

(a) small (\eqn{r = 0.10}), medium (\eqn{r = 0.30}), and large (\eqn{r = 0.50}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{n}

(c) significance level specified by \code{sig.level}

}


\value{
  The returned object of \code{zero.r.test} contains the following components:
  \item{correlation}{returns a product moment correlation, its' approximate confidence interval for population correlation, and standard error}
  \item{power}{returns statistical power for detecting 
small (\eqn{r = 0.10}), medium (\eqn{r = 0.30}), and large (\eqn{r = 0.50}) population effect sizes}
}




\references{
Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.
}




\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{zero.r.test}}, \code{\link{samplesize.r}}} 


\examples{
zero.r.test.second(r = 0.571, n = 9)
}

\keyword{design}
\keyword{htest}
