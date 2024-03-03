\encoding{UTF-8}
\name{multreg}
\alias{multreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A multiple regression analysis using individual data}

\description{
   \code{multreg} conducts a multiple regression analysis using individual data.
}

\usage{
multreg(formula, data, sig.level = 0.05, digits = 3)
}

\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of several independent variables containing a numeric variable}
  \item{data}{a data frame contains the variables in the \code{fomrmula}}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}
\details{
This function conducts a multiple regression analysis using individual data.
The dependent variable and independent variables should be a numeric vector. 
In this function, you cannot specify any interaction nor any curvilinear effect.
Statistical power is calculated using the following specifications: 

(a) small (\eqn{R^{2} = 0.02}), medium (\eqn{R^{2} = 0.13}), and large (\eqn{R^{2} = 0.26}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{data}

(c) significance level specified by \code{sig.level}
 
(d) numbers of independent variable specified by \code{formula}
}


\value{
  \item{samp.stat}{returns the means and unbiased standard deviations}
  \item{corr.partial.corr}{returns a product-moment correlation matrix (lower triangle) and a partial correlation matrix given all remaining variables (upper triangle)}
  \item{corr.confidence}{returns lower and upper confidence limits (lower and upper triangles, respectively)}
  \item{omnibus.es}{returns a coefficient of determination and its' confidence interval}
  \item{raw.estimates}{returns partial regression coefficients, their confidence intervals, and standard errors}
  \item{standardized.estimates}{returns standardized partial regression coefficients, their confidence intervals, and standard errors}
  \item{power}{returns statistical power for detecting small (\eqn{R^{2} = 0.02}), medium (\eqn{R^{2} = 0.13}), and large (\eqn{R^{2} = 0.26}) population effect sizes}
}

\references{
Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Cohen J, Cohen P, Aiken LS (2003) Applied multiple regression/correlation analysis for the behavioral sciences (3rd ed). Mahwah, NJ: Erlbaum.

Smithson M (2001) Correct confidence intervals for various regression effect sizes and parameters: The importance of noncentral distributions in computing intervals, 61, 605-632.
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} } 


\seealso{\code{\link{multreg.second}}, \code{\link{samplesize.rsq}}}

\examples{
##Cohen (2003) Table 3.5.1
dat <- data.frame(
salary = c(51876, 54511, 53425, 61863, 52926, 47034, 66432, 61100, 41934, 
  47454, 49832, 47047, 39115, 59677, 61458, 54528, 60327, 56600, 
  52542, 50455, 51647, 62895, 53740, 75822, 56596, 55682, 62091, 
  42162, 52646, 74199, 50729, 70011, 37939, 39652, 68987, 55579, 
  54671, 57704, 44045, 51122, 47082, 60009, 58632, 38340, 71219, 
  53712, 54782, 83503, 47212, 52840, 53650, 50931, 66784, 49751, 
  74343, 57710, 52676, 41195, 45662, 47606, 44301, 58582),
pubs  = c(18, 3, 2, 17, 11, 6, 38, 48, 9, 22, 30, 21, 
  10, 27, 37, 8, 13, 6, 12, 29, 29, 7, 6, 69, 11, 9, 
  20, 41, 3, 27, 14, 23, 1, 7, 19, 11, 31, 9, 12, 32, 
  26, 12, 9, 6, 39, 16, 12, 50, 18, 16, 5, 20, 50, 
  6, 19, 11, 13, 3, 8, 11, 25, 4),
cits = c(50, 26, 50, 34, 41, 37, 48, 56, 19, 29, 
    28, 31, 25, 40, 61, 32, 36, 69, 47, 29, 35, 
    35, 18, 90, 60, 30, 27, 35, 14, 56, 50, 25, 
    35, 1, 69, 69, 27, 50, 32, 33, 45, 54, 47, 29, 
    69, 47, 43, 55, 33, 28, 42, 24, 31, 27, 
    83, 49, 14, 36, 34, 70, 27, 28)   
)
multreg(salary~ pubs + cits, data=dat)
}


\keyword{design}
\keyword{htest}

