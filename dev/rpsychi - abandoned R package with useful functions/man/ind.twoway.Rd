\name{ind.twoway}
\alias{ind.twoway}

\title{
A two-way design with independent samples using individual data
}

\description{
    \code{ind.twoway} conducts a two-way design with independent samples, namely two-way randomized-group analysis of variance, using individual data. 
}
\usage{
ind.twoway(formula, data, sig.level=.05, digits=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{two-sided formula; the left-hand-side of which gives one dependent variable containing a numeric variable, 
  and the right-hand-side of two independent variables containing a factor with two or more levels}
  \item{data}{a data frame contains the variables in the \code{fomrmula}}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}


\details{
This function conducts a two-way design with independent samples, namely two-way randomized-group analysis of variance, using individual data.
The main effect of the first independent variable will be shown in "(row)" section, 
and the main effect of the second one will be shown in "(col)" section.
Always the interaction effect of the two independent variables will be calculated.
Statistical power is calculated using the following specifications: 

(a) small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), and large (\eqn{\eta^2 = 0.14}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{formula} and \code{data}

(c) significance level specified by \code{sig.level}
}


\value{
  The returned object of \code{ind.oneway} contains the following components:
  \item{anova.table}{returns a ANOVA table containing sums of squares, degrees of freedom, mean squares, \eqn{F} values}
  \item{omnibus.es}{returns a omnibus effect sizes which is a partial \eqn{\eta^2}, and its' confidence interval for each main and interaction effect}
  \item{power}{returns statistical power for detecting small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), 
  and large (\eqn{\eta^2 = 0.14}) population effect sizes}
}

\references{
Cohen B (2000) Calculating a factorial ANOVA from means and standard deviations. Understanding Statistics, 1, 191-203.

Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.

}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} 
} 


\seealso{\code{\link{ind.twoway.second}}}

\examples{
##Kline (2004) Table 7.5
dat <- data.frame(
           y = c(2,3,4,1,3,1,3,4,5,5,6,6,6,7),
           A = factor(c(rep("A1",5), rep("A2", 9))),
           B = factor(c(rep("B1",3), rep("B2",2), rep("B1",2), rep("B2",7)))
           )

ind.twoway(y~A*B, data=dat)
}

\keyword{design}
\keyword{htest}


