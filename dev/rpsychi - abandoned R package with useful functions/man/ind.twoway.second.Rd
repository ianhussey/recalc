\encoding{UTF-8}
\name{ind.twoway.second}
\alias{ind.twoway.second}

\title{
A two-way design with independent samples using published work
}
\description{
    \code{ind.twoway.second} conducts a two-way design with independent samples, namely two-way randomized-group analysis of variance, using published work. 
}
\usage{
    ind.twoway.second(m, sd, n, 
            unbiased = TRUE, sig.level = 0.05, digits = 3)
}

\arguments{
  \item{m}{a matrix contains the means}
  \item{sd}{a matrix contains the sample/unbiased standard deviations}
  \item{n}{a matrix contains the sample size}
  \item{unbiased}{\code{sd} contains unbiased standard deviations (\code{unbiased} = TRUE, default) or sample standard deviations (\code{unbiased} = FALSE)}
  \item{sig.level}{a numeric contains the significance level (default 0.05)}
  \item{digits}{the specified number of decimal places (default 3)}
}


\details{
This function conducts a two-way design with independent samples, namely two-way randomized-group analysis of variance, using published work.
Statistical power is calculated using the following specifications: 

(a) small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), and large (\eqn{\eta^2 = 0.14}) population effect sizes, 
 according to the interpretive guideline for effect sizes by Cohen (1992)

(b) sample size specified by \code{n}

(c) significance level specified by \code{sig.level}

}

\value{
  The returned object of \code{ind.oneway.second} contains the following components:
  \item{anova.table}{returns a ANOVA table containing sums of squares, degrees of freedom, mean squares, \eqn{F} values}
  \item{omnibus.es}{returns a omnibus effect sizes which is a partial \eqn{\eta^2}, and its' confidence interval for each main and interaction effect}
  \item{power}{returns statistical power for detecting small (\eqn{\eta^2 = 0.01}), medium (\eqn{\eta^2 = 0.06}), 
  and large (\eqn{\eta^2 = 0.14}) population effect sizes}
}

\references{
Cohen B (2000) Calculating a factorial ANOVA from means and standard deviations. Understanding Statistics, 1, 191-203.

Cohen J (1992) A power primer. Psychological Bulletin, 112, 155-159.

Kline RB (2004) Beyond significance testing: Reforming data analysis methods in behavioral research. Washington: American Psychological Association.

Tabachnick BG, Fidell LS (2007) Experimental designs using ANOVA. Belmont, CA: Thomson. 
}


\author{Yasuyuki Okumura\cr 
        Department of Social Psychiatry, \cr 
        National Institute of Mental Health, \cr 
        National Center of Neurology and Psychiatry \cr 
        \email{yokumura@blue.zero.jp} 
} 


\seealso{\code{\link{ind.twoway}}}

\examples{
##Cohen (2000) Table 1
m.mat  <- matrix(c(37.13, 39.31, 39.22, 32.71), ncol=2) #2 * 2
sd.mat <- matrix(c(13.82, 9.42, 9.43, 9.62), ncol=2)
n.mat <- matrix(c(9, 13, 8, 14), ncol=2)

ind.twoway.second(m = m.mat, sd = sd.mat, n = n.mat)


##Tabachnick and Fidell (2007) 
#5.7 Complete example of two-way randomized-groups ANOVA (p.221-236)
m.mat <- matrix(c(837.9, 573.6, 354.9, 699.0, 112.0, 
      852.2, 781.6, 683.3, 1193.9, 130.0), ncol=2)    #5 * 2
sd.mat <- matrix(c(189.87449, 61.31195, 147.93351, 128.51891, 43.36922, 
      227.17042, 104.81221, 116.25934, 198.36692, 37.64158), ncol=2) #5 * 2
n.mat <- matrix(rep(10, 10), ncol=2)

ind.twoway.second(m = m.mat, sd = sd.mat, n = n.mat)


##Kline (2004) Table 7.5
dat <- data.frame(
           y = c(2,3,4,1,3,1,3,4,5,5,6,6,6,7),
           A = factor(c(rep("A1",5), rep("A2", 9))),
           B = factor(c(rep("B1",3), rep("B2",2), rep("B1",2), rep("B2",7)))
           )
ind.twoway.second(m = tapply(dat$y, list(dat$A,dat$B), mean), 
                  sd = tapply(dat$y, list(dat$A,dat$B), sd), 
                  n = tapply(dat$y, list(dat$A,dat$B), length)
                    )
}


\keyword{design}
\keyword{htest}

