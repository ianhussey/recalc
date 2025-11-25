
# source of equations
# https://matthewbjane.quarto.pub/pre-post-correlations/#scenario-3-is-the-t-statistic-from-a-paired-t-test-available


# simulate data
set.seed(1)

# parameters
n      <- 200
rho     <- 0.60  # moderate pre–post correlation

# target means/SDs
mu      <- c(20, 16)
sigma   <- c(8, 8)

# covariance matrix
Sigma <- matrix(c(sigma[1]^2, rho*sigma[1]*sigma[2],
                  rho*sigma[1]*sigma[2], sigma[2]^2), nrow = 2)

# simulate underlying continuous scores
library(MASS)
cont <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)

# impose Likert bounds and integer responses
pre  <- cont[,1]
post <- cont[,2]

dat <- data.frame(id = 1:n, pre, post)

# check empirical moments
summary(dat)
mean(dat$pre); mean(dat$post)
# [1] 20.10891
# [1] 16.39969

sd(dat$pre); sd(dat$post)

# [1] 7.634948
# [1] 7.500166

nrow(dat)
# [1] 200

# empirical correlation
cor(dat)[2,3]
#[1] 0.5434789

# paired Student t-test
t_res <- t.test(dat$pre, dat$post, paired = TRUE, var.equal = TRUE)
t_res

# t = 7.2533, df = 199, p-value = 1.383e-11
# mean difference 3.70922

library(effectsize)
repeated_measures_d(dat$pre, dat$post, method = "z")
# 0.51



# option 1: using cohens d_z
mean_pre = 20.11
mean_post = 16.40
sd_pre = 7.63
sd_post = 7.50
n = 200
dz = 0.51

r <- (sd_pre^2 + sd_post^2 - ((mean_post - mean_pre)/dz)^2) /
  (2 * sd_pre * sd_post)

r
# .54




# option 2: Is the t-statistic from a paired t-test available?
prepost_cor_from_dependent_ttest <- function(paired_t, sd_pre, sd_post, n, mean_pre, mean_post){

  # old - wrong?
  # r <- 
  #   (t^2 * (SD1^2 + SD2^2)) - (n * m_diff) / 
  #   (2 * t^2 * SD1 * SD2)
  
  r <- 
    (paired_t^2*(sd_pre^2 + sd_post^2)-n*(mean_post-mean_pre)^2) /
    (2*paired_t^2*sd_pre*sd_post)
  
  return(r)
}

prepost_cor_from_dependent_ttest(paired_t = 7.25, 
                                 mean_pre = 20.11, 
                                 mean_post = 16.40, 
                                 sd_pre = 7.63, 
                                 sd_post = 7.50, 
                                 n = 200)
# r = 0.5425477




# option 3: Is the p-value from a paired t-test reported?
paired_p <- 8.856e-12 # requires highly precise p value
mean_pre = 20.11
mean_post = 16.40
sd_pre = 7.63
sd_post = 7.50
n = 200

# get paired t from p value
paired_t <- qt(paired_p/2, n-1, lower.tail = FALSE)

r <- (paired_t^2*(sd_pre^2 + sd_post^2)-n*(mean_post-mean_pre)^2) /
  (2*paired_t^2*sd_pre*sd_post)
r
