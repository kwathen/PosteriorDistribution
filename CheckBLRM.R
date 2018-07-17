
library(mcmc)

# BLRM model with a bivariate normal prior
blrm <- function(beta, prior_mean, prior_cov, x, y) {
    eta <- beta[1] + x * exp(beta[2])
    
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    
    log_likelihood <- sum(logp[y == 1]) + sum(logq[y == 0])
    log_prior = - 0.5 * t(beta - prior_mean) %*% solve(prior_cov) %*%  (beta - prior_mean) 
    return(log_likelihood + log_prior)
}

LogPrior<- function(beta, prior_mean, prior_cov) {
    
    log_prior = - 0.5 * t(beta - prior_mean) %*% solve(prior_cov) %*%  (beta - prior_mean) 
    return(log_prior)
}


# BLRM model with a bivariate normal prior
LogLike <- function(beta, x, y) {
    eta <- beta[1] + x * exp(beta[2])
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    log_likelihood <- sum(logp[y == 1]) + sum(logq[y == 0])
    return(log_likelihood )
}

# Prior distribution
prior_mean = c(log(0.05 / 0.95), 0)

sd = c(2, 1)
corr = 0
prior_cov = matrix(c(sd[1]**2, sd[1]*sd[2]*corr, 
                     sd[1]*sd[2]*corr, sd[2]**2), ncol = 2, byrow = TRUE)

prior_mean
prior_cov
beta = c( -1.2, -.5)
log_prior = -  t(beta - prior_mean) %*% solve(prior_cov) %*%  (beta - prior_mean) 
log_prior
# Three patients are given the same dose of 1 mg
x = c(0,0,0,1, 1, 1,2,2,2)
# There are no dose-limiting toxicities
y = c(0, 0, 0, 1,0,0,0,1,1)

LogPrior( beta, prior_mean, prior_cov  )
LogLike( beta, x, y)

set.seed(7898)   
beta.init = rnorm( 2,0,1)
# The Scale parameter is chosen to get a reasonable acceptance rate 
bayes_out = metrop(blrm, beta.init, 10000, prior_mean = prior_mean, prior_cov = prior_cov, x = x, y = y, scale = 2)

# Diagnostics
plot(ts(bayes_out$batch))

# Acceptance rate
bayes_out$accept

# Sample from the posterior distribution
#bayes_out$batch
dim( bayes_out$batch)
mean( bayes_out$batch[,1])
mean( bayes_out$batch[,2])

median( bayes_out$batch[,2])
mean( ifelse(bayes_out$batch[,2] > 0, 1,0 ))

bayes_out = metrop(blrm, beta.init, 10000, prior_mean = prior_mean, prior_cov = prior_cov, x = x, y = y, scale = 2)

mean( ifelse(bayes_out$batch[,2] > 0, 1,0 ))
