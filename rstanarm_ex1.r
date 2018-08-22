# Try using rstanarm
# https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library(rstanarm)
library(tidyverse)

# from Liu et al. sim code
generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=log.y, z1=z1, z2=z2)
  return(data)
}


#set cores
options(mc.cores = parallel::detectCores())

lev15<-generate.data.1(seed=45232,n=1000) %>% mutate(y=cut_number(y,15))

lev5<-generate.data.1(seed=45232,n=1000) %>% mutate(y=cut_number(y,5))

# takes a few minutes to run
sfit_1 <- stan_polr(ordered(y) ~ z1+z2, data = lev5, 
                   prior=NULL, prior_counts =dirichlet(1),
                   adapt_delta = 0.99)


# check different methods to summarize fit object
prior_summary(sfit_1)

# or just sfit_1; short summary of model info, coefs, cutpoints/intercepts
print(sfit_1) 

# plot point ests and intervals for coefs and cutpoints  
plot(sfit_1)

# point est of coefs
coef(sfit_1)


fitted(sfit_1)

residuals(sfit_1)

nobs(sfit_1)

se(sfit_1)

vcov(sfit_1)

# model info, estimates w/ mean, sd, quantiles (0.025,0.25,0.5,0.75,0.975),
# diagnostics (mcse, Rhat, n_eff)
summary(sfit_1)


# get leave-one-out (loo) cross-validation 
loo(sfit_1)

#diagnostic plots for loo
plot(loo(sfit_1))

#posterior predictive check
pp_check(sfit_1)


# get posterior interval estimates
posterior_interval(sfit_1)

# draw from posterior predictive distribution
# dimension # draws x # observations 
post_draws<-posterior_predict(sfit_1)
dim(post_draws)

#posterior draws of the linear predictor

nd <- data.frame(z1=1,z2=0)
post_lp <- posterior_linpred(sfit_1, newdata=nd)

# get probs from linear predictor
plogis(post_lp) %>% qplot()
