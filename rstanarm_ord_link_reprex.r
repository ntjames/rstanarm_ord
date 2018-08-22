# reprex
rm(list=ls())
library(dplyr)
library(rstanarm)

generate.data <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}

dat10levs <- generate.data(seed=45232,n=100) %>% mutate(y=cut_number(y,10))


fit_logit <- stan_polr(ordered(y) ~ z1+z2, data = dat10levs, 
                    prior=NULL, prior_counts = dirichlet(1),
                    method="logistic", adapt_delta = 0.999,
                    init_r=0.1)

fit_logit
foo<-summary(fit_logit)
posterior_vs_prior(fit_logit, regex_pars = "\\|")
#posterior_vs_prior(fit_logit, regex_pars = "\\|",color_by = "vs", group_by_parameter = TRUE)
posterior_predict(fit_logit)

fit_probit <- stan_polr(ordered(y) ~ z1+z2, data = dat10levs, 
                       prior=NULL, prior_counts =dirichlet(1),
                       method="probit", init_r=1, adapt_delta = 0.999)

plot(fit_probit)

fit_loglog <- stan_polr(ordered(y) ~ z1+z2, data = dat10levs, 
                       prior=NULL, prior_counts =dirichlet(1),
                       method="loglog", init_r=1, adapt_delta = 0.999)

plot(fit_loglog)

fit_cloglog <- stan_polr(ordered(y) ~ z1+z2, data = dat10levs, 
                       prior=NULL, prior_counts =dirichlet(1),
                       method="cloglog", init_r=1, adapt_delta = 0.999)

plot(fit_cloglog)
