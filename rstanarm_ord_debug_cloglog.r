rm(list=ls())

library(rstanarm)
#library(rms)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)

dir<-getwd()
source(file.path(dir,"rstanarm_ord_functions.r"))

# functions to generate data 

# modified from Liu et al. sim code
generate.data.norm <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  y <- rnorm(n, alpha+beta[1]*z1, sigma)
  data <- data.frame(y=y, z1=z1)
  return(data)
}

generate.data.logis <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  y <- rlogis(n, alpha+beta[1]*z1, sigma)
  data <- data.frame(y=y, z1=z1)
  return(data)
}

rGumbelMin <- function(n, mu=0, sigma=1){
  u<-runif(n,min=0,max=1)
  x<-mu+sigma*log(-log(1-u))
  return(x)
}

pGumbelMin <- function(x, mu=0, sigma=1){
  pf <- 1-exp(-exp((x-mu)/sigma))
  return(pf)
}

generate.data.ev1 <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  error <-rGumbelMin(n, mu=0, sigma=1)-digamma(1)
  y <- alpha+beta[1]*z1 +error
  
  data <- data.frame(y, z1=z1)
  return(data) 
}


#hist(pGumbelMin(rGumbelMin(1e6)))


rGumbelMax <- function(n, mu=0, sigma=1){
  u<-runif(n,min=0,max=1)
  x<-mu+sigma*(-log(-log(u)))
  return(x)
}

pGumbelMax <- function(x, mu=0, sigma=1){
  pf <- exp(-exp(-(x-mu)/sigma))
  return(pf)
}

# check
#hist(pGumbelMax(rGumbelMax(1e6)))

generate.data.ev2 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  error <- rGumbelMax(n, 0, 1)+digamma(1) # not sure why digamma() added here
  y <- alpha+beta[1]*z1 +error
  data <- data.frame(y=y, z1=z1)
  return(data)
}

# true beta for z1=1
bta<-3

### generate normal data
dat25levs_norm <- generate.data.norm(seed=2417,n=25,beta=c(bta,0))

## generate logistic data
dat40levs_logis <- generate.data.logis(seed=2417,n=40,beta=c(bta,0))

# generate extreme value type 1 (Gumbel min) ?
dat25levs_ev1 <- generate.data.ev1(seed=2458,n=25,beta=c(bta,0))
dat50levs_ev1 <- generate.data.ev1(seed=2458,n=50,beta=c(bta,0))  

# generate extreme value type 2
dat40levs_ev2 <- generate.data.ev2(seed=2417,n=40,beta=c(bta,0))


# call this once to distribute chains across cpu cores:
options(mc.cores=parallel::detectCores())

#possible methods (links): "logistic", "probit", "loglog", "cloglog", "cauchit"

# Note about parameterization from documentation: 
# The Dirichlet distribution is used in stan_polr for an implicit prior on the cutpoints 
# in an ordinal regression model. More specifically, the Dirichlet prior pertains to the 
# prior probability of observing each category of the ordinal outcome when the predictors 
# are at their sample means. Given these prior probabilities, it is straightforward to add 
# them to form cumulative probabilities and then use an inverse CDF transformation of the 
# cumulative probabilities to define the cutpoints.

# If a scalar is passed to the concentration argument of the dirichlet function, then it is 
# replicated to the appropriate length and the Dirichlet distribution is symmetric. If 
# concentration is a vector and all elements are 1, then the Dirichlet distribution is 
# jointly uniform. If all concentration parameters are equal but greater than 1 then the 
# prior mode is that the categories are equiprobable, and the larger the value of the 
# identical concentration parameters, the more sharply peaked the distribution is at the 
# mode. The elements in concentration can also be given different values to represent that 
# not all outcome categories are a priori equiprobable.


# start w/ loglog
if (0){
fit_tmp1<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
               prior=NULL, prior_counts = dirichlet(1), method="loglog",
               adapt_delta = 0.99, chains=1, iter=2000)

fit_tmp2<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=NULL, prior_counts = dirichlet(1.1), method="loglog",
                    adapt_delta = 0.99, chains=1, iter=2000)

fit_tmp3<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=R2(0.5,"mean"), prior_counts = dirichlet(1), method="loglog",
                    adapt_delta = 0.99, chains=1, iter=2000)


fit_tmp4<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=R2(0.5,"mean"), prior_counts = dirichlet(1), method="loglog",
                    adapt_delta = 0.99, chains=1, iter=2000, verbose=TRUE)

cbind(summary(fit_tmp1)[1:25,1],summary(fit_tmp2)[1:25,1],summary(fit_tmp3)[1:25,1])
}

# debug cloglog
# try just drawing from prior
fit_tmp5<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=NULL, prior_counts = dirichlet(1), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, 
                    prior_PD=TRUE)

attr(summary(fit_tmp5),"family")

# pass arg to sampling to get additional diagnostics (verbose=TRUE)
if(0){
fit_tmp6<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=NULL, prior_counts = dirichlet(1), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, verbose=TRUE)
}

# set refresh=1 to see if any samples are being drawn
if(0){
fit_tmp7<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=NULL, prior_counts = dirichlet(1), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, refresh=5)

#! runs very slowly and eventually hangs/freezes R session
}

# try with informative counts (dirichlet)
if(0){
fit_tmp8<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=NULL, prior_counts = dirichlet(1.5), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, refresh=5)
# hangs
}


# try reducing adapt_delta
if(0){
fit_tmp9<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=R2(0.5,"mean"), prior_counts = dirichlet(1), method="cloglog",
                    adapt_delta = NULL, chains=1, iter=2000, refresh=5)
# still very slow - not sure how long to finish
}

# try with informative R2
fit_tmp10<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    prior=R2(0.5,"mean"), prior_counts = dirichlet(1), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, refresh=1)


# try with both informative priors
fit_tmp11<-stan_polr(ordered(y) ~ z1, data = dat25levs_norm, 
                    priorR2(0.5,"mean"), prior_counts = dirichlet(1.5), method="cloglog",
                    adapt_delta = 0.99, chains=1, iter=2000, refresh=1)

# extract code and check
str(fit_tmp5$stanfit)

# control args?
#test_grad
#adapt_term_buffer?
#adapt_window?

str(fit_tmp5$stanfit@stanmodel)

write(fit_tmp4$stanfit@stanmodel@model_code[1],file=file.path(dir,"loglog_mod_code.txt"))
      
write(fit_tmp5$stanfit@stanmodel@model_code[1],
      file=file.path(dir,"cloglog_prior_PD_mod_code.txt"))
