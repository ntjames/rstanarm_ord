rm(list=ls())

library(rstanarm)
library(rms)
library(dplyr)
library(stringr)

dir<-getwd()
source(file.path(dir,"rstanarm_ord_functions.r"))

##  check functions  ##

# modified from Liu et al. sim code
generate.data.2 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  y <- rnorm(n, alpha+beta[1]*z1, sigma)
  data <- data.frame(y=y, z1=z1)
  return(data)
}


generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data)
}


dat25levs <- generate.data.2(seed=2458,n=25,beta=c(2,0))
dat50levs <- generate.data.2(seed=2458,n=50,beta=c(2,0))
dat100levs <- generate.data.2(seed=2458,n=100,beta=c(2,0))
dat200levs <- generate.data.2(seed=2458,n=200,beta=c(2,0))

dat50levs2 <- generate.data.1(seed=2458,n=50,beta=c(2,-0.5))


#dat100levs$y <- factor(dat100levs$y)

#with(subset(dat100levs,z1==1),hist( as.numeric(levels(dat100levs$y)) ))
#with(subset(dat100levs,z1==0),hist( as.numeric(levels(dat100levs$y)) ))
#hist(dat100levs$y)

# call this once to distribute chains across cpu cores:
options(mc.cores=parallel::detectCores())

#possible methods (links): "logistic", "probit", "loglog", "cloglog", "cauchit"
rstanarm_fit_probit <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
                                 prior=NULL, prior_counts = dirichlet(1), method="probit",
                                 adapt_delta = 0.99, chains=1)

rstanarm_fit_probit2 <- stan_polr(ordered(y) ~ z1, data = dat50levs2, 
                                  prior=NULL, prior_counts = dirichlet(1), method="probit",
                                  adapt_delta = 0.99, chains=1)

rstanarm_fit_probit3 <- stan_polr(ordered(y) ~ z1 + z2, data = dat50levs2, 
                                 prior=NULL, prior_counts = dirichlet(1), method="probit",
                                 adapt_delta = 0.99, chains=1)

#posterior_vs_prior(rstanarm_fit_probit)
#rstanarm_fit_probit2$stanfit@par_dims$cutpoints
#rstanarm_fit_probit2$stanfit@stan_args
#rstanarm_fit_probit2$stanfit@stanmodel

# prior is prior location of R^2 
# (prior for prop. of variance explained by model in underlying latent variable )
# prior_counts is prior counts of outcome when the predictors are at sample mean
#rstanarm_fit_logit <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
#                      prior=R2(0.5,"mean"), prior_counts = dirichlet(1), method="logistic",
#                      adapt_delta = 0.99)

# other links
if (0){
  rstanarm_fit_logit <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
                                  prior=NULL, prior_counts = dirichlet(1), method="logistic",
                                  adapt_delta = 0.99, chains=1)
  
  rstanarm_fit_loglog <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
                                  prior=NULL, prior_counts = dirichlet(1), method="loglog",
                                  adapt_delta = 0.99, chains=1)
  
  rstanarm_fit_cloglog <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
                                   prior=NULL, prior_counts = dirichlet(1), method="cloglog",
                                   adapt_delta = 0.99, chains=1)
  
  rstanarm_fit_cauchit <- stan_polr(ordered(y) ~ z1, data = dat50levs, 
                                    prior=NULL, prior_counts = dirichlet(1), method="cauchit",
                                    adapt_delta = 0.99, chains=1)
}

#true dist. for z1=0
norm_cdf_z10 <- curve(pnorm(x,0,1),-4.5,4.5)
n_cdf_0 <- data.frame(yval=norm_cdf_z10$x, cdf=norm_cdf_z10$y)

#true dist. for z1=1
norm_cdf_z11 <- curve(pnorm(x,2,1),-4.5,4.5)
n_cdf_1 <- data.frame(yval=norm_cdf_z11$x, cdf=norm_cdf_z11$y)


## conditional CDF ##

cond_cdf_samps <- getCDF(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),summ=FALSE)
cond_cdf_summ <- getCDF(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)))

cond_cdf_summ

#png("cdf1.png")
cond_cdf_summ %>% ggplot(aes(group=z1)) +
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5), fill="grey30", alpha=0.4) +
  geom_step(aes(x=yval,y=med_cdf)) +
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) +
  xlab("") + ylab("Conditional CDF")
#dev.off()

# empirical cdfs
Fnz1_0 <- dat50levs %>% filter(z1==0) %>% pull(y) %>% ecdf()
Fnz1_1 <- dat50levs %>% filter(z1==1) %>% pull(y) %>% ecdf()

# z1=0 
plot(Fnz1_0);curve(pnorm(x),-3,3, add=TRUE, col="blue")

cond_cdf_summ %>% filter(z1==0) %>% ggplot()+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=med_cdf))+
  geom_line(data=n_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  xlab("") + ylab("Conditional CDF")

# z1=1 
plot(Fnz1_1);curve(pnorm(x,2,1),-3,5, add=TRUE, col="blue")

cond_cdf_summ %>% filter(z1==1) %>% ggplot()+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=med_cdf))+
  geom_line(data=n_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  xlab("") + ylab("Conditional CDF")


# verify link is correct

cond_cdf_samps %>% filter(ndrow==1) %>% 
  ggplot(aes(x=yval,y=Freq,group=Var1)) + geom_line(alpha=0.01) +
  geom_abline(slope=1,intercept = 0)


# log

cond_cdf_samps2 <- getCDF(rstanarm_fit_probit2,newdata=data.frame(z1=c(0)),summ=FALSE)

cond_cdf_samps2 %>%
  ggplot(aes(x=log(yval),y=Freq,group=Var1)) + geom_line(alpha=0.01) +
  geom_abline(slope=1,intercept = 0)


cond_cdf_summ2 <- getCDF(rstanarm_fit_probit2,newdata=data.frame(z1=c(0,1)))


# two covariates
cond_cdf_summ3 <- getCDF(rstanarm_fit_probit3,newdata=data.frame(z1=c(0,1),z2=c(0,0)))




## -- conditional Mean -- ##

cond_mean_samps <- getMean(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)), summ=FALSE )
cond_mean_summ <- getMean(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)) )
cond_mean_summ

cond_mean_samps %>% ggplot(aes(x=mn))+geom_histogram()+facet_grid(~ndrow)

getMean(rstanarm_fit_probit2,newdata=data.frame(z1=c(0,1),z2=c(0,0)))


## -- conditional Quantiles -- ##

# head(posterior_samples(brm_fit_logit))

cond_med_samps<-getQuantile(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),q=0.5, summ=FALSE)
cond_med_summ<-getQuantile(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),q=0.5)
cond_med_summ

cond_med_summ2<-getQuantile(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),q=0.75)
cond_med_summ2

cond_med_samps3<-getQuantile(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),q=0.01, summ=FALSE)

cond_med_summ3<-getQuantile(rstanarm_fit_probit,newdata=data.frame(z1=c(0,1)),q=0.01)
cond_med_summ3
