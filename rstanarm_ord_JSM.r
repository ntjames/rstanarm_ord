rm(list=ls())
libs<-c("rstanarm","dplyr","stringr","tibble","tidyr","cowplot","ggalt")
invisible(lapply(libs,library,character.only=TRUE))

plotdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/JSM2018/speed session/fig")

dir<-getwd()
source(file.path(dir,"rstanarm_ord_functions.r"))


# functions to generate data 
if(1){
  
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

# ? also see generate.data.c w/ digamma??
#! old
#generate.data.ev1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
#  set.seed(seed)
#  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
#  y <- rGumbelMin(n, mu=alpha+beta[1]*z1, sigma=sigma)
#  data <- data.frame(y=y, z1=z1)
#  return(data)
#}


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

# for log outcome
if(0){
 
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
}


}

# true beta for z1=1
bta<-3

### generate normal data
dat25levs_norm <- generate.data.norm(seed=2417,n=25,beta=c(bta,0))
dat50levs_norm <- generate.data.norm(seed=2417,n=50,beta=c(bta,0))

#dat50levs_norm %>% group_by(z1) %>% dplyr::summarize(mn=mean(y)) 

# true dist. for z1=0
norm_cdf_z10 <- curve(pnorm(x,0,1),-6.5,6.5)
n_cdf_0 <- data.frame(yval=norm_cdf_z10$x, cdf=norm_cdf_z10$y)

# true dist. for z1=1
norm_cdf_z11 <- curve(pnorm(x,bta,1),-6.5,6.5)
n_cdf_1 <- data.frame(yval=norm_cdf_z11$x, cdf=norm_cdf_z11$y)

rm(norm_cdf_z10,norm_cdf_z11)

## generate logistic data

dat25levs_logis <- generate.data.logis(seed=2417,n=25,beta=c(bta,0))
dat50levs_logis <- generate.data.logis(seed=2417,n=50,beta=c(bta,0))

#dat50levs_logis %>% group_by(z1) %>% dplyr::summarize(mn=mean(y)) 

# check true mean
if(0){
  dat10000levs_logis <- generate.data.logis(seed=245,n=10000,beta=c(bta,0))
  dat10000levs_logis %>% group_by(z1) %>% dplyr::summarize(mn=mean(y)) 
  
  mean(rlogis(1e5,0,1))
  mean(rlogis(1e5,3,1))
}

#true dist. for z1=0
logis_cdf_z10 <- curve(plogis(x,0,1),-6.5,7)
l_cdf_0 <- data.frame(yval=logis_cdf_z10$x, cdf=logis_cdf_z10$y)

#true dist. for z1=1
logis_cdf_z11 <- curve(plogis(x,bta,1),-6.5,7)
l_cdf_1 <- data.frame(yval=logis_cdf_z11$x, cdf=logis_cdf_z11$y)

rm(logis_cdf_z10, logis_cdf_z11)


# generate extreme value type 1 (Gumbel min) ?
if (0){
dat25levs_ev1 <- generate.data.ev1(seed=2458,n=25,beta=c(bta,0))
dat50levs_ev1 <- generate.data.ev1(seed=2458,n=50,beta=c(bta,0))  

#true dist. for z1=0
ev1_cdf_z10 <- curve(pGumbelMin(x,0,1),-6.5,5.5)
ev1_cdf_0 <- data.frame(yval=ev1_cdf_z10$x, cdf=ev1_cdf_z10$y)

#true dist. for z1=1
ev1_cdf_z11 <- curve(pGumbelMin(x,bta,1),-6.5,5.5)
ev1_cdf_1 <- data.frame(yval=ev1_cdf_z11$x, cdf=ev1_cdf_z11$y)

rm(ev1_cdf_z10, ev1_cdf_z11)
}

# generate extreme value type 2
dat25levs_ev2 <- generate.data.ev2(seed=2417,n=25,beta=c(bta,0))
dat50levs_ev2 <- generate.data.ev2(seed=2417,n=50,beta=c(bta,0))

#dat50levs_ev2 %>% group_by(z1) %>% dplyr::summarize(mn=mean(y)) 
#dat50levs_ev2 %>% ggplot(aes(x=y))+ geom_histogram(bins=15) + facet_grid(~z1)

#!! check true mean and generate.data.ev2
if(0){
dat1e5levs_ev2 <- generate.data.ev2(seed=245,n=1e5,beta=c(bta,0))
dat1e5levs_ev2 %>% group_by(z1) %>% dplyr::summarize(mn=mean(y)) 

dat1e5levs_ev2%>% 
  ggplot(aes(x=y))+ geom_density() + facet_grid(~z1)

mean(rGumbelMax(1e5,0,1)+digamma(1))
mean(rGumbelMax(1e5,3,1)+digamma(1))

Fnz1_0 <- dat1e5levs_ev2 %>% filter(z1==0) %>% pull(y) %>% ecdf()
Fnz1_1 <- dat1e5levs_ev2 %>% filter(z1==1) %>% pull(y) %>% ecdf()

plot(Fnz1_0);curve(pGumbelMax(x-digamma(1),0,1),-6,5, add=TRUE, col="blue")
plot(Fnz1_1);curve(pGumbelMax(x-digamma(1),bta,1),-4,8, add=TRUE, col="blue")

}

#true dist. for z1=0
ev2_cdf_z10 <- curve(pGumbelMax(x-digamma(1),0,1),-6.5,6.5)
ev2_cdf_0 <- data.frame(yval=ev2_cdf_z10$x, cdf=ev2_cdf_z10$y)

#true dist. for z1=1
ev2_cdf_z11 <- curve(pGumbelMax(x-digamma(1),bta,1),-6.5,6.5)
ev2_cdf_1 <- data.frame(yval=ev2_cdf_z11$x, cdf=ev2_cdf_z11$y)

rm(ev2_cdf_z10, ev2_cdf_z11)

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


fitmod<-function(fitlink=NULL, moddata, chains=1, ...){
  nm<-paste0("fit_",fitlink,"_",moddata)
  print(nm)
  
  fit<-stan_polr(ordered(y) ~ z1, data = get(moddata), 
            prior=NULL, prior_counts = dirichlet(1), method=fitlink,
            adapt_delta = 0.99, chains=chains, ...)
  
  assign(nm,fit, pos=.GlobalEnv)
}




## - check invariance to (monotonic) transformation 
if (0){
dat1 <- generate.data.norm(seed=2417,n=50,beta=c(bta,0))
dat2 <- dat1 %>% mutate(y=exp(y))
dat3 <- dat1 %>% mutate(y=y^2)
dat4 <- dat1 %>% mutate(y=1/y) # not monotonic

tmp<-rnorm(50)
plot(tmp,exp(tmp))
plot(tmp,tmp^2)
plot(tmp,1/tmp)

fit1<-stan_polr(ordered(y) ~ z1, data = dat1, 
               prior=NULL, prior_counts = dirichlet(1), method="probit",
               adapt_delta = 0.99, chains=1, iter=2000)

fit2<-stan_polr(ordered(y) ~ z1, data = dat2, 
               prior=NULL, prior_counts = dirichlet(1), method="probit",
               adapt_delta = 0.99, chains=1, iter=2000)

fit3<-stan_polr(ordered(y) ~ z1, data = dat3, 
               prior=NULL, prior_counts = dirichlet(1), method="probit",
               adapt_delta = 0.99, chains=1, iter=2000)

fit4<-stan_polr(ordered(y) ~ z1, data = dat4, 
                prior=NULL, prior_counts = dirichlet(1), method="probit",
                adapt_delta = 0.99, chains=1, iter=2000)

coef(fit1)
coef(fit2)
coef(fit3)
coef(fit4) # transform not monotonic


# compare to regular linear model
fit4<-stan_lm(y ~ z1, data = dat1, prior=R2(0.5,"mean"),
                adapt_delta = 0.99, chains=1, iter=2000)

fit5<-stan_lm(y ~ z1, data = dat2, prior=R2(0.5,"mean"),
              adapt_delta = 0.99, chains=1, iter=2000)

fit6<-stan_lm(y ~ z1, data = dat3, prior=R2(0.5,"mean"),
              adapt_delta = 0.99, chains=1, iter=2000)

coef(fit4)
coef(fit5)
coef(fit6)
}

## - check with floor and ceiling effect
if (0){
dat5<- generate.data.norm(seed=2417,n=150,beta=c(1.5,0)) %>% arrange(y) %>% 
  mutate(ynew=ifelse(y<=-0.5,-0.5, y),
         ynew=ifelse(ynew>2.5,2.5,ynew))

fit7<-stan_polr(ordered(ynew) ~ z1, data = dat5, 
            prior=NULL, prior_counts = dirichlet(1), method="probit",
            adapt_delta = 0.99, chains=1, iter=2000)

fit7_comp<-stan_polr(ordered(y) ~ z1, data = dat5, 
                prior=NULL, prior_counts = dirichlet(1), method="probit",
                adapt_delta = 0.99, chains=1, iter=2000)

coef(fit7)
coef(fit7_comp)

cdf_150_norm_fc <- getCDF(fit7,newdata=data.frame(z1=c(0,1)))

cdf_150_norm_fc %>% ggplot(aes(group=z1)) +
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))

cdf_150_norm_fc_comp <- getCDF(fit7_comp,newdata=data.frame(z1=c(0,1)))

cdf_150_norm_fc_comp %>% ggplot(aes(group=z1)) +
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))

fit8<-stan_lm(ynew ~ z1, data = dat5, prior=R2(0.5,"mean"),
              adapt_delta = 0.99, chains=1, iter=2000)

coef(fit8)
}

## - annotated plot w/ Expectations and Quantiles

fitmod("logistic", "dat50levs_logis", chains=4, iter=5000)

prior_summary(fit_logistic_dat50levs_logis)

cdf_logistic_dat50levs_logis<-getCDF( fit_logistic_dat50levs_logis,
                                   newdata=data.frame(z1=c(0,1)) )

mean_logistic_dat50levs_samps <- getMean(fit_logistic_dat50levs_logis,
                                     newdata=data.frame(z1=c(0,1)), 
                                     summ=FALSE ) 

q75_logistic_dat50levs_samps <- getQuantile(fit_logistic_dat50levs_logis,
                                       newdata=data.frame(z1=c(0,1)), q=0.75,
                                       summ=FALSE ) 

q<-0.75
qy1_n<-qlogis(q,0,1)
qy2_n<-qlogis(q,bta,1)

post_mn<- mean_logistic_dat50levs_samps %>% group_by(ndrow) %>%
  mutate(x=ndrow-1,z1=x, ci95_lb=quantile(mn,0.025), ci95_ub=quantile(mn,0.975))

post_q75<-q75_logistic_dat50levs_samps %>% group_by(z1) %>%
  mutate(x=z1, ci95_lb=quantile(qtile,0.025), ci95_ub=quantile(qtile,0.975)) 


# make plots

cdf_plt<-cdf_logistic_dat50levs_logis %>% ggplot(aes(group=z1))+
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5),stat="stepribbon", direction="hv", fill="grey30", alpha=0.4)+
  geom_step(aes(x=yval,y=mn_cdf))+
  geom_line(data=l_cdf_0,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) + 
  geom_line(data=l_cdf_1,aes(x=yval,y=cdf),color="blue",inherit.aes = FALSE) +
  xlab("") + ylab("Conditional CDF") + 
  geom_segment(aes(x=-6.5 ,y=q, xend=qy1_n, yend=q),lty=3,color="blue")+
  geom_segment(aes(x=qy1_n,y=q, xend=qy1_n, yend=0),lty=3,color="blue")+
  geom_segment(aes(x=-6.5 ,y=q, xend=qy2_n, yend=q),lty=3,color="blue")+
  geom_segment(aes(x=qy2_n,y=q, xend=qy2_n, yend=0),lty=3,color="blue") + 
  geom_segment(data=post_mn, aes(x=ci95_lb ,y=-0.05, xend=ci95_ub, yend=-0.05),
               lty=6,color="forestgreen") +
  geom_segment(data=post_q75, aes(x=ci95_lb ,y=-0.13, xend=ci95_ub, yend=-0.13),
               lty=6,color="red") +
  annotate("text",x=-4,y=-0.05, label="Mean 95% PIs", size=4) +
  annotate("text",x=-4,y=-0.13, label="paste(Q[0.75],' ',95,'%','PIs')", parse=TRUE, size=4)+
  theme(text=element_text(size=20), axis.text = element_text(size=20))

#cdf_plt

post_mn_plt <- post_mn %>%
  ggplot(aes(x=mn))+geom_histogram(bins=40) +
  geom_vline(aes(xintercept=ci95_lb), lty=6, color="forestgreen")+
  geom_vline(aes(xintercept=ci95_ub), lty=6, color="forestgreen")+
  facet_grid(~x, labeller="label_both") + 
  xlab("")+
  ggtitle("Posterior Mean")

#post_mn_plt

post_q75_plt <-post_q75 %>%
  ggplot(aes(x=qtile))+geom_histogram(bins=40) +
  geom_vline(aes(xintercept=ci95_lb), lty=6, color="red")+
  geom_vline(aes(xintercept=ci95_ub), lty=6, color="red")+
  facet_grid(~x, labeller="label_both") +
  xlab("")+
  ggtitle(expression(paste("Posterior ", Q[0.75])) )

#post_q75_plt

pdf(file.path(plotdir,"cdf1.pdf"), height=4, width=4)
cdf_plt
dev.off()

pdf(file.path(plotdir,"mean_qtile1.pdf"), height=4, width=4)
plot_grid(post_mn_plt,post_q75_plt, nrow=2)
dev.off()

pdf(file.path(plotdir,"cdf2.pdf"), height=3.5, width=6)
cdf_plt
dev.off()



plot_grid(post_mn_plt,post_q75_plt, nrow=2)

post_mn_2<-post_mn %>% ungroup() %>% select(mn,x,ci95_lb,ci95_ub) %>% 
  mutate(stat="Mean", x=paste0("x=",x)) %>% rename(val=mn)

post_q75_2<-post_q75 %>% ungroup() %>% select(qtile,x,ci95_lb,ci95_ub) %>% 
  mutate(stat="Q 0.75", x=paste0("x=",x)) %>% rename(val=qtile)

combined<-rbind(post_mn_2,post_q75_2) 

pdf(file.path(plotdir,"mean_qtile2.pdf"), height=3.5, width=6)
combined %>% 
  ggplot(aes(x=val))+geom_histogram(bins=40) +
  geom_vline(aes(xintercept=ci95_lb), lty=6, color="forestgreen", data=post_mn_2)+
  geom_vline(aes(xintercept=ci95_ub), lty=6, color="forestgreen", data=post_mn_2)+
  geom_vline(aes(xintercept=ci95_lb), lty=6, color="red", data=post_q75_2)+
  geom_vline(aes(xintercept=ci95_ub), lty=6, color="red", data=post_q75_2)+
  xlab("")+
  facet_grid(stat~x, labeller="label_value")+
  theme(text=element_text(size=20), axis.text = element_text(size=20))
dev.off()



post_q75 %>% group_by(x) %>% select(x,qt_95ci_lb) %>% table()
post_q75 %>% group_by(x) %>% select(x,qt_95ci_ub) %>% table()

cdf_logistic_dat50levs_logis

## - model misspecification

if(0){
fitmod("probit", "dat25levs_norm"); attr(summary(fit_probit_dat25levs_norm),"family")
fitmod("logistic", "dat25levs_norm")
fitmod("loglog", "dat25levs_norm")
fitmod("cloglog", "dat25levs_norm", refresh=50)
}

nsamptag<-"dat50levs_"

# couldn't get cloglog link to work

lnks<-c("probit","logistic","loglog")
datasets<-paste0(nsamptag,c("norm","logis","ev2"))

grd<-expand.grid(lnks,datasets, stringsAsFactors = FALSE)

# run for all combinations of true model and fit
all_fits<-mapply(fitmod, grd[,1], grd[,2])

all_fits_nms<-paste0("fit_",grd[,1],"_",grd[,2])
mod_nms<-paste0("mod_",grd[,1],"_",grd[,2])
  
all_cdfs<-lapply(all_fits_nms, function(x)
  getCDF(get(x), newdata=data.frame(z1=c(0,1))) )

names(all_cdfs)<-mod_nms

# combine & add name as variable
# make fit_link and true_link vars
combined_cdfs <- do.call(rbind,all_cdfs) %>% rownames_to_column() %>%
  separate(rowname,c("mod","fit_link","sampsize","true_link","num")) %>%
  mutate(true_link=factor(true_link,levels=c("norm","logis","ev2")),
         fit_link=factor(fit_link,levels=c("probit","logistic","loglog")))

combined_truth<-
  rbind(  
    rbind(data.frame(n_cdf_0,z1=0), data.frame(n_cdf_1,z1=1)) %>% mutate(true_link=rep("norm",202)), 
    rbind(data.frame(l_cdf_0,z1=0), data.frame(l_cdf_1,z1=1)) %>% mutate(true_link=rep("logis",202)), 
    rbind(data.frame(ev2_cdf_0,z1=0), data.frame(ev2_cdf_1,z1=1)) %>% mutate(true_link=rep("ev2",202)) 
  ) %>%
  mutate(true_link=factor(true_link,levels=c("norm","logis","ev2")))

all_Means<-lapply(all_fits_nms, function(x)
  getMean(get(x), newdata=data.frame(z1=c(0,1))) )
names(all_Means)<-mod_nms

combined_Means<-do.call(rbind,all_Means) %>% rownames_to_column() %>%
  separate(rowname,c("mod","fit_link","sampsize","true_link","num")) %>%
  mutate(mn_plt=ifelse(z1==0,
           paste0("E[Y|x=0]=",round(mean_mn,3)),
           paste0("E[Y|x=1]=",round(mean_mn,3)))) %>%
  mutate(true_link=factor(true_link,levels=c("norm","logis","ev2")),
         fit_link=factor(fit_link,levels=c("probit","logistic","loglog")))


# plot

#single row for slides
combined_cdfs1<-combined_cdfs %>% filter(true_link=="norm")
combined_truth1<-combined_truth %>% filter(true_link=="norm")
combined_Means1<-combined_Means %>% filter(true_link=="norm")

pdf(file.path(plotdir,"misspec_plot1.pdf"), width=10.5, height=5)
combined_cdfs1 %>% ggplot(aes(group=z1)) +
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5), stat="stepribbon", direction="hv",fill="grey30", alpha=0.4) +
  geom_step(aes(x=yval,y=med_cdf)) +
  geom_line(data=combined_truth1,aes(x=yval,y=cdf, group=z1),color="blue",inherit.aes = FALSE)+ 
  geom_text(data=combined_Means1,aes(x=rep(-3.7,6),y=rep(c(0.95,0.85),3),label=
                                      mn_plt )) +
  xlab("") + ylab("Conditional CDF") +
  geom_segment(data=combined_Means1, aes(x=mn_q2.5 ,y=-0.07, xend=mn_q97.5, yend=-0.07),
               lty=1,color="red") +
  annotate("text",x=-4.1,y=-0.07, label="95% PIs for E[Y|x]") +
  facet_grid(true_link~fit_link)
dev.off()

#matrix for poster
pdf(file.path(plotdir,"misspec_plot2.pdf"), width=10, height=6)
combined_cdfs %>% ggplot(aes(group=z1)) +
  geom_ribbon(aes(x=yval, ymin=cdf_q2.5,ymax=cdf_q97.5), stat="stepribbon", direction="hv",fill="grey30", alpha=0.4) +
  geom_step(aes(x=yval,y=med_cdf)) +
  geom_line(data=combined_truth,aes(x=yval,y=cdf, group=z1),color="blue",inherit.aes = FALSE)+ 
  geom_text(data=combined_Means,aes(x=rep(-3.8,nrow(grd)*2),y=rep(c(0.95,0.8),nrow(grd)),label=
                                      mn_plt )) +
  xlab("") + ylab("Conditional CDF") +
  geom_segment(data=combined_Means, aes(x=mn_q2.5 ,y=-0.07, xend=mn_q97.5, yend=-0.07),
               lty=1,color="red") +
  annotate("text",x=-4.1,y=-0.07, label="E[Y|x] 95% PIs") +
  facet_grid(true_link~fit_link)
dev.off()



## addtl checks/code
if(0){
# is fit link is exactly the same?
combined_cdfs %>% filter(true_link=="norm") %>% 
  ggplot(aes(group=fit_link)) +
  geom_step(aes(x=yval,y=med_cdf)) + facet_grid(~z1)

combined_cdfs %>% filter(true_link=="logis") %>% 
  ggplot(aes(group=fit_link)) +
  geom_step(aes(x=yval,y=med_cdf)) + facet_grid(~z1)

combined_cdfs %>% filter(true_link=="ev2") %>% 
  ggplot(aes(group=fit_link)) +
  geom_step(aes(x=yval,y=med_cdf)) + facet_grid(~z1)
  
Fnz1_0 <- dat40levs_logis %>% filter(z1==0) %>% pull(y) %>% ecdf()
Fnz1_1 <- dat40levs_logis %>% filter(z1==1) %>% pull(y) %>% ecdf()

plot(Fnz1_0);curve(plogis(x),-6,5, add=TRUE, col="blue")
plot(Fnz1_1);curve(plogis(x,bta,1),-4,8, add=TRUE, col="blue")

}


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
