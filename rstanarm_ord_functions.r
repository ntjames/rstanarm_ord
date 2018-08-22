# rm(list=ls())

library(rstanarm)
library(rms)
library(dplyr)

# modified from Liu et al. sim code
generate.data.2 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  y <- rnorm(n, alpha+beta[1]*z1, sigma)
  data <- data.frame(y=y, z1=z1)
  return(data)
}

# call this once to distribute chains across cpu cores:
options(mc.cores=parallel::detectCores())

# get inverse link function ??
# for logistic use plogis()

if (0){
inv_logit<-function(x){
  exp(x)/(1+exp(x))
}
  
logit<-function(p){
  log(p/(1-p))
}

logit(0.3)

inv_logit(-0.8472979)
plogis(-0.8472979)
}

## -- estimate conditional CDF -- ##

getCDF <- function(spolrfit, newdata, summ=TRUE,...){
  require(dplyr)
  require(stringr)
  
  #check that cumulative model used
  #check for brmfit and newdata
  #check newdata is a data.frame
  # check that names in newdata match coefs from model
  if( !identical(sort(names(coef(spolrfit))), 
                sort(names(newdata))) ) stop("newdata vars must match model")
  
  # other checks?
  
  nsamps<-attr(summary(spolrfit),"posterior_sample_size")
  
  # get values of outcome from ordered factor to numeric
  
  #!old don't need real name of outcome, fit object always calls it y
  #fmla<-attr(summary(spolrfit),"formula")
  #outcome<-str_remove(as.character(fmla[2]),"ordered") %>% str_sub(2,-2)
  #truey0 <- as.numeric( levels( unlist(spolrfit[outcome]) ) ) %>% sort()
  
  truey0 <- as.numeric( levels(spolrfit$y) ) %>% sort()
  
  # prepend value less than min(y) for alpha_0=-Inf intercept
  truey<-c(-Inf,truey0) 
  
  # format newdata, betas, and intercepts
  cv_nms<-names(coef(spolrfit)) 
  ndr <- newdata %>% mutate(ndrow=1:n())
  nd <- ndr %>% select(-ndrow) %>% as.matrix()
 
  beta <- as.data.frame(spolrfit) %>% select(cv_nms) %>% as.matrix()
  int <- as.data.frame(spolrfit) %>% select(-cv_nms) %>% as.matrix()
  
  # get matrix of linear predictions Xb
  # (rxp)x (pxs) = rxs
  # r is rows in newdata, p is parameters (cols) in newdata, 
  # s is number of MCMC samples
  Xb <- nd %*% t(beta) 
  
  # add Xb to each intercept (4000xints) 
  #dim(int) => s x (ints-1)
  #dim(Xb) => r x s
  
  #use inverse function based on family
  #! add cauchit
  fam <- spolrfit$family
  if (fam=="probit") {
    inv_func <- pnorm
  } else if (fam == "logistic") {
    inv_func <- plogis
  } else if (fam == "loglog") {
    inv_func <- function(y) exp(-exp(-y))
  } else if (fam == "cloglog") {
    inv_func <- function(y) 1-exp(-exp(y))
  } else if (fam == "cauchit") {
    inv_func <- pcauchy #! not sure if this is right
  }
  
  #will have 1 for each row of nd
  # check model/doc to make sure values are being calculated correctly
  # are cutpoints y<= or y< ??
  for (i in 1:nrow(nd)){
      tmpcdf0 <- int - t(Xb[rep(i,ncol(int)),, drop=FALSE]) 
      tmpcdf1 <- cbind(`-Inf`=-Inf, tmpcdf0, `Inf`=Inf) # add alpha_0=-Inf and alpha_n = Inf
      tmpcdf <- tmpcdf1 %>% as.data.frame.table() %>% 
        mutate(cdf=inv_func(Freq), ndrow=i) %>%
        cbind(nd[i,,drop=FALSE])
      assign(paste0("cc",i), tmpcdf)
  }
   
#  F(y_1|X)=G^-1(alpha_i-betaX)
  
# combine conditional cdfs
nd_ds<-ls()[grep("cc",ls(),fixed=TRUE)] # list of all conditional cdf datasets
cdf_vals<-do.call(rbind, lapply(nd_ds, function(x) get(as.character(x))))

  if (summ){
    cdf_summ<-cdf_vals %>%
      ungroup() %>%
      group_by(ndrow, Var2) %>%
      dplyr::summarize(mn_cdf=mean(cdf),
                       med_cdf=median(cdf),
                       cdf_q2.5=quantile(cdf,probs=0.025),
                       cdf_q97.5=quantile(cdf,probs=0.975)) %>%
      ungroup() %>% mutate(yval=rep(truey,nrow(nd))) %>%  
      full_join(., ndr, by="ndrow")
    return(cdf_summ)
  } else {
   cdf_out <- cdf_vals %>%
      ungroup() %>%
      dplyr::arrange(ndrow, Var1) %>%
      mutate(yval=rep(truey,nrow(nd)*nsamps  ))
   return(cdf_out)
  }

}

#scratch
if (0){
#! need to get posterior dist. of the parameters (z1 & all intercepts)
# not just the linear predictor or outcome

# no. draws x 1 vector
posterior_linpred(rstanarm_fit_logit,newdata=data.frame(z1=0)) %>% dim()
posterior_predict(rstanarm_fit_logit,newdata=data.frame(z1=0)) %>% dim()

predict(rstanarm_fit_logit) %>% dim()

# gives err
predict(rstanarm_fit_logit, newdata=data.frame(z1=0)) %>% dim()

posterior_interval(rstanarm_fit_logit,newdata=data.frame(z1=0))

fitted(rstanarm_fit_logit,newdata=data.frame(z1=0))

# this is closest to what we want, but can't take newdata,
# have to manually calculate int + Xbeta 
as.data.frame(rstanarm_fit_logit) %>% dim()
dat25levs



bb1<-data.frame(a=c(1,2),b=c(3,4))
bb2<-data.frame(a=c(4,2,5),b=c(9,3,4))
rbind(bb1,bb2)

nd_ds<-ls()[grep("^bb*",ls())]

do.call(rbind, lapply(nd_ds, get) )

lapply(ls()[grep("^bb*",ls())], get)

fmla<-attr(summary(rstanarm_fit_probit2),"formula")
outcome<-str_remove(as.character(fmla[2]),"ordered") %>% str_sub(2,-2)

}

# manual coding for getCDF() function
if(0){
beta<-as.data.frame(rstanarm_fit_logit) %>% pull(z1) %>% as.matrix()
int<-as.data.frame(rstanarm_fit_logit) %>% select(-z1) %>% as.matrix()
nd<-data.frame(z1=c(0,1)) %>% as.matrix()

# (rxp)x (px4000) = rx4000

Xb <- nd %*% t(beta) 

# add Xb to each intercept (4000xints) 
#dim(int) => 4000 24
#dim(Xb) => 2 4000

#will have 1 for each row of nd
# add loop??
cdf0 <- int - t(Xb[rep(1,24),, drop=FALSE]) 
cdf1 <- int - t(Xb[rep(2,24),, drop=FALSE]) 

# add alpha_0=-Inf and alpha_n = Inf
cdf0_add<-cbind(`-Inf`=-Inf,cdf0, `Inf`=Inf)
cdf1_add<-cbind(-Inf,cdf1, Inf)


# need some way to combine cdf0, cdf1, etc. and add covar info

cc0<-cdf0_add %>% as.data.frame.table() %>% mutate(cdf=plogis(Freq)) 
cc1<-cdf1_add %>% as.data.frame.table() %>% mutate(cdf=plogis(Freq)) 

cc0 %>% arrange(Var1) %>% filter(Var1=="A")

#go into guts of fit object ??

rstanarm_fit_logit[["stanfit"]]@sim[["samples"]][[1]][["zeta[1]"]]
rstanarm_fit_logit[["stanfit"]]@sim[["samples"]][[1]][["zeta[2]"]]
class(rstanarm_fit_logit[["stanfit"]]@sim[["samples"]][[1]][["zeta[1]"]])

get_zeta<-function(i){
  zt<-paste0("zeta[",i,"]")
  rstanarm_fit_logit[["stanfit"]]@sim[["samples"]][[1]][[zt]]
}

#samples from cdf for x=0
do.call(cbind, lapply(1:24, function(x) get_zeta(x))) %>% data.frame %>% head() %>%
  mutate_all(funs(plogis))


#match coefs in sample to variable name
#i.e. beta for z1 is rstanarm_fit_logit[["stanfit"]]@sim[["samples"]][[1]][["beta[1]"]]


as.data.frame.table(fitvals)

}


## -- estimate conditional mean -- ##

getMean<-function(spolrfit,newdata,summ=TRUE,...){
  require(dplyr)
  
  cdf_samps <- getCDF(spolrfit, newdata, summ=FALSE)
  
  mn_vals<- cdf_samps %>% filter(cdf!=0) %>% group_by(ndrow, Var1) %>%
    mutate(pdf0=lag(cdf), pdf=ifelse(is.na(pdf0),cdf,cdf-pdf0), fy_Py=pdf*yval) %>% 
    dplyr::summarize(n=n(),mn=sum(fy_Py)) %>% ungroup()
  
  #! pass this from getCDF rather than recreating 
  ndr<- newdata %>% mutate(ndrow=1:n())
  
  if (summ){
    mn_summ<-mn_vals %>%
      ungroup() %>%
      group_by(ndrow) %>%
      dplyr::summarize(mean_mn=mean(mn),
                       med_mn=median(mn),
                       sd_mn=sd(mn),
                       mn_q2.5=quantile(mn,probs=0.025),
                       mn_q97.5=quantile(mn,probs=0.975)) %>%
      full_join(., ndr, by="ndrow")
    return(mn_summ)
  } else {
    #! append variable info
    return(mn_vals)
  }
  
}

#scratch
if(0) {
d<-cond_cdf_samps %>% filter(cdf!=0) %>%
  group_by(ndrow, Var1) %>%
  mutate(pdf0=lag(cdf), pdf=ifelse(is.na(pdf0),cdf,cdf-pdf0), fy_Py=pdf*yval) %>% 
  dplyr::summarize(n=n(),mn=sum(fy_Py))  %>% ungroup()
}


## -- estimate conditional quantiles -- ##

getQuantile<-function(spolrfit,newdata,q,summ=TRUE,...){
  require(dplyr)
  
  cdf_samps <- getCDF(spolrfit, newdata, summ=FALSE)

  #! pass these rather than recreating
  
  #! check if truey is correct
  # need to deal with quantiles below lowest obs and greater than highest obs
  truey0 <- cdf_samps %>% filter(ndrow==1 & Var1=="A") %>% select(yval) %>% pull()
  truey<-c(truey0[-1],Inf)
  
  ndr <- newdata %>% mutate(ndrow=1:n())
  cv_nms <- names(newdata)
  
  qtile_vals <- cdf_samps %>% group_by(ndrow, Var1) %>% 
    mutate( idx.1 = max(which(cdf<=q)), idx.2 = min(which(cdf>=q)),
            cdf.1 = cdf[idx.1], cdf.2 = cdf[idx.2], n=1:n()) %>%
    filter(n==1) %>%
    mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
           idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2), 
           idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
           idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
           qtile=ifelse(idx.1==idx.2,idx.y1,
                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1)) %>%
          ungroup()
  
  
  if (summ){
    qtile_summ <- qtile_vals %>%
      group_by(ndrow) %>% 
      dplyr::summarize(mean_qtile=mean(qtile),
                       med_qtile=median(qtile),
                       sd_qtile=sd(qtile),
                       qtile_q2.5=quantile(qtile,probs=0.025),
                       qtile_q97.5=quantile(qtile,probs=0.975)) %>%
      full_join(., ndr, by="ndrow")
    return(qtile_summ)
  } else {
    #! append nd vars, keep other intermediary vars??
    
    qtile_vals_out<-qtile_vals %>% select(ndrow, cv_nms, 
                                          idx.y1.cdf, idx.y2.cdf,
                                          idx.1, idx.2, qtile) 
    return(qtile_vals_out)
  }
  
}

# scratch
if (0){

  q<-0.5
  
  truey<-cdf_samps %>% filter(ndrow==1 & Var1=="A") %>% select(yval) %>% pull()
  
  
  
  foo<-cdf_samps %>% group_by(ndrow, Var1) %>% 
    mutate( idx.1 = max(which(cdf<=q)), idx.2 = min(which(cdf>=q)),
            cdf.1 = cdf[idx.1], cdf.2 = cdf[idx.2], n=1:n()) %>%
    filter(n==1) %>%
    mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
           idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2), 
           idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
           idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
           qtile=ifelse(idx.1==idx.2,idx.y1,
                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1))  %>%
    ungroup()
  
  
  foo %>%
    group_by(ndrow) %>% 
    dplyr::summarize(mean_qtile=mean(qtile),
                     med_qtile=median(qtile),
                     sd_qtile=sd(qtile),
                     qtile_q2.5=quantile(qtile,probs=0.025),
                     qtile_q97.5=quantile(qtile,probs=0.975)) 
    
  cv_nms<-"z1"
  
  foo %>% select(ndrow, cv_nms, idx.y1.cdf , idx.y2.cdf, idx.1, idx.2, qtile) %>% head()
  
  qtile_vals<-fv_tab %>% group_by(Var2, Var1) %>% 
    mutate(cdf=cumsum(Freq),
           idx.1 = max(which(cdf<=q)),
           idx.2 = min(which(cdf>=q)),
           cdf.1 = cdf[idx.1],
           cdf.2 = cdf[idx.2]) %>%
    filter(Var3=="A") %>%
    mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
           idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2), 
           idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
           idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
           qtile=ifelse(idx.1==idx.2,idx.y1,
                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1)) 
  
  if (summ){
    qtile_summ <- qtile_vals %>%
      group_by(Var2) %>% 
      dplyr::summarize(mean_qtile=mean(qtile),
                       med_qtile=median(qtile),
                       sd_qtile=sd(qtile),
                       qtile_q2.5=quantile(qtile,probs=0.025),
                       qtile_q97.5=quantile(qtile,probs=0.975)) %>%
      full_join(., nd, by="Var2")
    return(qtile_summ)
  } else {
    #! append nd vars, keep other intermediary vars??
    qtile_vals_out<- qtile_vals %>% select(Var1, Var2, Var3, qtile)
    return(qtile_vals_out)
  }
  
  
  
}


## -- Exceedance Probability -- ##

if(0){
?rms::ExProb

set.seed(1)
x1 <- runif(100)
yvar <- x1 + runif(100)
f <- orm(yvar ~ x1)
d <- ExProb(f)
lp <- predict(f, newdata=data.frame(x1=c(.2,.8)))
w <- d(lp)

w

s1 <- abs(x1 - .2) < .1
s2 <- abs(x1 - .8) < .1
plot(w, data=data.frame(x1=c(rep(.2, sum(s1)), rep(.8, sum(s2))),
                        yvar=c(yvar[s1], yvar[s2])))
}
