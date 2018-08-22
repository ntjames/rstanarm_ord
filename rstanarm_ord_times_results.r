
rm(list=ls())
library(tidyverse)
library(rstanarm)
library(ggplot2)

plotdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/JSM2018/speed session/fig")

if(0){
# from Liu et al. sim code
generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}

# call this once to distribute chains across cpu cores:
options(mc.cores=parallel::detectCores())

# n=2000, vary levels of y from 200 to 2000
#levs<-seq(200,2000,by=100)
levs<-seq(200,500,by=100)
times<-vector("list",length(levs))

#!! takes a while to run
for(i in seq_along(levs)){
  
  assign(paste0("dat",i), generate.data.1(seed=47232, n=levs[i]))
  times[[i]]<-system.time( assign(paste0("mod",i), stan_polr(ordered(y) ~ z1+z2, data = get(paste0("dat",i)),
                                                             prior=NULL, prior_counts = dirichlet(1), method="logistic") ) )
  # lm for debugging
  # times[[i]]<-system.time(  assign(paste0("mod",i), lm(y~z1+z2,data=get(paste0("dat",i))) ) )
}


#! load times from sim on statcomp2
simtimes1<-do.call(rbind,times) %>% as_tibble() %>% mutate(levels=levs)

#simtimes1_0<-do.call(rbind,times) %>% as_tibble() %>% mutate(levels=levs[1:3])

save.image(file="rstanarm_ord_sim1_pc.RData")

load("/home/nathan/Dropbox/njames/school/PhD/ra/rstanarm/rstanarm_ord_sim1_pc.RData")
}


simlog<-readLines("/home/nathan/Dropbox/njames/school/PhD/ra/rstanarm/rstanarm_ord_time_logit_20180722_0842.out")

lns<-grep("Elapsed Time",simlog)

#foo[sort(c(lns,lns+1,lns+2))]
vals0<-simlog[lns+2] 
vals<- as.numeric( gsub("seconds (Total)"," ",vals0,fixed=TRUE) )
valmat<-matrix(vals, ncol=4, byrow=TRUE)

elapsed<-apply(valmat,1,mean)
levels<-seq(200,2000,by=100)[1:length(elapsed)]

simtimes_from_log<- data.frame(elapsed,levels) %>% as_tibble() 


rm(list=ls()[grep("mod",ls())])
rm(list=ls()[grep("dat",ls())])
rm(times)
rm(levs)
rm(i)

if(0){
load("/home/nathan/Dropbox/njames/school/PhD/ra/brms/brms_ord_sim3.RData")
simtimes3 <- do.call(rbind,times) %>% as_tibble() %>% mutate(levels=levs)

#extend model to look at more intercepts, see if trend continues to be linear
rm(list=ls()[grep("mod",ls())])
rm(list=ls()[grep("dat",ls())])
rm(list=c("times","levs","i"))
}

load("/home/nathan/Dropbox/njames/school/PhD/ra/brms/brms_ord_sim4.RData")
simtimes4<-do.call(rbind,times) %>% as_tibble() %>% mutate(levels=levs)

rm(list=ls()[grep("mod",ls())])
rm(list=ls()[grep("dat",ls())])

#rstanarm_times<-simtimes1_0 %>% mutate(package="rstanarm")

rstanarm_times<-simtimes_from_log %>% mutate(package="rstanarm") %>% select(elapsed, levels, package)
#brms_times<-rbind(simtimes3, simtimes4) %>% mutate(package="brms")
brms_times<- simtimes4 %>% mutate(package="brms") %>% select(elapsed, levels, package)


#rbind(rstanarm_times,brms_times)

pdf(file.path(plotdir,"comptime_out.pdf"),width=7, height=3.75)
rbind(rstanarm_times,brms_times) %>% filter(levels<=1500) %>%
  ggplot(aes(x=levels, y=elapsed, color=package))+
  geom_line() +
  geom_point() +
  xlab("# distinct y values")+
  ylab("computation time (seconds)")
dev.off()

#pdf(file.path(plotdir,"comptime.pdf"),width=7, height=4)
#rbind(rstanarm_times,brms_times) %>%
#  ggplot(aes(x=levels, y=elapsed, color=package))+
#  geom_line() +
#  geom_point() +
#  ylab("computation time (seconds)")
#dev.off()


#! add model description
# iterations, chains, model spec, etc