
rm(list=ls())
library(tidyverse)
library(rstanarm)
library(ggplot2)

plotdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/JSM2018/speed session/fig")
dir<-file.path(getwd(),"sim_results")

# rstanarm_accre sims
if(1){
objs0<-dir(file.path(dir),pattern="rstanarm_ord_sim_logit", full.names=TRUE)
objs<-objs0[grep("levs",objs0)]

for (i in seq(objs)){ load(objs[i]) }

times<-ls()[grep("time_rstanarm",ls())]
rstan_simtimes0<-do.call(rbind, lapply(times, get))

levs0<-gsub("time_rstanarm_","",times)
levs<-ifelse(nchar(levs0)>4,substring(levs0,19),levs0) %>% as.numeric()

rstanarm_times0<-rstan_simtimes0 %>% as.tibble() %>% mutate(levels=levs) %>% 
  select(elapsed, user.child, levels) %>% mutate(package="rstanarm")

rstanarm_times <- rstanarm_times0 %>% group_by(levels, package) %>% 
  summarize(elapsed=mean(elapsed), user.child=mean(user.child), nsims=n())

rm(list=ls()[grep("^time",ls())])
rm(list=ls()[grep("dat",ls())])
rm(i, lev, levs, objs, objs0, rstan_simtimes0)
}

# brms_accre sims
if (1){
objs0<-dir(file.path(dir),pattern="brms_ord_sim_logit", full.names=TRUE)
objs<-objs0[grep("levs",objs0)]

for (i in seq(objs)){ load(objs[i]) }

times<-ls()[grep("time_brms",ls())]
brms_simtimes0<-do.call(rbind, lapply(times, get))

levs0<-gsub("time_brms_","",times)
levs<-ifelse(nchar(levs0)>4,substring(levs0,20),levs0) %>% as.numeric()

brms_times0<-brms_simtimes0 %>% as.tibble() %>% mutate(levels=levs) %>% 
  select(elapsed, user.child, levels) %>% mutate(package="brms")

brms_times <- brms_times0 %>% group_by(levels, package) %>% 
  summarize(elapsed=mean(elapsed), user.child=mean(user.child), nsims=n()) 

rm(list=ls()[grep("^time",ls())])
rm(list=ls()[grep("dat",ls())])
rm(i, lev, levs, objs, objs0, brms_simtimes0)
}

pdf(file.path(plotdir,"comptime_accre.pdf"),width=7.5, height=3.25)
rbind(rstanarm_times,brms_times) %>%
  ggplot(aes(x=levels, y=elapsed, color=package))+
  geom_line() +
  geom_point() +
  xlab("# distinct y values")+
  ylab("computation time (seconds)")
dev.off()

# plot individual runs
rbind(rstanarm_times0,brms_times0) %>% group_by(levels) %>%
    mutate(run=1:n()) %>%
  ggplot(aes(x=levels, y=elapsed, color=package, group=run))+
  geom_line(alpha=0.5) +
  geom_point(alpha=0.75) +
  xlab("# distinct y values")+
  ylab("computation time (seconds)")

if(0){
# use user.child instead of elapsed?
# ! no, this is cumulative time including parallel jobs
rbind(rstanarm_times,brms_times) %>%
  ggplot(aes(x=levels, y=user.child, color=package))+
  geom_line() +
  geom_point() +
  xlab("# distinct y values")+
  ylab("computation time (seconds)")
}

#! add model description
# iterations, chains, model spec, etc