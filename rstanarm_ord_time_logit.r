rm(list=ls())
library(tidyverse)
library(rstanarm)
library(ggplot2)

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
sessionInfo()
parallel::detectCores()

# vary n from 200 to 10000
#levs<-2^(1:8)*100

# vary levels of y from 200 to 2000
levs<-seq(200, 2000,by=100)
times<-vector("list",length(levs))

#!! takes a while to run
for(i in seq_along(levs)){

assign(paste0("dat",i), generate.data.1(seed=47232, n=levs[i]))
times[[i]]<-system.time( assign(paste0("mod",i), stan_polr(ordered(y) ~ z1+z2, data = get(paste0("dat",i)),
                                 prior=NULL, prior_counts = dirichlet(1), method="logistic") ) )

# lm for debugging
# times[[i]]<-system.time(  assign(paste0("mod",i), lm(y~z1+z2,data=get(paste0("dat",i))) ) )
}

timetag <- format(Sys.time(), "_%Y%m%d_%H%M_%Z")
save.image(file=paste0("rstanarm_ord_sim_logit", timetag, ".RData"))

#pdf("rstanarm_ord_logit.pdf")
#do.call(rbind,times) %>% as_tibble() %>% mutate(levels=levs)%>%
#  ggplot(aes(x=levels,y=elapsed))+geom_smooth(method="lm",formula=y~x+I(x^2),se=FALSE)+
#  geom_smooth(method="lm",formula=y~x,col="red",se=FALSE)+geom_point()
#dev.off()
