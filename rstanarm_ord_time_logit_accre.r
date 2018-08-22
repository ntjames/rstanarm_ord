Sys.time()

rm(list=ls())

# define library directory on ACCRE
.libPaths("~/R/rlib_3.4.3")

library(rstanarm)

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

# get number of levels from command line
arg<-commandArgs(trailing=TRUE)
lev<-as.integer(arg[1])

# make dataset
assign(paste0("dat_rstanarm_",lev), generate.data.1(seed=47232, n=lev))

# run model and save fit and time
#times<-system.time( fit<-stan_polr(ordered(y) ~ z1+z2, data = get(paste0("dat_rstanarm_",lev)),
#                          prior=NULL, prior_counts = dirichlet(1), method="logistic") )

# don't save fit, just time
times<-system.time( stan_polr(ordered(y) ~ z1+z2, data = get(paste0("dat_rstanarm_",lev)),
                          prior=NULL, prior_counts = dirichlet(1), method="logistic",
                        do_residuals = FALSE) )
times

timetag <- format(Sys.time(), "_%Y%m%d_%H%M_%Z")

# uncomment to save fits
# assign(paste0("mod_rstanarm_",lev), fit)
assign(paste0("time_rstanarm",timetag,"_",lev), times)

#rm(fit, arg, times)
rm(arg, times)

save.image(file=paste0("rstanarm_ord_sim_logit_",lev,"levs",timetag, ".RData"))

Sys.time()
