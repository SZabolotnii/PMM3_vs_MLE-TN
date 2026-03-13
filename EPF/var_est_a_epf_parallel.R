# install.packages("normalp")
install.packages("parallel")

library(normalp)
library(parallel)

# library(L1pack)

setwd("ZSV/EPF/")


rm(list = ls()) 

source("moment_epf.R")
source("a_pmm3ex.R")
source("min_var_est_a_epf.R")
source("min_var_est_a_epf_apply.R")


w<-10000

x.list<- list(w)

reg_ef <- function(x) {
  p_start <- 1.5 # pochatok p
  p_delta <- 0.05 # krok p
  n_start <- 20 # pochatok n
  n_delta <- 20 # krok n
  i_p<-20
  j_n<-15
  ZZ <- array(0, dim=c(i_p,j_n))
  for (i in 1:i_p)
  { 
    for (j in 1:j_n)
    {
      ZZ[i, j] <- min_var_est_a_epf_apply(p_start+p_delta*(i-1), n_start+n_delta*(j-1), x)
    }
  }
  return(ZZ)
}


# detect the number of cores
n.cores <- detectCores()

# single core
# system.time(a <- lapply(x.list, reg_ef))
# a

system.time(aa <- mclapply(x.list, reg_ef, mc.cores = n.cores))
aa

write.csv(aa, file = "obl_ef.csv")
