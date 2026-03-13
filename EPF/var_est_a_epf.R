# install.packages("normalp")
library(normalp)
# library(L1pack)

setwd("D:/Dropbox/_My_Doc/Regression/EPF/")



rm(list = ls()) 

source("moment_epf.R")
source("a_pmm3ex.R")
source("min_var_est_a_epf.R")

p_start <- 1.5 # pochatok p
p_delta <- 0.1 # krok p
n_start <- 20 # pochatok n
n_delta <- 20 # krok n

i_p<-20
j_n<-10
w<-10000

ZZ <- array(0, dim=c(i_p,j_n))


for (i in 1:i_p)
{ 
  for (j in 1:j_n)
  {
    ZZ[i, j] <- min_var_est_a_epf(p_start+p_delta*(i-1), n_start+n_delta*(j-1), w)
  }
}

ZZ

