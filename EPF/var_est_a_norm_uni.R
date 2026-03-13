
# setwd("D:/Dropbox/_My_Doc/Regression/EPF/")
setwd("/home/docsa/Dropbox/_My_Doc/Regression/EPF")



rm(list = ls()) 

source("moment_epf.R")
source("a_pmm3ex.R")
source("min_var_est_a_norm_uni.R")

p_start <- 0.05 # pochatok p
p_delta <- 0.05 # krok p
n_start <- 20 # pochatok n
n_delta <- 20 # krok n

i_p<-20
j_n<-25
w<-10000

ZZ <- array(0, dim=c(i_p,j_n))


for (i in 1:i_p)
{ 
  for (j in 1:j_n)
  {
    ZZ[i, j] <- min_var_est_a_norm_uni(p_start + p_delta*(i-1), n_start + n_delta*(j-1), w)
  }
}

ZZ

