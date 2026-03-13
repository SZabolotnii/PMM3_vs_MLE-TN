
library(GLDEX)
library(nleqslv)

setwd("ZSV/EPF/")

rm(list = ls())

source("solvNL2.R")
source("solvNL2_nleqslv.R")

n <- 100
X <- runif(n, min = 0, max = 20)
a0 <- 10
b0 <- -0.6
err <- (rchisq(n, 3, ncp = 0)-3)/sqrt(3)
# mean(err)
# var(err)
# skewness(err)
# kurtosis(err)
y<-a0*exp(b0*X) + err

nonlin_mod=nls(y~a*exp(b*X),start=list(a=5,b=-1))

B.ols <- coef(nonlin_mod)

res <- residuals(nonlin_mod)
k1 <- mean(res)
k2 <- var(res)
g3 <- skewness(res)
g4 <- kurtosis(res)

# summary(nonlin_mod)
B.ols
solvNL2(B.ols,X,y,k2,g3,g4)
B.pmm2 <- solvNL2_nleqslv(B.ols,X,y,k2,g3,g4)
B.pmm2$x

# plot(X,y)
#
# plot(X,predict(nonlin_mod))
# plot(X,residuals(nonlin_mod))


n <- 50 #!!!!!!!!!!!!!!!!!!!!!!
w <- 10000
B.t <- c(10, -0.6)
L <- length(B.t)



BB.ols <- B.t
BB.pmm2 <- B.t


for (j in 1:w)
{

  X <- runif(n, min = 0, max = 20)
  err <- (rchisq(n, 3, ncp = 0)-3)/sqrt(3)

  Y.t <- B.t[1]*exp(B.t[2]*X)

  Y <- Y.t+err
  nonlin_mod=nls(Y~a*exp(b*X),start=list(a=11,b=-0.5))
  
  B.ols <- coef(nonlin_mod)
  BB.ols <- rbind(BB.ols, B.ols)
  
  res <- residuals(nonlin_mod)
  
  k2 <- var(res)
  g3 <- skewness(res)
  g4 <- kurtosis(res)
  

  # B.pmm2 <- solvNL2(B.ols,X,Y,k2,g3,g4)
  B.pmm2 <- solvNL2_nleqslv(B.ols,X,Y,k2,g3,g4)
  BB.pmm2 <- rbind(BB.pmm2, B.pmm2)

}

BB.ols <- BB.ols[1:w+1,1:L]
BB.pmm2 <- BB.pmm2[1:w+1,1:L]



V.ols <- var(BB.ols)
V.pmm2 <- var(BB.pmm2)


M.ols <- rbind(mean(BB.ols[,1:1]), mean(BB.ols[,2:2]))
M.pmm2 <- rbind(mean(BB.pmm2[,1:1]), mean(BB.pmm2[,2:2]))


g2ef.t <- 1-g3^2/(g4+2)

g2ef.e <- V.pmm2/V.ols


hist(BB.pmm2[,1])





boxplot(BB.ols[,1:1], BB.pmm2[,1:1], col = "lightgray", horizontal = FALSE)
boxplot(BB.ols[,2:2], BB.pmm2[,2:2], col = "lightgray", horizontal = FALSE)
# 
# 
# 
# 
