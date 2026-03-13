# install.packages("normalp")
library(normalp)
# library(L1pack)

setwd("D:/Dropbox/_My_Doc/Regression/EPF/")

rm(list = ls()) 

source("moment_epf.R")
source("a_pmm3ex.R")

n <- 100
w <- 10000
A.t <- 1

AA.mean <- A.t
AA.median <- A.t
AA.midrange <- A.t
AA.mle <- A.t
AA.pmm3 <- A.t

p_epf <- 5

for (j in 1:w)
{ 
 
  err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
  X <- A.t + err
  A.mean <- mean(X)
  res <- X - A.mean
  A.median <- median(X)
  A.midrange <- (max(X)+min(X))/2
  # MLE <- paramp(X, p=p_epf)
  MLE <- paramp(X)
  A.mle <- as.numeric(MLE["mp"])
  A.pmm3 <- a_pmm3ex(X)
  
    
  AA.mean <- rbind(AA.mean, A.mean)
  AA.median <- rbind(AA.median, A.median)
  AA.midrange <- rbind(AA.midrange, A.midrange)
  AA.mle <- rbind(AA.mle, A.mle)
  AA.pmm3 <- rbind(AA.pmm3, A.pmm3)
}

AA.mean <- AA.mean[1:w+1]
AA.median <- AA.median[1:w+1]
AA.midrange <- AA.midrange[1:w+1]
AA.mle <- AA.mle[1:w+1]
AA.pmm3 <- AA.pmm3[1:w+1]

V.mean <- var(AA.mean)
V.median <- var(AA.median)
V.midrange <- var(AA.midrange)
V.mle <- var(AA.mle)
V.pmm3 <- var(AA.pmm3)

M.mean <- mean(AA.mean)
M.median <- mean(AA.median)
M.midrange <- mean(AA.midrange)
M.mle <- mean(AA.mle)
M.pmm3 <- mean(AA.pmm3)


g3ef.e.mean <- V.pmm3/V.mean
g3ef.e.mle <- V.pmm3/V.mle






V.pmm3/V.mean
V.pmm3/V.median
V.pmm3/V.midrange
V.pmm3/V.mle


g4.t <- moment.epf(4,p_epf,1)/moment.epf(2,p_epf,1)^2-3
g6.t <- moment.epf(6,p_epf,1)/moment.epf(2,p_epf,1)^3-15*(g4.t+1)

g3ef.t <- 1-g4.t^2/(6+9*g4.t+g6.t)

x <- c(-1, 1)
# f <- dnormp(x, p=p_epf)
# print(f)
plot(function(x) dnormp(x, p=p_epf) , -4, 4,
     main = "Exponential power distribution density function", ylab="f(x)")

boxplot(AA.mean, AA.median, AA.midrange, AA.mle, AA.pmm3, col = "lightgray", horizontal = TRUE)
