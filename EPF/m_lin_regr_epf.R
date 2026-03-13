library(normalp)
library(L1pack)

setwd("D:/Dropbox/_My_Doc/R/EPF/")

rm(list = ls())

source("moment_epf.R")
source("solvMRs3ex.R")
source("step_pol.R")


n <- 200 #!!!!!!!!!!!!!!!!!!!!!!
w <- 1000
B.t <- c(1, 3, -1)
L <- length(B.t)


BB.ols <- B.t
BB.mle <- B.t
BB.lad <- B.t
BB.pmm3 <- B.t

p_epf <- 6.4 #!!!!!!!!!!!!!!!!!!!!

m2.t<-moment.epf(2,p_epf,1)
m4.t<-moment.epf(4,p_epf,1)
m6.t<-moment.epf(6,p_epf,1)

X0 <- replicate(n, 1)
  
for (j in 1:w)
{ 
  
  X1 <- runif(n, min = 0, max = 1)
  X2 <- runif(n, min = 0, max = 1)
  
  
  err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
    
  Y.t <-  B.t[1]+B.t[2]*X1+B.t[3]*X2

  Y <- Y.t+err
  
  reg.ols <- lm(Y ~ X1 + X2)
  # reg.ols <- lm(Y ~ X1 + I(X1^2)) # добвка "+ I(X^2)" це для квадратичної регресії (три коефіцієнти) !!!
  
  
  res <- residuals(reg.ols)
  B.ols <- coefficients(reg.ols)
  BB.ols <- rbind(BB.ols, B.ols)

  
  # reg.mle <- lmp(Y ~ X1 + X2, p=p_epf)
  reg.mle <- lmp(Y ~ X1 + X2)
  # reg.mle <- lmp(Y ~ X1 + I(X1^2))
  
  
  B.mle <- coefficients(reg.mle)
  BB.mle <- rbind(BB.mle, B.mle)
  
  # reg.lad <- lad(Y ~ X)
  # B.lad <- coefficients(reg.lad)
  # BB.lad <- rbind(BB.lad, B.lad)
  
  

  # m2 <-m2.t; m4 <- m4.t; m6 <- m6.t
  m2 <- mean(res^2); m4 <- mean(res^4); m6 <- mean(res^6)

  # g4 <- m4/m2^2-3; g6 <- m6/m2^3-15*(g4+1)

  A=1;
  B=-Y*3;
  C=3*Y^2-(m6-3*m4*m2)/(m4-3*m2^2);
  D=Y*(m6-3*m4*m2)/(m4-3*m2^2)-Y^3;

  B.pmm3 <- t(solvMRs3ex(B.ols,cbind(X0, X1, X2),Y,m2,m4,m6))
  # B.pmm3 <- t(solvMRs3ex(B.ols,cbind(X0, X1, X1^2),Y,m2,m4,m6))
  
  BB.pmm3 <- rbind(BB.pmm3, B.pmm3)
  
}

BB.ols <- BB.ols[1:w+1,1:L]
BB.pmm3 <- BB.pmm3[1:w+1,1:L]
BB.mle <- BB.mle[1:w+1,1:L]
# BB.lad <- BB.lad[1:w+1,1:L]


V.ols <- var(BB.ols)
V.pmm3 <- var(BB.pmm3)
V.mle <- var(BB.mle)
# V.lad <- var(BB.lad)


M.ols <- rbind(mean(BB.ols[,1:1]), mean(BB.ols[,2:2]))
M.pmm3 <- rbind(mean(BB.pmm3[,1:1]), mean(BB.pmm3[,2:2]))
M.mle <- mean(BB.mle[,1:1])
# M.lad <- mean(BB.lad[,1:1])


g4.t <- moment.epf(4,p_epf,1)/moment.epf(2,p_epf,1)^2-3
g6.t <- moment.epf(6,p_epf,1)/moment.epf(2,p_epf,1)^3-15*(g4.t+1)

g3ef.t <- 1-g4.t^2/(6+9*g4.t+g6.t)

g3.pmm3 <- V.pmm3/V.ols
# g3.lad  <- V.pmm3/V.lad
g3.mle <- V.mle/V.ols

(g3.mle[1,1]*g3.mle[2,2]*g3.mle[3,3])^(1/3)
(g3.pmm3[1,1]*g3.pmm3[2,2]*g3.pmm3[3,3])^(1/3)

g3ef.t

