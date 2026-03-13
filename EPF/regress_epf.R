library(normalp)
library(L1pack)

setwd("D:/Dropbox/_My_Doc/Regression/EPF/")

rm(list = ls())

source("moment_epf.R")
source("solvMRs3ex.R")
source("step_pol.R")


n <- 500 #!!!!!!!!!!!!!!!!!!!!!!
w <- 10000
B.t <- c(-1, 2)
L <- length(B.t)


BB.ols <- B.t
BB.mle <- B.t
BB.lad <- B.t
BB.pmm3 <- B.t

p_epf <- 5 #!!!!!!!!!!!!!!!!!!!!

for (j in 1:w)
{ 
  
  X <- runif(n, min = 0, max = 1)
    err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
    
  Y.t <-  step_pol(B.t,X)

  Y <- Y.t+err
  
  ### reg.ols <- lm(Y ~ X + I(X^2)) # добвка "+ I(X^2)" це для квадратичної регресії (три коефіцієнти) !!!
  reg.ols <- lm(Y ~ X)
  
  res <- residuals(reg.ols)
  B.ols <- coefficients(reg.ols)
  BB.ols <- rbind(BB.ols, B.ols)

  
  # reg.mle <- lmp(Y ~ X, p=p_epf)
  reg.mle <- lmp(Y ~ X)
  
  B.mle <- coefficients(reg.mle)
  BB.mle <- rbind(BB.mle, B.mle)
  
  # reg.lad <- lad(Y ~ X)
  # B.lad <- coefficients(reg.lad)
  # BB.lad <- rbind(BB.lad, B.lad)
  
  
  m2 <- mean(res^2); m4 <- mean(res^4); m6 <- mean(res^6)

  g4 <- m4/m2^2-3; g6 <- m6/m2^3-15*(g4+1)

  A=1;
  B=-Y*3;
  C=3*Y^2-(m6-3*m4*m2)/(m4-3*m2^2);
  D=Y*(m6-3*m4*m2)/(m4-3*m2^2)-Y^3;

  B.pmm3 <- t(solvMRs3ex(B.ols,cbind(X^0, X),Y,m2,m4,m6))
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

g3.ols <- V.pmm3/V.ols

# g3.lad  <- V.pmm3/V.lad

g3.mle <- V.pmm3/V.mle

(g3.mle[1,1]*g3.mle[2,2])^(1/2)
(g3.ols[1,1]*g3.ols[2,2])^(1/2)
g3ef.t












#(g3.mle[1,1]*g3.mle[2,2]*g3.mle[3,3])^(1/3)
#(g3.ols[1,1]*g3.ols[2,2]*g3.ols[3,3])^(1/3)
# (g3.lad[1,1]*g3.lad[2,2]*g3.lad[3,3])^(1/3)


Y.ols <- step_pol(B.ols,X)
Y.pmm3 <- step_pol(B.pmm3,X)
Y.mle <- step_pol(B.mle,X)
#
#
plot(X,Y, pch = 20, cex = 2)
lines(X,Y.t, type = "b", pch = 1, lwd = 2)
lines(X,Y.ols, type = "b", pch = 15, cex = 1, col = "red", lwd = 2, lty = 2)
lines(X,Y.pmm3, type = "b", pch = 17, cex = 1, col = "green", lwd = 2, lty = 3)
lines(X,Y.mle, type = "b", pch = 18, cex = 1, col = "blue", lwd = 2, , lty = 4)

# Определим положение, названия и цвета:

location = "topleft"
labels = c("Teoret", "OLS", "PMM3", "MLE")
colors = c("black", "red", "green", "blue")

# Если цвет передать в параметр fill, то по умолчанию
# нарисуются цветовые плашки:
# legend(location, labels, fill=colors)

pts = c(1 , 15, 17, 18) # каждый элемент показывается точкой типа 20
lns = c(1, 1, 1, 1) # каждый элемент показывается линией толщиной 1

# теперь посмотрим на легенду (она нарисуется поверх старой)
legend(location, labels, col = colors, pch = pts, lwd = lns)

# x <- c(-1, 1)
# f <- dnormp(x, p=p_epf)
# print(f)
# plot(function(x) dnormp(x, p=p_epf) , -4, 4,
# main = "Exponential power distribution density function", ylab="f(x)")

boxplot(BB.ols[,1:1], BB.pmm3[,1:1], BB.mle[,1:1], col = "lightgray", horizontal = FALSE)
boxplot(BB.ols[,2:2], BB.pmm3[,2:2], BB.mle[,2:2], col = "lightgray", horizontal = FALSE)




