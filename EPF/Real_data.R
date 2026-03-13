library(normalp)
library(sn)
library (moments)

rm(list = ls())

getwd()
# setwd("/home/docsa/Dropbox/_My_Doc/Regression/EPF/")
setwd("D:/Dropbox/_My_Doc/Regression/EPF/")

source("moment_epf.R")
source("solvMRs3ex.R")
source("step_pol.R")

# data <- read.csv(file = 'Movie.csv', header = FALSE, sep = ";")

Y <- RealData$V5
X1 <- RealData$V2
X2 <- RealData$V3
X3 <- RealData$V4


reg.ols <- lm(Y ~ X1 + X2)

B.ols <- coefficients(reg.ols)

summary(reg.ols)

res_OLS <- residuals(reg.ols) 


hist(res_OLS)

skewness(res_OLS, na.rm = FALSE)
kurtosis(res_OLS, na.rm = FALSE)

agostino.test(res_OLS, alternative = c("two.sided", "less", "greater"))
jarque.test(res_OLS)

param_res <- paramp(res_OLS)

ks.test(res_OLS, pnormp, mu = 0, sigmap = param_res$sp, p = param_res$p)

reg.mle <- lmp(Y ~ X1 + X2 + X3)

B.mle <- coefficients(reg.mle)


summary(reg.mle)

paramp(res_OLS)

estimatep(res_OLS, 1, 2)





m2 <- mean(res_OLS^2); m4 <- mean(res_OLS^4); m6 <- mean(res_OLS^6)
g4 <- m4/m2^2-3; g6 <- m6/m2^3-15*(g4+1)

A=1;
B=-Y*3;
C=3*Y^2-(m6-3*m4*m2)/(m4-3*m2^2);
D=Y*(m6-3*m4*m2)/(m4-3*m2^2)-Y^3;
B.pmm3 <- t(solvMRs3ex(B.ols,cbind(X^0, X),Y,m2,m4,m6))

g3ef.e <- 1-g4^2/(6+9*g4+g6)

X.t <-min(X):max(X)

Y.ols <- step_pol(B.ols,X.t)
Y.pmm3 <- step_pol(B.pmm3,X.t)
Y.mle <- step_pol(B.mle,X.t)

par(mfrow = c(1, 1))
plot(t(X),t(Y))

lines(X.t,Y.ols, type = "l", col = "red") 
lines(X.t,Y.pmm3, type = "l", col = "green") 
lines(X.t,Y.mle, type = "l", col = "blue") 



par(mfrow = c(2, 2))
plot.lmp(reg.mle)

plot(Y ~ X)
abline(lm(Y ~ X))
abline(lmp(Y ~ X), col = "red")

