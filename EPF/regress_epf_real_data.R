# install.packages("normalp")
library(normalp)
library(L1pack)
# library(parallel)

setwd("D:/Dropbox/_My_Doc/Regression/EPF/")

rm(list = ls())

library(readr)
bluegills <- read_delim("D:/Dropbox/_My_Doc/Regression/EPF/bluegills.csv",
                        ";", escape_double = FALSE, trim_ws = TRUE)
# library(readxl)
# Bmi <- read_excel("D:/Dropbox/_My_Doc/Regression/EPF/Bmi.xlsx",
#                   col_names = FALSE)


source("moment_epf.R")
source("solvMRs3ex.R")
source("step_pol.R")

X <- bluegills[,1]
Y <- bluegills[,2]

# X <- Bmi[,3]
# Y <- Bmi[,4]

reg.ols <- lm(length ~ age + I(age^2), data = bluegills)
# reg.ols <- lm(X3 ~ X2 + I(X2^2), data = Bmi)
# reg.ols <- lm(Y ~ X + I(X^2))

summary(reg.ols)
B.ols <- coefficients(reg.ols)
res <- residuals(reg.ols)

hist(res)
param_res <- paramp(res)

ks.test(res,pnormp, mu = 0, sigmap = param_res$sp, p = param_res$p)

reg.mle <- lmp(length ~ age + I(age^2), data = bluegills)
# reg.mle <- lmp(X3 ~ X2 + I(X2^2), data = Bmi)
summary(reg.mle)
B.mle <- coefficients(reg.mle)

par(mfrow = c(2, 2))
plot.lmp(reg.mle)


m2 <- mean(res^2); m4 <- mean(res^4); m6 <- mean(res^6)
g4 <- m4/m2^2-3; g6 <- m6/m2^3-15*(g4+1)

A=1;
B=-Y*3;
C=3*Y^2-(m6-3*m4*m2)/(m4-3*m2^2);
D=Y*(m6-3*m4*m2)/(m4-3*m2^2)-Y^3;
B.pmm3 <- t(solvMRs3ex(B.ols,cbind(X^0, X, X^2),Y,m2,m4,m6))

g3ef.e <- 1-g4^2/(6+9*g4+g6)


X.t <-min(X):max(X)

Y.ols <- step_pol(B.ols,X.t)
Y.pmm3 <- step_pol(B.pmm3,X.t)
Y.mle <- step_pol(B.mle,X.t)

par(mfrow = c(1, 1))
# plot(t(Bmi[,3]),t(Bmi[,4]))
plot(t(X),t(Y))

lines(X.t,Y.ols, type = "l", col = "red") 
lines(X.t,Y.pmm3, type = "l", col = "green") 
lines(X.t,Y.mle, type = "l", col = "blue") 



regr <- function(data, indices) {
  # вектор indices будет формироваться функцией boot()
  dat <- data[indices, ]
  fit <- lmp(length ~ age + I(age^2), data = dat)
  return(coefficients(fit))
}
# Теперь подадим regr() на функцию boot():
  library(boot)
results <- boot(data = bluegills, statistic = regr, R = 1000)
results
