# Встановлення бібліотек
library(MASS)
library(normalp)


rm(list = ls())


# Параметри моделі
n <- 500
w <- 1000
p_epf <- 50
B.t <- c(1, -2) # b0, b1, ...
L <- length(B.t)
set.seed(123) # для відтворюваності

# Попереднє розміщення пам'яті
BB.ols <- matrix(NA, nrow = w, ncol = L)
BB.mle <- matrix(NA, nrow = w, ncol = L)
BB.robust <- matrix(NA, nrow = w, ncol = L)
MSE.ols <- numeric(w)
MSE.mle <- numeric(w)
MSE.robust <- numeric(w)
Pred.diff.ols <- numeric(w)
Pred.diff.mle <- numeric(w)
Pred.diff.robust <- numeric(w)

# Функція для генерації U-квадратичних випадкових величин
runifq <- function(n) {
  u <- runif(n, min=-sqrt(3), max=sqrt(3))
  return(u - (u^3)/sqrt(3))
}

X <- matrix(runif(n * L), nrow = n, ncol = L)
X[,1] <- 1

for (j in 1:w) {
  indices <- sample(1:n, size = 0.8 * n)
  X_train <- X[indices,]
  X_test <- X[-indices,]
  
  # err <- runifq(n)
  err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
  
  Y <- as.vector(X %*% B.t) + err
  Y_train <- Y[indices]
  Y_test <- Y[-indices]
  
  # OLS оцінювання
  reg.ols <- lm(Y_train ~ X_train[, -1, drop = FALSE])
  BB.ols[j,] <- coefficients(reg.ols)
  Y_pred_ols <- predict(reg.ols, newdata = data.frame(X_test[, -1, drop = FALSE]))
  MSE.ols[j] <- mean((Y_pred_ols - Y_test)^2)
  Pred.diff.ols[j] <- mean(abs(Y_pred_ols - Y_test))
  
  # MLE оцінювання
  reg.mle <- lmp(Y_train ~ X_train[, -1, drop = FALSE])
  BB.mle[j,] <- coefficients(reg.mle)
  Y_pred_mle <- predict(reg.mle, newdata = data.frame(X_test[, -1, drop = FALSE]))
  MSE.mle[j] <- mean((Y_pred_mle - Y_test)^2)
  Pred.diff.mle[j] <- mean(abs(Y_pred_mle - Y_test))
  
  # Робастна регресія
  reg.robust <- rlm(Y_train ~ X_train[, -1, drop = FALSE])
  BB.robust[j,] <- coefficients(reg.robust)
  Y_pred_robust <- predict(reg.robust, newdata = data.frame(X_test[, -1, drop = FALSE]))
  MSE.robust[j] <- mean((Y_pred_robust - Y_test)^2)
  Pred.diff.robust[j] <- mean(abs(Y_pred_robust - Y_test))
}


# Середні MSE та середні різниці
mean_MSE_ols <- mean(MSE.ols)
mean_MSE_robust <- mean(MSE.robust)
mean_MSE_mle <- mean(MSE.mle)
mean_Pred_diff_ols <- mean(Pred.diff.ols)
mean_Pred_diff_robust <- mean(Pred.diff.robust)
mean_Pred_diff_mle <- mean(Pred.diff.mle)

# Результати
list(
  mean_MSE_ols = mean_MSE_ols, 
  mean_MSE_robust = mean_MSE_robust,
  mean_MSE_mle = mean_MSE_mle,
  mean_Pred_diff_ols = mean_Pred_diff_ols, 
  mean_Pred_diff_robust = mean_Pred_diff_robust,
  mean_Pred_diff_mle = mean_Pred_diff_mle
)

MSE.ols <- sapply(1:L, function(i) mean((BB.ols[-1, i] - B.t[i])^2))
MSE.mle <- sapply(1:L, function(i) mean((BB.mle[-1, i] - B.t[i])^2))
MSE.robust <- sapply(1:L, function(i) mean((BB.robust[-1, i] - B.t[i])^2))

g3.mle <- MSE.mle / MSE.ols
g3.robust <- MSE.robust / MSE.ols

hist(err)

