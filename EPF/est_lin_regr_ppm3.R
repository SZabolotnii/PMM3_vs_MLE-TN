library(normalp)
library(L1pack) # Закоментовано, оскільки не використовується в прикладі

library(quantreg)
library(MASS)
library(robustbase)
library(gamlss)

setwd("D:/Dropbox/_My_Doc/R/EPF/")

rm(list = ls())

source("moment_epf.R")
source("solvMRs3ex.R")
source("step_pol.R")

# Встановлення параметрів моделі
n <- 500
w <- 2000
B.t <- c(1, 3) # b0, b1, b2, ...
L <- length(B.t)
p_epf <- 6

# Попереднє розміщення пам'яті для зберігання результатів
BB.ols <- matrix(NA, nrow = w+1, ncol = L)
BB.mle <- matrix(NA, nrow = w+1, ncol = L)
BB.lad <- matrix(NA, nrow = w+1, ncol = L)
BB.pmm3 <- matrix(NA, nrow = w+1, ncol = L)

BB.qr <- matrix(NA, nrow = w+1, ncol = L) # Для квантильної регресії
BB.rlm <- matrix(NA, nrow = w+1, ncol = L) # Для M-оцінок
BB.lmrob <- matrix(NA, nrow = w+1, ncol = L) # Для регресії на основі відсічених середніх
BB.glm <- matrix(NA, nrow = w+1, ncol = L) # Попереднє розміщення пам'яті для зберігання результатів GLM оцінок



# Заповнення першого рядка істинними значеннями B.t
BB.ols[1,] <- B.t
BB.mle[1,] <- B.t
BB.lad[1,] <- B.t
BB.pmm3[1,] <- B.t

BB.qr[1,] <- B.t
BB.rlm[1,] <- B.t
BB.lmrob[1,] <- B.t
BB.glm[1,] <- B.t

# Попереднє обчислення моментів
m2.t <- moment.epf(2, p_epf, 1)
m4.t <- moment.epf(4, p_epf, 1)
m6.t <- moment.epf(6, p_epf, 1)

X <- matrix(NA, nrow = n, ncol = L) # Матриця дизайну
X[,1] <- 1 # Перший стовпець для вільного члена


# Функція для генерації U-квадратичних випадкових величин
runifq <- function(n) {
  u <- runif(n, min=-sqrt(3), max=sqrt(3))  # U(-sqrt(3), sqrt(3)) має дисперсію 1
  return(u - (u^3)/sqrt(3))  # U-квадратичний розподіл
}

# Функція для генерації випадкових чисел з арксинусного розподілу
generate_arcsine_samples <- function(n) {
  # Визначаємо теоретичну дисперсію арксинусного розподілу
  theoretical_variance <- 1/2 - (1/pi)^2
  # Генеруємо рівномірні випадкові числа з інтервалу [0, 1]
  uniform_samples <- runif(n, 0, 1)
  # Використовуємо обернену функцію розподілу для отримання випадкових чисел з арксинусного розподілу
  arcsine_samples <- sin(pi * (uniform_samples - 0.5)) / sqrt(theoretical_variance)
  return(arcsine_samples)
}

# Функція для генерації випадкових величин з симетричного V-подібного розподілу
generate_V_shapered <- function(n, a, b) {
  # Генеруємо n рівномірно розподілених випадкових величин
  u <- runif(n)
  
  # Використовуємо обернену CDF для перетворення u в x
  x <- ifelse(u <= 0.5,
              a + sqrt(2) * b * sqrt(u),         # Для значень u від 0 до 0.5
              a - sqrt(2) * b * sqrt(1 - u))     # Для значень u від 0.5 до 1
  return(x)
}

# Функція для генерації випадкової величини з нормалізованого розподілу Kumaraswamy
rkumar_norm <- function(n, a, b) {
  u <- runif(n)  # Генеруємо n рівномірно розподілених випадкових величин
  x <- (1 - (1 - u)^(1/b))^(1/a)  # Обернена функція розподілу Kumaraswamy
  y <- x - 0.5  # Зсув розподілу
  
  # Обчислення стандартного відхилення розподілу Kumaraswamy
  mean_x <- beta(a + 1, b) / beta(a, b)
  var_x <- beta(a + 1, b) * (beta(a + 1, b) + beta(a, b + 1) - beta(a + 1, b)^2) / beta(a, b)^2
  sigma <- sqrt(var_x)
  
  # Нормалізація дисперсії
  z <- y / sigma
  
  return(z)
}


# Моделювання та оцінювання параметрів
for (j in 1:w) {
  for (k in 2:L) {
    X[,k] <- runif(n) # Заповнюємо матрицю дизайну випадковими числами
  }
  
  # err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
  # err <- runifq(n)
  # err <- rkumar_norm(n, 0.5, 0.5)
  # err <- generate_arcsine_samples(n)
  err <- generate_V_shapered(n,0,1)
  
  # err <- err / sd(err)  # Нормування помилки так, щоб її дисперсія була 1
  
  Y <- as.vector(X %*% B.t) + err # Вектор відгуку
  
  # OLS оцінювання
  reg.ols <- lm(Y ~ X[, -1]) # Виключаємо перший стовпець (вільний член)
  BB.ols[j+1,] <- coefficients(reg.ols)
  
  # MLE оцінювання
  reg.mle <- lmp(Y ~ X[, -1]) # Функція lmp має бути визначена
  BB.mle[j+1,] <- coefficients(reg.mle)
  
  # LAD оцінювання
  reg.lad <- lad(Y ~ X[, -1]) # Функція lad має бути визначена
  BB.lad[j+1,] <- coefficients(reg.lad)
  
  # PMM3 оцінювання
  res <- residuals(reg.ols)
  m2 <- mean(res^2)
  m4 <- mean(res^4)
  m6 <- mean(res^6)
  B.pmm3 <- t(solvMRs3ex(coefficients(reg.ols), X, Y, m2, m4, m6))
  BB.pmm3[j+1,] <- B.pmm3
  
  # Квантильна регресія
  qr_model <- rq(Y ~ X[, -1], tau = 0.5)
  BB.qr[j+1,] <- coefficients(qr_model)
  
  # M-оцінки
  rlm_model <- rlm(Y ~ X[, -1])
  BB.rlm[j+1,] <- coefficients(rlm_model)
  
  # Регресія на основі відсічених середніх
  lmrob_model <- lmrob(Y ~ X[, -1])
  BB.lmrob[j+1,] <- coefficients(lmrob_model)
  
  # GLM оцінювання
  # Припустимо, що ми використовуємо сімейство Гаусса з ідентифікаційною лінк-функцією
  # Це схоже на OLS, але з можливістю зміни розподілу помилок
  reg.glm <- glm(Y ~ X[, -1], family = gaussian(link = "identity"))
  # reg.glm <- glm(Y ~ X[, -1], family = binomial(link = "logit"))
  BB.glm[j+1,] <- coefficients(reg.glm)
}

BB.ols <- BB.ols[1:w+1,1:L]
BB.pmm3 <- BB.pmm3[1:w+1,1:L]
BB.mle <- BB.mle[1:w+1,1:L]
BB.lad <- BB.lad[1:w+1,1:L]

BB.qr<- BB.qr[1:w+1,1:L]
BB.rlm<- BB.rlm[1:w+1,1:L]
BB.lmrob<- BB.lmrob[1:w+1,1:L]
BB.glm<- BB.glm[1:w+1,1:L]

# Обчислення дисперсії, MSE та середніх значень
# V.ols <- apply(BB.ols, 2, var)
# V.mle <- apply(BB.mle, 2, var)
# V.lad <- apply(BB.lad, 2, var)
# V.pmm3 <- apply(BB.pmm3, 2, var)


# Розрахунок MSE для PMM3 оцінок
# Розрахунок MSE для кожного параметра PMM3 оцінок

MSE.ols <- sapply(1:L, function(i) mean((BB.ols[-1, i] - B.t[i])^2))
MSE.mle <- sapply(1:L, function(i) mean((BB.mle[-1, i] - B.t[i])^2))
MSE.lad <- sapply(1:L, function(i) mean((BB.lad[-1, i] - B.t[i])^2))
MSE.pmm3 <- sapply(1:L, function(i) mean((BB.pmm3[-1, i] - B.t[i])^2))
MSE.qr <- sapply(1:L, function(i) mean((BB.qr[-1, i] - B.t[i])^2))
MSE.rlm <- sapply(1:L, function(i) mean((BB.rlm[-1, i] - B.t[i])^2))
MSE.lmrob <- sapply(1:L, function(i) mean((BB.lmrob[-1, i] - B.t[i])^2))
MSE.glm <- sapply(1:L, function(i) mean((BB.glm[-1, i] - B.t[i])^2))


# M.ols <- colMeans(BB.ols)
# M.mle <- colMeans(BB.mle)
# M.lad <- colMeans(BB.lad)
# M.pmm3 <- colMeans(BB.pmm3)
# M.qr <- colMeans(BB.qr[-1,], na.rm = TRUE)
# M.rlm <- colMeans(BB.rlm[-1,], na.rm = TRUE)
# M.lmrob <- colMeans(BB.lmrob[-1,], na.rm = TRUE)
# M.glm <- colMeans(BB.glm[-1,], na.rm = TRUE)


# Обчислення та порівняння ефективності
g4.t <- m4.t / m2.t^2 - 3
g6.t <- m6.t / m2.t^3 - 15 * (g4.t + 1)
g3ef.t <- 1 - g4.t^2 / (6 + 9 * g4.t + g6.t)

g3.pmm3 <- MSE.pmm3 / MSE.ols
g3.mle <- MSE.mle / MSE.ols
g3.lad <- MSE.lad / MSE.ols
g3.qr <- MSE.qr / MSE.ols
g3.rlm <- MSE.rlm / MSE.ols
g3.lmrob <- MSE.lmrob / MSE.ols
g3.glm <- MSE.glm / MSE.ols

efficiency_mle <- prod(g3.mle)^(1/L)
efficiency_lad <- prod(g3.lad)^(1/L)
efficiency_pmm3 <- prod(g3.pmm3)^(1/L)
efficiency_qr <- prod(g3.qr)^(1/L)
efficiency_rlm <- prod(g3.rlm)^(1/L)
efficiency_lmrob <- prod(g3.lmrob)^(1/L)
efficiency_glm <- prod(g3.glm)^(1/L)


efficiency_mle
efficiency_lad
efficiency_pmm3
efficiency_qr
efficiency_rlm
efficiency_lmrob
efficiency_glm

g3ef.t

# Припускаємо, що BB.ols, BB.mle, BB.pmm3 мають однаковий формат і розміри

# Конвертуємо матриці в список датафреймів для кожного параметра
# Кожен датафрейм містить колонки для OLS, MLE, PMM3 оцінок
# Зауважте, що рядки, які містять NA, будуть виключені
parameters_estimates <- lapply(1:L, function(i) {
  data.frame(
    OLS = BB.ols[-1, i],
    MLE = BB.mle[-1, i],
    LAD = BB.lad[-1, i],
    PMM3 = BB.pmm3[-1, i],
    QR = BB.qr[-1, i],
    RLM = BB.rlm[-1, i],
    LMROB = BB.lmrob[-1, i],
    GLM = BB.glm[-1, i]
  )
})

# Боксплот для кожного параметра, включаючи нові методи
par(mfrow = c(1, L)) # Змінюємо розміщення боксплотів для включення нових методів
for (i in 1:L) {
  boxplot(parameters_estimates[[i]], main = paste("Parameter", i), las = 2, names = c("OLS", "MLE", "LAD", "PMM3", "QR", "RLM", "LMROB", "GLM"))
}



