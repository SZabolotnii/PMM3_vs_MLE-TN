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

# Приклад виклику функції
set.seed(123) # Для відтворюваності результатів
random_vars <- generate_random_variables(1000, 0, 1) # Генерація 1000 випадкових величин з a = 0 і b = 1
hist(random_vars, breaks=50, main='Histogram of V-shaped distribution')


# Функція для обчислення 2k-го центрального моменту
central_moment_2k <- function(a, b, k) {
  # Інтегрування здійснюється від a до a+b і результат множиться на 2
  moment <- 2 * integrate(function(x) ((x - a)^(2*k) * (x - a) / b^2), lower = a, upper = a + b)$value
  return(moment)
}

# Обчислення 2-го і 4-го центральних моментів
moment_2 <- central_moment_2k(a = 0, b = 1, k = 1)
moment_4 <- central_moment_2k(a = 0, b = 1, k = 2)
# Обчислення 6-го центрального моменту
moment_6 <- central_moment_2k(a = 0, b = 1, k = 3)

# Виведення результатів
print(paste("2nd central moment:", moment_2))
print(paste("4th central moment:", moment_4))
print(paste("6th central moment:", moment_6))
