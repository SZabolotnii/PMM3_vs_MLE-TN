# функция нахождения решения системы уравнений ММПл при S=2 (для модели экспоненциальной регресии y ~ a*exp(b*X))
solvNL2_nleqslv <- function(V.ols,X,Y,k2,g3,g4) 
{
  A <- g3; B <- -((2*g3*Y-k2^(1/2)*(g4+2))); C <- ((g3*Y^2-k2^(1/2)*(g4+2)*Y-k2*g3));

  syst_expregress_S2 <- function(x) {
    y <- numeric(2)
    y[1] <- sum((A * (x[1] * exp(x[2]*X))^2 + B * (x[1] * exp(x[2]*X)) + C) * exp(x[2] * X))
    y[2] <- sum((A * (x[1] * exp(x[2]*X))^2 + B * (x[1] * exp(x[2]*X)) + C) * x[1] * (exp(x[2] * X) * X))
    y
  }
  V.pmm2 <- nleqslv(V.ols, syst_expregress_S2)
  
  if (abs(V.pmm2$x[1]) > 2*abs(V.ols[1]) | abs(V.pmm2$x[2]) > 2*abs(V.ols[2])) V.pmm2$x <- V.ols
  
  return(V.pmm2$x)
}