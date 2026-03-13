# функция нахождения решения системы уравнений ММПл при S=2 (для модели экспоненциальной регресии y ~ a*exp(b*X))
solvNL2 <- function(V.ols,X,Y,k2,g3,g4) 
{
  P <- length(V.ols)
  A <- g3; B <- -((2*g3*Y-k2^(1/2)*(g4+2))); C <- ((g3*Y^2-k2^(1/2)*(g4+2)*Y-k2*g3));

  V.pmm2 <- V.ols 
  JZs <- matrix(data = NA, nrow = P, ncol=P)

  a <- V.ols[1];   b <- V.ols[2]; 

  for (i in 1:100)
  {  
      Yx <- a * exp(b*X)
      E <- A * Yx^2 + B * Yx + C
      Z1 <- E * exp(b * X)
      Z2 <- E * a * (exp(b * X) * X)
      
      Z <- c(sum(Z1), sum(Z2))
     
      JZ11 <-exp(b * X) * (A * (a * exp(b * X))^2 + B * (a * exp(b * X)) + 
                             C) + a * exp(b * X) * (A * (2 * (exp(b * X) * (a * exp(b * 
                                                                                      X)))) + B * exp(b * X))
      JZ12 <- a * (exp(b * X) * X) * (A * (a * exp(b * X))^2 + B * (a * exp(b * 
                                                                              X)) + C) + a * exp(b * X) * (A * (2 * (a * (exp(b * X) * 
                                                                                                                            X) * (a * exp(b * X)))) + B * (a * (exp(b * X) * X)))
      JZ21 <- (exp(b * X) * X) * (A * (a * exp(b * X))^2 + B * (a * exp(b * 
                                                                          X)) + C) + (a * (exp(b * X) * X)) * (A * (2 * (exp(b * X) * 
                                                                                                                           (a * exp(b * X)))) + B * exp(b * X))
      JZ22 <- a * (exp(b * X) * X * X) * (A * (a * exp(b * X))^2 + B * (a * 
                                                                          exp(b * X)) + C) + (a * (exp(b * X) * X)) * (A * (2 * (a * 
                                                                                                                                   (exp(b * X) * X) * (a * exp(b * X)))) + B * (a * (exp(b * 
                                                                                                                                                                                           X) * X)))
      
      JZs[1,1] <- sum(JZ11)
      JZs[1,2] <- sum(JZ12)
      JZs[2,1] <- sum(JZ21)
      JZs[2,2] <- sum(JZ22)
      
     
      Q <- solve(JZs, Z)

      V.pmm2 <- V.pmm2 - Q
      a <- V.pmm2[1]; b <- V.pmm2[2];
      
      if (norm(matrix(Q)) < 10^-3) break      
      
  }
    return(V.pmm2)
}