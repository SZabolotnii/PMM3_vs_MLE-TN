# функция нахождения решения системы уравнений ММПл при S=3 (симетричная модель ошибки)
solvMRs3ex <- function(B.ols,X,Y,m2,m4,m6) 
{
  P <- length(B.ols)
  A <- 1; B <- -Y*3; C <- 3*Y^2-(m6-3*m4*m2)/(m4-3*m2^2); D <- Y*(m6-3*m4*m2)/(m4-3*m2^2)-Y^3;
  B.pmm3 <- B.ols 
  JZs <- matrix(data = NA, nrow = P, ncol=P)

  for (i in 1:100)
  {  
      Yx <- B.pmm3[1]*X[,1]  
      for (r in 2:P)
         {Yx <- Yx+B.pmm3[r]*X[,r]}
      Z1 <- A*Yx^3+B*Yx^2+C*Yx+D
      Z <- sum(Z1)
      for (r in 2:P)
         {
        Zr <- Z1*X[,r]
        Z <- rbind(Z, sum(Zr))
         }
      JZ11 <- 3*A*Yx^2+2*B*Yx+C
      JZs[1,1] <- sum(JZ11)
      
      for (ii in 2:P)
      {
        JZ1i <- JZ11*X[,ii]
        JZs[1,ii] <- sum(JZ1i);
        JZs[ii,1] <- JZs[1,ii];
          
        for (jj in 2:P)
            {
             JZij=JZ11*X[,ii]*X[,jj]
             JZs[ii,jj]=sum(JZij)
            }
      }
     
      Q <- solve(JZs,Z)
      B.pmm3 <- B.pmm3 - Q
      
      if (norm(matrix(Q)) < 10^-4) break      
      
  }
  
  if (norm((B.ols - B.pmm3), "M" ) > 10) B.pmm3 <- B.ols
  return(B.pmm3)
}