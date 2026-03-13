min_var_est_a_epf_apply <- function(p_epf,n,w) 
{
  
  A.t <- 1
  
 
  err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
  for (j in 1:(w-1))
  {
    err.j <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
    err <- rbind(err, err.j)
  }
  
  
  X <- A.t + err
  
  AA.mean <-apply(X, MARGIN = 1, FUN = mean)
  
  AA.mle <- apply(X, 1, # Вызов нашей функции для каждой строчки
                  FUN = function(X)  # Определяем нашу функцию прямо в вызове apply
                  {
                    MLE <- paramp(X)
                    A.mle <- as.numeric(MLE["mp"])
                    (A.mle)
                  }
  )
  
  AA.pmm3 <-apply(X, MARGIN = 1, FUN = a_pmm3ex)
  
  
  
  V.mean <- var(AA.mean)
  # V.median <- var(AA.median)
  # V.midrange <- var(AA.midrange)
  V.mle <- var(AA.mle)
  V.pmm3 <- var(AA.pmm3)
  
  ind_min_var<-which.min(c(V.mean, V.mle, V.pmm3))
  
  return(ind_min_var)
}