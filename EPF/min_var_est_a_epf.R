min_var_est_a_epf <- function(p_epf,n,w) 
{
  
  A.t <- 1
  
  AA.mean <- A.t
  # AA.median <- A.t
  # AA.midrange <- A.t
  AA.mle <- A.t
  AA.pmm3 <- A.t
  

  for (j in 1:w)
  { 
    
    err <- rnormp(n, mu = 0, sigmap = 1, p = p_epf, method = "def")
    X <- A.t + err
    A.mean <- mean(X)
    res <- X - A.mean
    # A.median <- median(X)
    # A.midrange <- (max(X)+min(X))/2
    MLE <- paramp(X)
    A.mle <- as.numeric(MLE["mp"])
    A.pmm3 <- a_pmm3ex(X)
    AA.mean <- rbind(AA.mean, A.mean)
    # AA.median <- rbind(AA.median, A.median)
    # AA.midrange <- rbind(AA.midrange, A.midrange)
    AA.mle <- rbind(AA.mle, A.mle)
    AA.pmm3 <- rbind(AA.pmm3, A.pmm3)
  }
  
  AA.mean <- AA.mean[1:w+1]
  # AA.median <- AA.median[1:w+1]
  # AA.midrange <- AA.midrange[1:w+1]
  AA.mle <- AA.mle[1:w+1]
  AA.pmm3 <- AA.pmm3[1:w+1]
  
  V.mean <- var(AA.mean)
  # V.median <- var(AA.median)
  # V.midrange <- var(AA.midrange)
  V.mle <- var(AA.mle)
  V.pmm3 <- var(AA.pmm3)
  
  ind_min_var<-which.min(c(V.mean, V.mle, V.pmm3))
  
  return(ind_min_var)
}