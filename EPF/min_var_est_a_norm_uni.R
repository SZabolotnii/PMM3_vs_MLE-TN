min_var_est_a_norm_uni <- function(q,n,w) 
{
  
  A.t <- 1
  
  AA.mean <- A.t
  # AA.median <- A.t
  AA.midrange <- A.t
  AA.pmm3 <- A.t
  
  n1 <- round(n * (1-q), digits = 0)
  n2 <- n-n1
  
  for (j in 1:w)
  { 
    err <- c(rnorm(n1, 0, 1), runif(n2, -sqrt(3), sqrt(3)))
    X <- A.t + err
    A.mean <- mean(X)
    res <- X - A.mean
    # A.median <- median(X)
    A.midrange <- (max(X)+min(X))/2
    A.pmm3 <- a_pmm3ex(X)
    AA.mean <- rbind(AA.mean, A.mean)
    # AA.median <- rbind(AA.median, A.median)
    AA.midrange <- rbind(AA.midrange, A.midrange)
    AA.pmm3 <- rbind(AA.pmm3, A.pmm3)
  }
  
  AA.mean <- AA.mean[1:w+1]
  # AA.median <- AA.median[1:w+1]
  AA.midrange <- AA.midrange[1:w+1]
  AA.pmm3 <- AA.pmm3[1:w+1]
  
  V.mean <- var(AA.mean)
  # V.median <- var(AA.median)
  V.midrange <- var(AA.midrange)
  V.pmm3 <- var(AA.pmm3)
  
  ind_min_var<-which.min(c(V.mean, V.midrange, V.pmm3))
  
  return(ind_min_var)
}