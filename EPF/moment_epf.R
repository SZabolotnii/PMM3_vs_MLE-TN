# функция вычисления значения моментов EPD 
moment.epf <- function(k,p,sigma) {
  # mk <- gamma((k+1)/p)*gamma(1/p)^((k-2)/2)*gamma(3/p)^(-k/2)*sigma^k
  mk <- (sigma^k)*(p^(k/p))*gamma((k+1)/p)/(gamma(1/p))
  
  return(mk)
}