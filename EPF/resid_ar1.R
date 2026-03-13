resid_ar1 <- function(y,a0,a1) 
{
  y0 <- c(y,0)
  y1 <- c(0,y)
  res <- y0 - (a0 + a1 * y1)
  return(res[1:length(y)])
}