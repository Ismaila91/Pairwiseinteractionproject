Cosine.Coef <- function(r, K, R) {
  Coeff <- rep(NA, K+1)
  Coeff[1] <- 1/sqrt(R)
  for(k in 1:K) {
    Coeff[k+1] <- sqrt(2/R)*cos(k*pi*r/R)
  }
  return(Coeff)
}

A <- Cosine.Coef(r=3,K=10,R=5)
r=as.matrix(c(1,2,3))
apply(r, 1, Cosine.Coef, K=5, R=5)



