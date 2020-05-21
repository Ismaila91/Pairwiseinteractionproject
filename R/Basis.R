Cosine.Basis <- function(r, R, k) {
  res <- rep(NA,length(r))
  cond <- (k==1)
  for(i in 1:length(r)){
    res[i] <- ifelse(cond, 1/sqrt(R), sqrt(2/R)*cos((k-1)*pi*r[i]/R))
  }
  return(res)
}
curve(Cosine.Basis(x,R=0.125,k=1),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=2),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=3),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=4),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=5),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=6),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=7),0,0.12,xlab="r",ylab="phi")
curve(Cosine.Basis(x,R=0.125,k=8),0,0.12,xlab="r",ylab="phi")

N <- 100
bzeros <- CircularDDM::besselzero(0,N,1)

FourierBessel.Basis <- function(r, R, k){
  num <- besselJ((r*bzeros[k]/R), nu = 0)
  denom <- sqrt(0.5)*(R*besselJ(bzeros[k], nu = 1))
  return(num/denom)
}
curve(FourierBessel.Basis(x,R=0.125,k=1),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=2),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=3),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=4),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=5),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=6),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=7),0,0.12,xlab="r",ylab="phi")
curve(FourierBessel.Basis(x,R=0.125,k=8),0,0.12,xlab="r",ylab="phi")




