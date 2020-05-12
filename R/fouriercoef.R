## Soft core interaction, 0<kappa<1 and sigma>0
Soft.int <- function(r, kappa=0.1, sigma=0.1){
  return( exp( -(sigma/r)^(2/kappa) ) )
}
curve(Soft.int(x), 0, 5)
abline(h=1,lty="dashed")

## kth element of the Cosine basis
kthcoeff.Cosine <- function(r, R, k) {
  if(k==1) {
    res <- rep(1/sqrt(R),length(r))
  }
  else {
    res <- sqrt(2/R)*cos(pi*(k-1)*r/R)
  }
  return(res)
}

## kth element of the Fourier-Bessel basis
N <- 100
bzeros <- CircularDDM::besselzero(0,N,1)

kthcoeff.FourierBessel <- function(r, R, rmin=0.001, k){
  num <- besselJ((r*bzeros[k]/R), nu = 0)
  denom <- sqrt(0.5)*(R*besselJ(bzeros[k], nu = 1))
  return(num/denom)
}

## Weight function
Wght <- function(r, d=2) {
  return(r^(d-1))
}

IIntegrand <- function(r, R=0.8, rmin=0.001, k, kappa=0.1, sigma=0.1, basis='cosine'){
  intg <- switch(basis,
                 'cosine' = Soft.int(r+rmin)*kthcoeff.Cosine(r,R=R,k=k),
                 'fourierbessel' = Soft.int(r+rmin)*kthcoeff.FourierBessel(r,R=R,k=k)*Wght(r))
  return(intg)
}

FourierCoef <- function(K, R=1, rmin=0.001,kappa=0.1,sigma=0.1,basis='cosine') {
  Fcoef <- rep(NA, K)
  for(k in 1:K){
    Fcoef[k] <- switch(basis,
                       'cosine' = integrate(IIntegrand,lower=0,upper=R,k=k,R=R)$value,
                       'fourierbessel' = integrate(IIntegrand,lower=0,upper=R,R=R,rmin=rmin,k=k,
                                                   kappa=kappa,sigma=sigma,basis=basis)$value)
  }
  return(Fcoef)
}

Phi.u_Approx <- function(r, K, R=1, rmin=0.001, kappa=0.1, sigma=0.1, basis='cosine') {
  Fcoef <- switch(basis,
                  'cosine' = FourierCoef(K,R=R),
                  'fourierbessel' = FourierCoef(K,R=R,rmin=rmin,kappa=kappa,sigma=sigma,basis=basis))
  ssum <- 0
  for(k in 1:K){
    ssum <- switch(basis,
                   'cosine' = { Fcoef[k]*kthcoeff.Cosine(r,R=R,k=k) + ssum },
                   'fourierbessel' = { Fcoef[k]*kthcoeff.FourierBessel(r,R=R,rmin=rmin,k=k) + ssum } )
  }
  return(ssum)
}

curve(Soft.int(x,sigma=.1), 0, 5, main="softcore interaction")
curve(Phi.u_Approx(x,K=10,R=5), from=0, to=5, lty="dashed", add=TRUE)
curve(Phi.u_Approx(x,K=10,sigma=1,R=5,basis='fourierbessel'), from=0, to=5, lty="dashed", add=TRUE)
curve(Phi.u_Approx(x,K=20,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=30,sigma=1,R=5), from=0, to=5, lty="twodash", add=TRUE)
curve(Phi.u_Approx(x,K=40,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=50,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=100,sigma=1,R=5,basis='fourierbessel'), from=0, to=5, lty="longdash", add=TRUE)
