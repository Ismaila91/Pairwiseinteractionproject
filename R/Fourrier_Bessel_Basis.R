require("spatstat")
window = owin(xrange = c(0,2), yrange = c(0,2))
Nbre_points <- 300

X <- rpoispp(lambda = 5, win = window)
Y <- gridcentres(window = window, 2, 2)
Y <- ppp(Y$x,Y$y,window = window)

plot(X)
points(Y,col="red",pch=16)


#ddd$i #index X
#ddd$j #index Y
#ddd$d #distance

Create.basis <- function(X, Y, Rmax = 1 , rmin = 0.001, n, xrange = c(0,1), yrange = c(0,1)){
  ddd <- crosspairs(X, Y, Rmax+rmin)
  bzeros <- CircularDDM::besselzero(0,n,1) #change 1
  eval.basis <- function(i,n0){
    flag <- ddd$j == i
    ratio <- 0
    if(sum(flag) != 0){
      use.dist <- ddd$d[flag]
      use.dist <- use.dist[use.dist>= rmin]  #change 2
      num <- besselJ(((use.dist-rmin)*bzeros[n0]/Rmax), nu = 0)
      denom <- sqrt(0.5)*(Rmax*besselJ(bzeros[n0], nu = 1))
      ratio <- num/denom
    }
    return(cbind(i,sum(ratio))) 
  }
  Y.indexes <- matrix( 1:Y$n , nrow = Y$n )
  basis.img <- list()
  names <- NULL
  for(i in 1:n){
    Phi.u <- apply(Y.indexes,1,eval.basis , n0 = i)
    basis.img[[i]] <- im(matrix(Phi.u[2,],ncol = sqrt(Y$n), nrow = sqrt(Y$n)), xrange = xrange, yrange = yrange)
    names[i] <- paste("basis",i,collapse = "",sep="")
  }
  names(basis.img) <- names
  return(basis.img)
}

my.basis <- Create.basis(X,Y,Rmax = 0.1, rmin = 0.0001, 10)
plot(my.basis$basis10)

## Example 1: Strauss model
Y <- gridcentres(window = window, 20, 20)
Y <- ppp(Y$x, Y$y, window = window)
#X.strauss <- rStrauss(beta = 100, gamma = 1, R=0.05, W = window)
strauss.param <- list(cif="strauss", par=list(beta=Nbre_points/area(window), gamma=.1, r=0.05), w=window)
strauss.mod <- rmh(strauss.param, start=list(n.start=Nbre_points), control=list(nrep=1e6), verbose=FALSE)
Sk_ux <- Create.basis(strauss.mod, Y, Rmax = 0.10, rmin = 0.0001, 10, xrange = c(0,2), yrange = c(0,2))
theta.est <- ppm(strauss.mod ~. + 1, covariates = Sk_ux)
theta.true <- ppm(strauss.mod ~ 1, Strauss(r = 0.15))

par(mfrow = c(1,2))
plot(predict(theta.est, Y, type = "cif")); plot(predict(theta.true, Y, type ="cif"))
par(mfrow = c(1,1))

Y.res <- gridcentres(window = window, 200, 200)
Y.res <- ppp(Y.res$x,Y.res$y,window = window)
par(mfrow = c(1,2))
plot(predict(theta.est,Y.res, type ="cif")); plot(predict(theta.true,Y.res, type ="cif"))
par(mfrow = c(1,1))

###############################################################################################
############### Fourier-Bessel basis
## Example: Soft core point process
softcore.param <- list(cif="sftcr", par=list(beta=Nbre_points/area(window), kappa=0.1, sigma=0.1), w=window)
softcore.mod <- rmh(softcore.param, start=list(n.start=Nbre_points), control=list(nrep=1e6), verbose=FALSE)
Sk_ux <- Create.basis(softcore.mod, Y, Rmax = 0.10, rmin = 0.0001, 20, xrange = c(0,2), yrange = c(0,2))
theta.est <- ppm(softcore.mod ~ ., covariates = Sk_ux)
theta.true <- ppm(softcore.mod ~ 1, Strauss(r = 0.15))
par(mfrow = c(1,2))
plot(predict(theta.est, Y, type = "cif")); plot(predict(theta.true, Y, type ="cif"))
par(mfrow = c(1,1))




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

Phi.hat <- function(r, R, rmin=0.001) {
  Theta.hat <- coef(theta.est)[-1]
  K <- length(Theta.hat); som <- 0
  for(k in 1:K) {
    som <- Theta.hat[k]*kthcoeff.FourierBessel(r,R=R,rmin=rmin,k=k) + som
  }
  return(som)
}


curve(Soft.int(x,sigma=.1), 0, 5, main="softcore interaction")
curve(Phi.hat(x,R=5), from=0, to=5, lty="dashed", add=TRUE)
curve(Phi.u_Approx(x,K=10,R=5), from=0, to=5, lty="dashed", add=TRUE)
curve(Phi.u_Approx(x,K=10,sigma=1,R=5,basis='fourierbessel'), from=0, to=5, lty="dashed", add=TRUE)
curve(Phi.u_Approx(x,K=20,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=30,sigma=1,R=5), from=0, to=5, lty="twodash", add=TRUE)
curve(Phi.u_Approx(x,K=40,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=50,sigma=1,R=5), from=0, to=5, lty="longdash", add=TRUE)
curve(Phi.u_Approx(x,K=100,sigma=1,R=5,basis='fourierbessel'), from=0, to=5, lty="longdash", add=TRUE)


