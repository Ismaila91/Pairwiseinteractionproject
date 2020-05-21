require("spatstat")
window = owin(xrange = c(0,1), yrange = c(0,1))
Nbre_points <- 100

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
Y <- gridcentres(window = window, 10, 10)
Y <- ppp(Y$x, Y$y, window = window)
#X.strauss <- rStrauss(beta = 100, gamma = 1, R=0.05, W = window)
strauss.param <- list(cif="strauss", par=list(beta=Nbre_points/area(window), gamma=.1, r=0.125), w=window)
strauss.mod <- rmh(strauss.param, start=list(n.start=Nbre_points), control=list(nrep=1e6), verbose=FALSE)
Sk_ux.strauss <- Create.basis(strauss.mod, Y, Rmax = 0.125, rmin = 0.001, 50, xrange=window$xrange, yrange=window$yrange)
theta.est.strauss <- ppm(strauss.mod ~ ., covariates = Sk_ux.strauss)
theta.true.strauss <- ppm(strauss.mod ~ 1, Strauss(r = 0.05))

par(mfrow = c(1,2))
plot(predict(theta.est.strauss, Y, type = "cif")); plot(predict(theta.true.strauss, Y, type ="cif"))
par(mfrow = c(1,1))

Y.res <- gridcentres(window = window, 200, 200)
Y.res <- ppp(Y.res$x,Y.res$y,window = window)
par(mfrow = c(1,2))
plot(predict(theta.est,Y.res, type ="cif")); plot(predict(theta.true,Y.res, type ="cif"))
par(mfrow = c(1,1))

###############################################################################################
############### Fourier-Bessel basis

## Example 2: Soft core point process
softcore.param <- list(cif="sftcr", par=list(beta=Nbre_points/area(window), kappa=0.5, sigma=0.1), w=window)
softcore.mod <- rmh(softcore.param, start=list(n.start=Nbre_points), control=list(nrep=1e6), verbose=FALSE)
Sk_ux.softcore <- Create.basis(softcore.mod, Y, Rmax=0.06, rmin=0.001, 50, xrange=window$xrange, yrange=window$yrange)
theta.est.softcore <- ppm(softcore.mod ~ ., covariates = Sk_ux.softcore)
theta.true.softcore <- ppm(softcore.mod ~ 1, Softcore(kappa=0.2,sigma=0.01), rbord=0.001)
par(mfrow = c(1,2))
plot(predict(theta.est.softcore, Y, type = "cif")); plot(predict(theta.true.softcore, Y, type ="cif"))
par(mfrow = c(1,1))

## Example 3: Diggle-Gates-Stibbard point process
dgs.param <- list(cif="dgs", par=list(beta=Nbre_points/area(window), rho=0.01), w=window)
dgs.mod <- rmh(dgs.param, start=list(n.start=Nbre_points), control=list(nrep=1e6), verbose=FALSE) 
Sk_ux.dgs <- Create.basis(dgs.mod, Y, Rmax = 10, rmin = 0.001, 50, xrange=window$xrange, yrange=window$yrange)
theta.est.dgs <- ppm(dgs.mod ~ ., covariates = Sk_ux.dgs)
theta.true.dgs <- ppm(dgs.mod ~ 1, DiggleGatesStibbard(rho=0.01), rbord=0.15)
par(mfrow = c(1,2))
plot(predict(theta.est.dgs, Y, type = "cif")); plot(predict(theta.true.dgs, Y, type ="cif"))
par(mfrow = c(1,1))

## Example 4: 


####################################################################################################

####------- Pairwise interaction functions

## 1. Hard-core, Strauss, and hybrid Strauss-hard core processes
StrHc.pr <- function(t, hc=1, gamma=0.5, R=2, inter='hardcore') {
  T <- length(t); e.t <- rep(NA, T)
  for(i in 1:T) {
    e.t[i] <- switch(inter,
                     'hardcore' = { ifelse(t[i]>hc, 1, 0) },
                     'strauss' = { ifelse(t[i]<=R, gamma, 1) },
                     'strhc' = { ifelse(t[i]<=hc, 0, ifelse(((hc<t[i])&(t[i]<=R)), gamma, 1)) }
    )}
  return(e.t)
}

curve(StrHc.pr(x), 0, 3) 
curve(StrHc.pr(x, inter='strauss'), 0, 3)
curve(StrHc.pr(x, inter='strhc'), 0, 3)

## 2. Diggle-Gratton potential
DG.int <- function(t, rho, delta, kappa) {
  T <- length(t); e.t <- rep(NA,T)
  for(i in 1:T){
    e.t[i] <- ifelse(t[i]<delta, 0, ifelse(((delta<=t[i])&(t[i]<rho)), ((t[i]-delta)/(rho-delta))^kappa, 
                                           1))
  }
  return(e.t)
}
curve(DG.int(x,rho=3,delta=1,kappa=2), 0, 5)

## 3. Diggle-Gates-Stibbard potential
DGS.int <- function(t, rho=0.01) {
  T <- length(t); e.t <- rep(NA,T)
  for(i in 1:T) {
    if(t[i]<rho) { e.t[i] <- (sin(pi*t[i]/(2*rho)))^2 }
    else { e.t[i] <- 1 }
  }
  return(e.t)
}
curve(DGS.int(x,rho=1), 0, 2)

## 4. Soft core interaction, 0<kappa<1 and sigma>0
Soft.int <- function(r, kappa=0.1, sigma=0.1){
  return( exp( -(sigma/r)^(2/kappa) ) )
}
curve(Soft.int(x,sigma=1), 0, 5)
abline(h=1,lty="dashed")

## 5. Fiksel potential
Fik.int <- function(t, hc, r, kappa, a) {
  T <- length(t); e.t <- rep(NA,T)
  for(i in 1:T) {
    e.t[i] <- ifelse(t[i]<hc, 0, ifelse(((hc<=t[i])&(t[i]<r)), exp(a*exp(-kappa*t[i])), 1))  
  }
  return(e.t)
}
curve(Fik.int(x,hc=1,r=3,kappa=1,a=-2), 0, 5)

####################################################################################################

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

kthcoeff.FourierBessel <- function(r, R, k){ # r>rmin=0.001
  num <- besselJ((r*bzeros[k]/R), nu = 0)
  denom <- sqrt(0.5)*(R*besselJ(bzeros[k], nu = 1))
  return(num/denom)
}

## Orthogonal series estimation of the pairwise interaction function
Phi.hat <- function(r, R, Theta.hat) {
  K <- length(Theta.hat); som <- 0
  for(k in 1:K) {
    som <- Theta.hat[k]*kthcoeff.FourierBessel(r,R=R,k=k) + som
  }
  return(som)
}

curve(StrHc.pr(x, R=0.10, inter='strauss'), 0.001, 0.15)
curve(Phi.hat(x,R=9,Theta.hat=coef(theta.est.strauss)[-1]), from=0.001, to=0.15, lty="dashed", add=TRUE)

curve(Soft.int(x,sigma=.01,kappa=0.2), 0.001, 0.15, main="softcore interaction")
curve(Phi.hat(x,R=12,Theta.hat=coef(theta.est.softcore)[-1]), from=0.001, to=0.15, lty="dashed", add=TRUE)

curve(DGS.int(x,rho=0.01), 0.001, 0.15)
curve(Phi.hat(x,R=2,Theta.hat=coef(theta.est.dgs)[-1]), from=0.001, to=0.15, lty="dashed", add=TRUE)

