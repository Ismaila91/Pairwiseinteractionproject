require("spatstat")


## kth element of the Fourier-Bessel basis
N <- 100
bzeros <- CircularDDM::besselzero(0,N,1)

kthcoeff.FourierBessel <- function(r, R, k){ # r>rmin=0.001; k <= N=100
  num <- besselJ((r*bzeros[k]/R), nu = 0)
  denom <- sqrt(0.5)*(R*besselJ(bzeros[k], nu = 1))
  return(num/denom)
}


## Covariates created with the Fourier-Bessel basis
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

## Orthogonal series estimation of the pairwise interaction function
Phi.hat <- function(r, R, Theta.hat) {
  K <- length(Theta.hat)
  return(sum(Theta.hat*kthcoeff.FourierBessel(r,R=R,k=1:K))) 
}

## Monte carlo local estimation of Phi.hat
FSim <- function(Model="strauss",Y,r,R,nsim=1,kap=NULL,sig=NULL,gam=NULL,window=square(1),npts,K,f.dummy=2,int.rg=NULL,
                 Rho=NULL,Delta=NULL,Hc=NULL,a=NULL) {
  mod.par <- switch(Model,
                    'strauss'= { list(cif=Model,par=list(beta=npts/area(window),gamma=gam,r=int.rg),w=window) },
                    'sftcr'= { list(cif=Model,par=list(beta=npts/area(window),kappa=kap,sigma=sig),w=window) },
                    'dgs'= { list(cif=Model,par=list(beta=npts/area(window),rho=Rho),w=window) },
                    'diggra'= { list(cif=Model,par=list(beta=npts/area(window),kappa=kap,delta=Delta,rho=Rho),w=window) },
                    'fiksel'= { list(cif=Model,par=list(beta=npts/area(window),r=int.rg,hc=Hc,kappa=kap,a=a),w=window) })
  PhiEst <- rep(NA, nsim)
  for(i in 1:nsim){
    mod <- rmh(mod.par, start=list(n.start=npts), control=list(nrep=1e6), verbose=FALSE)
    Sk_ux <- Create.basis(mod,Y=Y,Rmax=R,rmin=0.001,K,xrange=window$xrange,yrange=window$yrange)
    n <- npoints(mod); n.dummy <- round(f.dummy*sqrt(n))
    ThetaEst <- ppm(mod ~ ., covariates=Sk_ux, nd=n.dummy)
    PhiEst[i] <- Phi.hat(r=r, R=R, Theta.hat=coef(ThetaEst)[-1])
    cat(i, " ")
    gc()
  }
  return(PhiEst)
}




