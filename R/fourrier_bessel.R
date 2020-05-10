require("spatstat")
require("CircularDDM")
N0 <- 100
bzeros <- besselzero(0, N0, 1) # alpha_(nu=0,k=N0) and first kind

# notice that nu = (d-2)/2, hence d = 2 => nu = 0
True_coef <- function(ggg,N0,hhh){
  bzeros<-besselzero(0,N0,1)
  beta <- rep(0,N0)
  for(i in 1:N0){
    integrand <- function(x) ggg(x)*besselJ((x*bzeros[i]/hhh),nu = 0)*x
    denom <- 0.5*(hhh*besselJ(bzeros[i],nu = 1))^2
    beta[i] <- (integrate(integrand,0,hhh)$value)/(denom)
  }
  return(beta)
}

True_coef_bis <- function(r, rmin, K, R){
  bzeros <- besselzero(0, K, 1)
  beta <- rep(0, K)
  for(i in 1:K){
    num <- besselJ(((r-rmin)*bzeros[i]/R), nu = 0)
    denom <- sqrt(0.5)*(R*besselJ(bzeros[i], nu = 1))
    beta[i] <- num/denom
  }
  return(beta)
}

## Strauss point process with gamma=0.5 on W=square(2)
X.strauss <- rStrauss(beta = 100, gamma = 0.5, W = square(2))
X.grid <- gridcentres(window = square(2), nx = 10, ny = 10)
X.grid.ppp <- ppp(X.grid$x,X.grid$y, window = square(2))
plot(X.strauss)
points(X.grid, col = "red", pch = 16)


ddd <- crosspairs(X.strauss, X.grid.ppp, 0.5)
rr <- closepairs(X.strauss, 2, what = "all")

## Euclidean distance for two vectors
euclid.dist <- function(x,y){
  return(sqrt(sum((c(x[1],x[2]) - c(y[1],y[2]))^2)))
}

## Distance of a point u to all points of X
dist.all <- function(X, u){ # u=(u_1,u_2)
  n <- npoints(X)
  dt.u <- rep(NA, n) 
  for(i in 1:n){
    v <- c(X$x[i],X$y[i])
    dt.u[i] <- euclid.dist(v,u)
  }
  return(dt.u)
}

## ith coeffient of Fourier-Bessel basis
coef.ith <- function(r, rmin, bi, R) { # r is a vector
  lr <- length(r)
  s <- rep(NA, lr) 
  for(i in 1:lr) {
    num <- besselJ(((r[i]-rmin)*bi/R), nu = 0)
    denom <- sqrt(0.5)*(R*besselJ(bi, nu = 1))
    s[i] <- num/denom
  }
  return(sum(s))
}

## Covariates with Fourier-Bessel basis
covariates.fourierbessel <- function(X, rmin, K){
  n <- ceiling(sqrt(npoints(X))); Xw <- X$window
  X.grid <- gridcentres(window = Xw, nx = n, ny = n)
  m <- length(X.grid$x); dist.vu <- matrix(NA, nrow = m, ncol = npoints(X))
  R <- .5*max((Xw$xrange[2]-Xw$xrange[1]),(Xw$yrange[2]-Xw$xrange[1]))
  bzeros <- besselzero(0, K, 1); Sk.ux <- matrix(NA, nrow = K, ncol = 1)
  for(i in 1:m){
    u <- c(X.grid$x[i],X.grid$y[i])
    dist.vu[i,] <- dist.all(X, u)
    for(k in 1:K) {
      Sk.ux[k,] <- coef.ith(dist.vu[i,], rmin, bzeros[k], R)
    }
  }
  return(list(l1 = dist.vu, l2 = Sk.ux))
}


## Estimation
Sk_ux <- covariates.fourierbessel(X.strauss, 0, 100)
Sk_ux.im <- im(Sk_ux$l2, xrange = c(0,2), yrange = c(0,2))
theta.est <- ppm(X.strauss ~ Sk_ux.im)

#XX0 is a ppp object and ddd is the maximum distance
Test <- closepairs(XX0, ddd, what = "all")


window = owin(xrange = c(0,2), yrange = c(0,2))
X <- rpoispp(lambda = 10)
Y <- gridcentres(window = owin(xrange = c(0,1), yrange = c(0,1)), 1, 1)
Y <- ppp(Y$x,Y$y,window = owin(xrange = c(0,1), yrange = c(0,1)))

plot(X)
points(Y,col="red",pch=16)

ddd <- crosspairs(X, Y, 0.25)
ddd$i #index X
ddd$j #index Y
ddd$d #distance

sqrt(sum((c(X$x[6],X$y[6]) - c(Y$x[51],Y$y[51]))^2))

flag = ddd$i[(ddd$j==5)]
c(Y$x[5],Y$y[5])
cbind(c(X$x[flag]),c(X$y[flag]))

points(Y$x[5],Y$y[5],col = "green",pch = 16)


window = owin(xrange = c(0,2), yrange = c(0,2))
X <- rpoispp(lambda = 30, win = window)
Y <- gridcentres(window = window, 10, 10)
Y <- ppp(Y$x,Y$y,window = window)


#ddd$i #index X
#ddd$j #index Y
#ddd$d #distance


Create.basis <- function(X, Y, Rmax = 1, rmin = 0.001, n, xrange = c(0,1), yrange = c(0,1)){
  ddd <- crosspairs(X, Y, Rmax+rmin)
  bzeros <- besselzero(0,n,1)
  eval.basis <- function(i,n0){
    flag <- ddd$j == i
    ratio <- 0
    if(sum(flag) != 0){
      num <- besselJ(((ddd$d[flag]-rmin)*bzeros[n0]/Rmax), nu = 0)
      denom <- sqrt(0.5)*(Rmax*besselJ(bzeros[n0], nu = 1))
      ratio <- num/denom
    }
    return(cbind(i,sum(ratio))) 
  }
  Y.indexes <- matrix(1:Y$n, nrow = Y$n)
  basis.img <- list()
  names <- NULL
  for(i in 1:n){
    Phi.u <- apply(Y.indexes, 1, eval.basis, n0 = i)
    basis.img[[i]] <- im(matrix(Phi.u[2,], ncol = sqrt(Y$n), nrow = sqrt(Y$n)), xrange = xrange, yrange = yrange)
    names[i] <- paste("basis", i, collapse = "", sep="")
  }
  names(basis.img) <- names
  return(basis.img)
}

my.basis <- Create.basis(X, Y, Rmax = 0.1, rmin = 0.0001, 10)
plot(my.basis$basis10)

## Example
X_strauss <- rStrauss(beta = 100, gamma = 0.5, W = square(2))
window = owin(xrange = c(0,2), yrange = c(0,2))
Y <- gridcentres(window = window, 20, 20)
Y <- ppp(Y$x,Y$y,window = owin(xrange = c(0,2), yrange = c(0,2)))
Sk_ux <- Create.basis(X.strauss,Y,Rmax = 0.1, rmin = 0.0001, 10, xrange = c(0,2), yrange = c(0,2))
theta.est <- ppm(X_strauss ~. + 1, covariates = Sk_ux)

plot(predict(theta.est,Y))
