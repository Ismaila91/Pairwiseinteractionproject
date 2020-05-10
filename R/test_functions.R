####------------- Test functions on [0,1]

##----1. Uniform
curve(dunif(x))
##----2. Normal
curve(dnorm(x,0.5,0.15),0,1)
##----3. Bimodal
curve(0.5*dnorm(x,0.4,0.12)+0.5*dnorm(x,0.7,0.08),0,1)
##----4. Strata
curve(0.5*dnorm(x,0.2,0.06)+0.5*dnorm(x,0.7,0.08),0,1)
##----5. Delta
curve(dnorm(x,0.5,0.02),0,1)
##----6. Angle
f6 <- function(x){
  cond <- ((x>=0)&(x<=0.5))
  return (ifelse(cond,(1/0.16095)*dnorm(x,1,0.7),(1/0.16095)*dnorm(x,0,0.7)))
}
curve(f6(x),0,1)
##----7. Monotone
integrand <- function(u) {
  return (dnorm(u, mean = 2, sd = 0.8))
}
f7 <- function(x) {
  denom <- integrate(integrand,0,1)$value
  return (dnorm(x,2,0.8)/denom)
}
curve(f7(x),0,1)
##----8. Steps
f8 <- function(x) {
  cond1 <- ((x>=0)&(x<1/3))
  cond2 <- ((x>=1/3)&(x<3/4))
  return (ifelse(cond1, 0.6, ifelse(cond2, 0.9, 204/120)))
}
curve(f8(x),0,1)

##################################################################################################

####----------------- Cosine basis on [0,1]

##---- 1. First four elements of the cosine system
par(mfrow=c(1,4))
curve(dunif(x),xlab="",ylab="", main="0"); curve(sqrt(2)*cos(pi*x),0,1,xlab="",ylab="", main="1")
curve(sqrt(2)*cos(2*pi*x),0,1,xlab="",ylab="", main="2"); curve(sqrt(2)*cos(3*pi*x),0,1,xlab="",ylab="", main="3")
par(mfrow=c(1,1))

##---- 2. Elements of the cosine basis
jth_el_cosine_basis <- function(x,j) { 
  cond <- (j==0)
  return (ifelse(cond, 1, sqrt(2)*cos(pi*j*x)))
}
elements_cosine_basis <- function(x, j) { 
  n <- length(x)
  el <- rep(NA, n)
  for(i in 1:n){
    el[i] <- jth_el_cosine_basis(x[i],j)
  }
  return(el)
}

##---- 3. Approximation of the Normal density by the cosine basis
Integrand_Normal <- function(x, j){
  if(j==0) return(dnorm(x,0.5,0.15))
  else{
    return(dnorm(x,0.5,0.15)*sqrt(2)*cos(j*pi*x))
  }
}

##---- 4. The Fourier coefficients for Normal density and cosine basis
FourierCoef_Normal_cosine_basis <- function(J) {
  f.coef <- rep(NA, J+1)
  for(i in 0:J){
    f.coef[i+1] <- integrate(Integrand_Normal, lower=0, upper=1, j=i)$value
  }
  return(f.coef)
}

##---- 5. Approximation of the Normal density via cosine basis expansion
Normal_Approx_cosine_basis <- function(x, J) {
  F.coef <- FourierCoef_Normal_cosine_basis(J)
  ssum <- 0
  for(i in 0:J){
    ssum <- F.coef[i+1]*elements_cosine_basis(x,i) + ssum
  }
  return(ssum)
}
curve(dnorm(x,0.5,0.15), main="Normal")
curve(Normal_Approx_cosine_basis(x,J=3),0,1,lty="dashed",add=TRUE)
curve(Normal_Approx_cosine_basis(x,J=5),0,1,lty="dotted",add=TRUE)
curve(Normal_Approx_cosine_basis(x,J=10),0,1,lty="longdash",add=TRUE)


