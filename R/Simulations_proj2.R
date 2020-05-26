Y <- gridcentres(window = window, 10, 10)
Y <- ppp(Y$x, Y$y, window = window)
W1 <- square(1); W2 <- square(2); W3 <- square(3)

## Strauss with Fourier-Bessel basis
 Sim.Strauss_FB_W1 <- FSim(Y=Y,r=0.01,R=0.05,nsim=1000,gam=0.2,npts=500,K=10,f.dummy=16,
                        int.rg=0.07,window=W1)
 Sim.Strauss_FB_W2 <- FSim(Y=Y,r=0.01,R=0.05,nsim=1000,gam=0.2,npts=500,K=10,f.dummy=16,
                        int.rg=0.07,window=W2)
 Sim.Strauss_FB_W3 <- FSim(Y=Y,r=0.01,R=0.05,nsim=1000,gam=0.2,npts=500,K=10,f.dummy=16,
                           int.rg=0.07,window=W3)

## Softcore with Fourier-Bessel basis
Sim.Softcore_FB_W1 <- FSim(Model='sftcr',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=0.5,sig=0.1,
                        window=W1)
Sim.Softcore_FB_W2 <- FSim(Model='sftcr',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=0.5,sig=0.1,
                           window=W2)
Sim.Softcore_FB_W3 <- FSim(Model='sftcr',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=0.5,sig=0.1,
                           window=W3)

## Diggle, Gates, and Stibbard with Fourier-Bessel basis 
Sim.DGS_FB_W1 <- FSim(Model='dgs',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,Rho=0.08,window=W1)
Sim.DGS_FB_W2 <- FSim(Model='dgs',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,Rho=0.08,window=W2)
Sim.DGS_FB_W3 <- FSim(Model='dgs',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,Rho=0.08,window=W3)

## Diggle-Gratton with fouier-Bessel basis
Sim.DigGrat_FB_W1 <- FSim(Model='diggra',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=3,Delta=0.02,
                       Rho=0.04,window=W1)
Sim.DigGrat_FB_W2 <- FSim(Model='diggra',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=3,Delta=0.02,
                          Rho=0.04,window=W2)
Sim.DigGrat_FB_W3 <- FSim(Model='diggra',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=3,Delta=0.02,
                          Rho=0.04,window=W3)

## Fiksel with Fourier-Bessel basis
Sim.Fiksel_FB_W1 <- FSim(Model='fiksel',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=2,Hc=0.07,
                      int.rg=0.5,a=-1.0,window=W1)
Sim.Fiksel_FB_W2 <- FSim(Model='fiksel',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=2,Hc=0.07,
                         int.rg=0.5,a=-1.0,window=W2)
Sim.Fiksel_FB_W3 <- FSim(Model='fiksel',Y=Y,r=0.01,R=0.05,nsim=100,npts=500,K=5,f.dummy=16,kap=2,Hc=0.07,
                         int.rg=0.5,a=-1.0,window=W3)
