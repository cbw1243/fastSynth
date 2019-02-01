# The codes below generate data for doing program evaluation.
rm(list = ls())
library(dplyr)
library(data.table)
set.seed(12345)
# DGP

{
itr = 1000
#case 1
J   = 30
N   = 1+J
t0  = 30
t1  = 10
T   = t0+t1
b1=1
b2=2
# Suppose there are two factors
# factor loading
lam1 = rnorm(N,mean = 0,sd = 1)
lam2 = rnorm(N,mean = 0,sd = 1)
lam3 = rnorm(N,mean = 0,sd = 1)
rho = runif(N, min = 0, max = 1)
rho1 = runif(N, min = 0.1, max = 0.9)
rho2 = runif(N, min = 0.1, max = 0.9)
c1 = runif(N, min = 1, max = 2)
c2 = runif(N, min = 1, max = 2)
sigeta = 0.5*rchisq(N,1)+0.5

#SCM by Abadie (2010)
yhat.scm=matrix(0,nrow=itr,ncol=t1 )
mse.scm = matrix(0,nrow=itr,ncol=1)
mer.scm = matrix(0,nrow=itr,ncol=1)
map.scm = matrix(0,nrow=itr,ncol=1)

#true outcome
y1ture   =matrix(0,nrow=itr,ncol=t1)


# Simulation.
it <- 1

  # factors
  f1   = matrix(rnorm(T,mean = 0,sd = 1),nrow = T,ncol = 1)
  f2   = matrix(rnorm(T,mean = 0,sd = 1),nrow = T,ncol = 1)
  f3   = matrix(rnorm(T,mean = 0,sd = 1),nrow = T,ncol = 1)
  # The idiosyncratic error
  eit    = matrix(0,nrow = T, ncol = N)
  eta    = matrix(rnorm(N*T,mean = 0,sd = sigeta),nrow = T,ncol = N)
  eit[,1] = (1+1)*eta[,1]+1*eta[,2]
  for (i in 2:(N-1)){
    eit[,i] = (1+1)*eta[,i]+1*eta[,i-1]+1*eta[,i+1]
  }
  eit[,N] = (1+1)*eta[,N]+1*eta[,N-1]

  # The time varying coveriate
  x1  = matrix(0,nrow = T,ncol = N)
  x1[1,] = rnorm(N,mean = 0,sd = 1)
  x2  = matrix(0,nrow = T,ncol = N)
  x2[1,] = rnorm(N,mean = 0,sd = 1)
  for (i in 1:N){
    for (s in 2:T){
      x1[s,i] = 1+rho1[i]*x1[s-1,i]+rnorm(1,mean = 0,sd = 1)
      x2[s,i] = 1+rho2[i]*x2[s-1,i]+rnorm(1,mean = 0,sd = 1)
    }
  }
  # The variable y
  # assume y = beta*x+lambda*f+e, where beta=1
  y  = x1*b1+x2*b2+cbind(f1,f2,f3)%*%rbind(t(lam1),t(lam2),t(lam3))+eit
  #control unit
  y0.post      = as.matrix(y[(t0+1):T,2:(J+1)])
  x10.post      = as.matrix(x1[(t0+1):T,2:(J+1)])
  x20.post      = as.matrix(x2[(t0+1):T,2:(J+1)])
  #treatment unit
  y1.post      = as.matrix(y[(t0+1):T,1])
  x11.post      = as.matrix(x1[(t0+1):T,1])
  x21.post      = as.matrix(x2[(t0+1):T,1])

  y1ture[it,]  = t(y1.post)

  y1.pre      = as.matrix(y[1:t0, 1])
  ytau.pre    = as.matrix(y[1:t0, 2:(J+1)])
  ytau.post   = as.matrix(y[(t0+1):T, 2:(J+1)])
  x11.pre     = as.matrix(x1[1:t0, 1])
  x11.post    = as.matrix(x1[(t0+1):T, 1])
  x1tau.pre   = as.matrix(x1[1:t0, 2:(J+1)])

  x21.pre     = as.matrix(x2[1:t0, 1])
  x21.post    = as.matrix(x2[(t0+1):T, 1])
  x2tau.pre   = as.matrix(x2[1:t0, 2:(J+1)])


  ########################################
  ##SCM
  ydata = as.matrix(y)
  nam   = as.matrix(seq(1,(J+1),1))%x%rep(1,T)      # Identifier Column
  yr0   = as.matrix(seq(1,T,1))%*%rep(1,(J+1))      # Year Matrix
  yr    = matrix(yr0,nrow = T*(J+1),ncol=1)         # Year Column

  ycol  = matrix(ydata,nrow = T*(J+1), ncol = 1)    # Dependent Column
  x1col  = matrix(x1,nrow = T*(J+1), ncol = 1)    # X1
  x2col  = matrix(x2,nrow = T*(J+1), ncol = 1)    # X2
  xx=cbind(x1col,x2col)

  data           = as.data.frame(cbind(nam,yr,ycol,xx))
  colnames(data) = c("id","time","ycol","x1", 'x2')
}
write.csv(data, 'test.csv', row.names = F)


#
rm(list = ls())

data <- fread('huayanDGP.csv')
J   = 30
t0  = 30
s.data  = dataprep(foo = data,
                     predictors = c("x1", "x2"),
                     predictors.op = "mean",
                     time.predictors.prior = 1:t0,
                     dependent = "ycol",
                     unit.variable = "id",
                     time.variable = "time",
                     special.predictors = list(list("ycol", t0, "mean")),
                     treatment.identifier = 1,
                     controls.identifier = c(2:(J+1)),
                     time.optimize.ssr = 1:t0,
                     time.plot = 1:T)

# Compare the computation time.
microbenchmark(
  scm <- synth(data.prep.obj = s.data, optimxmethod = "BFGS", quadopt = "ipop", Margin.ipop = 0.001), # default is 0.0005.
  scm2 <- fastSynthV1(data.prep.obj = s.data, optimxmethod = "NLOPT_LN_NELDERMEAD", Margin.ipop = 0.001),
  times = 1)

system.time( scm2 <- fastSynthV1(data.prep.obj = s.data, optimxmethod = "NLOPT_LN_NELDERMEAD", Margin.ipop = 0.001))
microbenchmark(scm2 <- fastSynthV1(data.prep.obj = s.data, optimxmethod = "NLOPT_LN_NELDERMEAD", Margin.ipop = 0.001), times = 10)



