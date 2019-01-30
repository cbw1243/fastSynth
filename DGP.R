# The codes below generate data for doing program evaluation.

# DGP
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
