# Test on the new function for synthetical control.
# --------------------------------------------------------------------------------------------------------------#
# Part 1. Prepare for running the script.
# --------------------------------------------------------------------------------------------------------------#
rm(list = ls())
# Prepare the library
library(data.table)
library(nloptr)
library(kernlab)
library(microbenchmark)
library(Synth)

source('fastSynthV1_NLopt.r')

data <- fread('test.csv')
J   = 30
t0  = 30

s.data <- dataprep(foo = data,
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

# --------------------------------------------------------------------------------------------------------------#
# Part 2. Set up the parameters, and run the codes within the function.
# --------------------------------------------------------------------------------------------------------------#
{
# Benchmark testing.
quadopt = "ipop"
Margin.ipop = 0.0005
Sigf.ipop = 5
Bound.ipop = 10
custom.v = NULL
optimxmethod = c("Nelder-Mead")
data.prep.obj <- dataprep(foo = data,
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

X1 <- data.prep.obj$X1
Z1 <- data.prep.obj$Z1
X0 <- data.prep.obj$X0
Z0 <- data.prep.obj$Z0

# Normalize X
nvarsV <- dim(X0)[1]
big.dataframe <- cbind(X0, X1)
divisor <- sqrt(apply(big.dataframe, 1, var))
scaled.matrix <-
  t(t(big.dataframe) %*% ( 1/(divisor) *
                             diag(rep(dim(big.dataframe)[1], 1)) ))

X0.scaled <- scaled.matrix[,c(1:(dim(X0)[2]))]
if(is.vector(X0.scaled)==TRUE)
{X0.scaled <- t(as.matrix(X0.scaled))}
X1.scaled <- scaled.matrix[,dim(scaled.matrix)[2]]


SV1 <- rep(1/nvarsV,nvarsV)

microbenchmark(
rgV.optim.1 <- nloptr( x0=SV1,
                       eval_f=fn.V,
                       opts=list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                 "xtol_rel"=1.0e-5),
                       X0.scaled = X0.scaled,
                       X1.scaled = X1.scaled,
                       Z0 = Z0,
                       Z1 = Z1,
                       quadopt = quadopt,
                       margin.ipop = Margin.ipop,
                       sigf.ipop = Sigf.ipop,
                       bound.ipop = Bound.ipop),

rgV.optim.2 <- optimx::optimx(par=SV1, fn=fn.V,
                              gr=NULL, hess=NULL,
                              method=optimxmethod, itnmax=NULL, hessian=FALSE,
                              X0.scaled = X0.scaled,
                              X1.scaled = X1.scaled,
                              Z0 = Z0,
                              Z1 = Z1,
                              quadopt = quadopt,
                              margin.ipop = Margin.ipop,
                              sigf.ipop = Sigf.ipop,
                              bound.ipop = Bound.ipop),
times = 1)
}

# --------------------------------------------------------------------------------------------------------------#
# Part 3. Compare the new function with the benchmark function: synth.
# --------------------------------------------------------------------------------------------------------------#
# Compare the computation time.
microbenchmark(
scm <- synth(data.prep.obj = s.data, optimxmethod = "BFGS", quadopt = "ipop", Margin.ipop = 0.001), # default is 0.0005.
scm2 <- fastSynthV1(data.prep.obj = s.data, optimxmethod = "NLOPT_LN_NELDERMEAD", quadopt = "ipop", Margin.ipop = 0.001),
times = 1)

# Run the new function alone.
scm2 <- fastSynthV1(data.prep.obj = s.data, optimxmethod = "NLOPT_LN_NELDERMEAD", quadopt = "ipop", Margin.ipop = 0.001)
