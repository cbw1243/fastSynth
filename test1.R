library(DEoptim)
library(nloptr)

## Rosenbrock Banana function and gradient in separate functions
eval_f <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}

eval_grad_f <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
             200 * (x[2] - x[1] * x[1])) )
}


# initial values
x0 <- c( -1.2, 1 )

opts <- list("algorithm"="NLOPT_LD_LBFGS",
             "xtol_rel"=1.0e-8)


# solve Rosenbrock Banana function
res <- nloptr( x0=x0,
               eval_f=eval_f,
               eval_grad_f=eval_grad_f,
               opts=opts)
print( res )

res$solution

res2 <- nloptr( x0=x0,
               eval_f=eval_f,
               opts= list("algorithm"="NLOPT_LN_SBPLX",
                          "xtol_rel"=1.0e-8))
res2$solution





rastrigin <- function(x) 10*length(x)+sum(x^2-10*cos(2*pi*x))

est.ras <- DEoptim(rastrigin, lower = c(-5, -5), upper = c(5, 5), control = list(storepopfrom = 1, trace = FALSE))

est.ras$optim

rastrigin(c(1, 1))


genrose.f <- function(x){
    n <- length(x)
    fval <- 1.0 + sum (100 * (x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
    return(fval)
    }

n <- 10

microbenchmark(
ans <- DEoptim(fn=genrose.f, lower=rep(-5, n), upper=rep(5, n),
                control=list(NP=100, itermax=4000,trace=FALSE)),

ans1 <- optim(par=runif(10,-5,5), fn=genrose.f, method="BFGS",
               control=list(maxit=4000)),
res <- nloptr( x0=x0,
               eval_f=eval_f,
               eval_grad_f=eval_grad_f,
               opts=opts)
times = 2)

ans1$par

c(1, 1, 1)

microbenchmark
  ans1 <- optim(par=c(1, 1, 1), fn=genrose.f, method="Nelder-Mead",
                control=list(maxit=4000, reltol = 1.0e-8))
  res <- nloptr( x0=c(1, 1, 1),
                 eval_f=genrose.f,
                 opts=list("algorithm"="NLOPT_LN_NELDERMEAD",
                           "xtol_rel"=1.0e-8))
  times = 2)

ans1$par
ans1$value
res$solution
res$objective



nloptr.print.options()


optimFunc <- function(x) x[1]^2 + x[2]^2 + x[1] + 2*x[2]

ans1 <- optim(par=c(1, 1), fn=optimFunc, method="Nelder-Mead",
              control=list(maxit=4000, reltol = 1.0e-8))
ans1$par

res <- nloptr( x0=c(1, 1),
               eval_f=optimFunc,
               opts=list("algorithm"="NLOPT_LN_NELDERMEAD",
                         "xtol_rel"=1.0e-8))
res$solution

microbenchmark(
ans1 <- optim(par=c(1, 1), fn=optimFunc, method="Nelder-Mead",
              control=list(maxit=4000, reltol = 1.0e-8)),
res <- nloptr( x0=c(1, 1),
               eval_f=optimFunc,
               opts=list("algorithm"="NLOPT_LN_NELDERMEAD",
                         "xtol_rel"=1.0e-8)),
times = 1)
