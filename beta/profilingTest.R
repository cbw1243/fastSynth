profvis::profvis({
rgV.optim.1 <- optimx::optimx(par=SV1, fn=fn.V,
                              gr=NULL, hess=NULL,
                              method=optimxmethod, itnmax=NULL, hessian=FALSE,
                              control=list(kkt=FALSE,
                                           starttests=FALSE,
                                           dowarn=FALSE,
                                           all.methods=all.methods),
                              X0.scaled = X0.scaled,
                              X1.scaled = X1.scaled,
                              Z0 = Z0,
                              Z1 = Z1,
                              quadopt = quadopt,
                              margin.ipop = Margin.ipop,
                              sigf.ipop = Sigf.ipop,
                              bound.ipop = Bound.ipop)

test <- optim(par=SV1, fn=fn.V,
              X0.scaled = X0.scaled,
              X1.scaled = X1.scaled,
              Z0 = Z0,
              Z1 = Z1,
              quadopt = quadopt,
              margin.ipop = Margin.ipop,
              sigf.ipop = Sigf.ipop,
              bound.ipop = Bound.ipop)
})

source('fnV.r')
fn.V(c(test$par), X0.scaled = X0.scaled,
     X1.scaled = X1.scaled,
     Z0 = Z0,
     Z1 = Z1)
