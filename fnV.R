# An objective function to be optimized, directly from the Synth package with minor modification.

fn.V <-
  function (variables.v = stop("variables.v missing"), X0.scaled = stop("X0.scaled missing"),
          X1.scaled = stop("X1.scaled missing"), Z0 = stop("Z0 missing"),
          Z1 = stop("Z1 missing"), margin.ipop = 5e-04, sigf.ipop = 5,
          bound.ipop = 10, quadopt = "ipop")
{
  Check <- sum(quadopt %in% c("ipop", "LowRankQP"))
  if (Check != 1) {
    stop("option quadopt must be one of ipop or LowRankQP")
  }
  V <- diag(x = as.numeric(abs(variables.v)/sum(abs(variables.v))),
            nrow = length(variables.v), ncol = length(variables.v))
  H <- t(X0.scaled) %*% V %*% (X0.scaled)
  a <- X1.scaled
  c <- -1 * c(t(a) %*% V %*% (X0.scaled))
  A <- t(rep(1, length(c)))
  b <- 1
  l <- rep(0, length(c))
  u <- rep(1, length(c))
  r <- 0
  if (quadopt == "ipop") {
    res <- kernlab::ipop(c = c, H = H, A = A, b = b, l = l, u = u,
                r = r, bound = bound.ipop, margin = margin.ipop,
                maxiter = 1000, sigf = sigf.ipop)
    solution.w <- as.matrix(kernlab::primal(res))
    #print(solution.w)
  }
  else {
    if (quadopt == "LowRankQP") {
      require(LowRankQP)
      res <- LowRankQP(Vmat = H, dvec = c, Amat = A, bvec = 1,
                       uvec = rep(1, length(c)), method = "LU")
      solution.w <- as.matrix(res$alpha)
    }
  }
  loss.w <- as.numeric(t(X1.scaled - X0.scaled %*% solution.w) %*%
                         (V) %*% (X1.scaled - X0.scaled %*% solution.w))
  loss.v <- as.numeric(t(Z1 - Z0 %*% solution.w) %*% (Z1 -
                                                        Z0 %*% solution.w))
  loss.v <- loss.v/nrow(Z0)
  #return(invisible(loss.v))
  return(loss.v)
}


