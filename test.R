#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void function01(){

  // Obtain environment containing function
  Rcpp::Environment package_env("package:package_name_here");

  // Make function callable from C++
    Rcpp::Function rfunction = package_env["function_name"];

    // Call the function and receive output (might not be list)
    Rcpp::List test_out = rfunction(....);

}



test <- '
Rcpp::Environment package_env("package:package_name_here");

return(wrap(u));
'
#int b = 1;
#int r = 0
x <- cxxfunction(signature(variablesv = 'numeric', X0scaled = 'matrix', X1scaled = 'matrix'), fn.VRcpp, "Rcpp")

set.seed(42)
x <- rnorm(1e5)
fivenum(x)

res <- kernlab::ipop(c = c, H = H, A = A, b = b, l = l, u = u,
                     r = r, bound = bound.ipop, margin = margin.ipop,
                     maxiter = 1000, sigf = sigf.ipop)

callFunction(x, testFunc)
str(x)
testFunc(x)
testFunc <- function(x) sd(x)

sourceCpp("ipopRcpp.cpp")



src <- '
arma::mat cRcpp = Rcpp::as<arma::mat>(c);
arma::mat HRcpp = Rcpp::as<arma::mat>(H);
arma::mat ARcpp = Rcpp::as<arma::mat>(A);
arma::vec bRcpp = Rcpp::as<arma::vec>(b);
arma::mat lRcpp = Rcpp::as<arma::mat>(l);
arma::mat uRcpp = Rcpp::as<arma::mat>(u);
arma::vec rRcpp = Rcpp::as<arma::vec>(r);

List out = ipopInRcpp(wrap(cRcpp), wrap(HRcpp), wrap(ARcpp), wrap(bRcpp), wrap(lRcpp), wrap(uRcpp), wrap(rRcpp));
return(wrap(out))
'

x <- cxxfunction(signature(c = 'matrix',  H = 'matrix', A = 'matrix', b = 'numeric',
                           l = 'matrix', u = 'matrix', r = 'numeric'), src, "RcppArmadillo")


?rnorm

src <- '
int nRcpp = Rcpp::as<int>(n);
double meanRcpp = Rcpp::as<double>(mean);
double sdRcpp = Rcpp::as<double>(sd);

Rcpp::Environment stats("package:stats");
Rcpp::Function Rfunc = stats["rnorm"];

double opt_results = Rfunc(Rcpp::_["n"] = nRcpp, Rcpp::_["mean"] = meanRcpp, Rcpp::_["sd"] = sdRcpp);

return wrap(opt_results);
'

x <- cxxfunction(signature(n = 'integer',  mean = 'numeric', sd = 'numeric'), src, "RcppArmadillo")





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

# Evaluate on the synth(...) function. Return character Error if it encounters error.
microbenchmark(
scm1          = synth(data.prep.obj = s.data, optimxmethod = "BFGS", quadopt = "ipop", Margin.ipop = 0.005),
scm          = fastSynth(data.prep.obj = s.data, optimxmethod = "BFGS"),
times = 5)

scm1$solution.v

scm$solution.v


w.scm        = as.matrix(scm$solution.w)  ##weigting vector for SCM
y0.post      = as.matrix(ydata[(t0+1):T,2:(J+1)])
y.scm     = y0.post%*%w.scm   ##counterfactuls for SCM
diff.scm    =abs(y1.post-y.scm)
mean(diff.scm^2); mean(diff.scm); mean(abs(y1.post-y.scm)/abs(y1.post))

w.scm1        = as.matrix(scm1$solution.w)  ##weigting vector for SCM
y0.post1      = as.matrix(ydata[(t0+1):T,2:(J+1)])
y.scm1     = y0.post1%*%w.scm1   ##counterfactuls for SCM
diff.scm1    =abs(y1.post-y.scm1)
mean(diff.scm1^2); mean(diff.scm1); mean(abs(y1.post-y.scm1)/abs(y1.post))

