fn.VRcppFunc <- '
arma::vec objectiveFunc(arma::vec variablesvRcpp, arma::mat X0scaledRcpp, arma::mat X1scaledRcpp, arma::mat Z0Rcpp, arma::mat Z1Rcpp){
arma::vec diagElement = abs(variablesvRcpp)/sum(abs(variablesvRcpp));
arma::mat V = arma::diagmat(diagElement, 0);

arma::mat H = trans(X0scaledRcpp) * V * X0scaledRcpp;
//Rcpp::Rcout << H << std::endl;

arma::mat a = X1scaledRcpp;
// Rcpp::Rcout << a.size() << std::endl;
// Rcpp::Rcout << a << std::endl;

arma::mat c2 = -1* (trans(a) * V * X0scaledRcpp);
arma::colvec c1 = arma::conv_to<arma::colvec>::from(c2);
int k = c2.n_cols;

arma::rowvec A = arma::ones<arma::rowvec>(k);

arma::vec l = arma::zeros<arma::vec>(k);
arma::vec u = arma::ones<arma::vec>(k);
int b = 1;
int r = 0;
Rcpp::Environment kernlab("package:kernlab");
Rcpp::Function ipopFunc = kernlab["ipop"];

Rcpp::S4 out = ipopFunc(_["c"] = c1, _["H"] = H, _["A"] = A, _["b"] = b, _["l"] = l, _["u"] = u,_["r"] = r, _["sigf"] = 5,
_["maxiter"] = 1000, _["margin"] = 0.005, _["bound"] = 10, _["verb"] = 0);

Rcpp::NumericVector out1 = Rcpp::NumericVector(out.slot("primal"));

arma::vec out2 = Rcpp::as<arma::vec>(out1);

//Rcpp::Rcout << "Okay" << std::endl;

arma::vec lossw = arma::trans(a - X0scaledRcpp * out2) * V * (a - X0scaledRcpp * out2);

arma::vec lossv = arma::trans(Z1Rcpp - Z0Rcpp * out2) * (Z1Rcpp - Z0Rcpp * out2);

lossv = lossv/Z0Rcpp.n_rows;
return lossv;
//return(wrap(lossv));
//return Rcpp::List::create(Rcpp::Named("c")=c1,
//                          Rcpp::Named("H")=H,
//                          Rcpp::Named("A")=A,
//                          Rcpp::Named("b")=b, Rcpp::Named("l")=l,Rcpp::Named("u")=u, Rcpp::Named("r")=r,
//                          Rcpp::Named("k")=k);
    }
'

maxFunc <- '
arma::vec variablesvRcpp = Rcpp::as<arma::vec>(variablesv);
arma::mat X0scaledRcpp = Rcpp::as<arma::mat>(X0scaled);
arma::mat X1scaledRcpp = Rcpp::as<arma::mat>(X1scaled);
arma::mat Z0Rcpp = Rcpp::as<arma::mat>(Z0);
arma::mat Z1Rcpp = Rcpp::as<arma::mat>(Z1);

Rcpp::Environment stats("package:stats");
Rcpp::Function optim = stats["optim"];
Rcpp::Environment optimx("package:optimx");
Rcpp::Function optimxFunc = optimx["optimx"];


Rcpp::List opt_results = optim(Rcpp::_["par"]    = variablesv,
                               Rcpp::_["fn"]     = Rcpp::InternalFunction(&objectiveFunc),
                               Rcpp::_["method"] = "BFGS",
                               Rcpp::_["X0scaledRcpp"] = X0scaledRcpp,
                               Rcpp::_["X1scaledRcpp"] = X1scaledRcpp,
                               Rcpp::_["Z0Rcpp"] = Z0Rcpp,
                               Rcpp::_["Z1Rcpp"] = Z1Rcpp);

//Rcpp::List opt_results = optimxFunc(Rcpp::_["par"]    = variablesv,
//                                    Rcpp::_["fn"]     = Rcpp::InternalFunction(objectiveFunc),
//                                    Rcpp::_["method"] = "BFGS",
//                                    Rcpp::_["X0scaledRcpp"] = X0scaledRcpp,
//                                    Rcpp::_["X1scaledRcpp"] = X1scaledRcpp,
//                                    Rcpp::_["Z0Rcpp"] = Z0Rcpp,
//                                   Rcpp::_["Z1Rcpp"] = Z1Rcpp);
// arma::vec opt_results = objectiveFunc(variablesvRcpp, X0scaledRcpp, X1scaledRcpp, Z0Rcpp, Z1Rcpp);

// Extract out the estimated parameter values
// arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);

// Return estimated values
return wrap(opt_results);
'
optimFunc <- cxxfunction(signature(variablesv = 'numeric', X0scaled = 'matrix', X1scaled = 'matrix',
                            Z0 = 'matrix', Z1 = 'matrix'), maxFunc, includes = fn.VRcppFunc,  "RcppArmadillo")

test1 <- optimFunc(variablesv = c(0.33, 0.33, 0.33), X0scaled = as.matrix(X0.scaled), X1scaled = as.matrix(X1.scaled), as.matrix(Z0), as.matrix(Z1))
test1



# library(microbenchmark)
#
# microbenchmark(
# test1 <- x2(variablesv = c(0.33, 0.33, 0.33), X0scaled = as.matrix(X0.scaled), X1scaled = as.matrix(X1.scaled), as.matrix(Z0), as.matrix(Z1)),
# rgV.optim.1 <- optimx::optimx(par=c(0.33, 0.33, 0.33), fn=fn.V,
#                               gr=NULL, hess=NULL,
#                               method='BFGS', itnmax=NULL, hessian=FALSE,
#                               X0.scaled = X0.scaled,
#                               X1.scaled = X1.scaled,
#                               Z0 = Z0,
#                               Z1 = Z1,
#                               quadopt = "ipop",
#                               margin.ipop = 0.0005,
#                               sigf.ipop = 5,
#                               bound.ipop = 10),
# times  = 1)


