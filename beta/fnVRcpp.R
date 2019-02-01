fn.VRcpp <- '
arma::vec variablesvRcpp = Rcpp::as<arma::vec>(variablesv);
arma::vec diagElement = abs(variablesvRcpp)/sum(abs(variablesvRcpp));
arma::mat V = arma::diagmat(diagElement, 0);

arma::mat X0scaledRcpp = Rcpp::as<arma::mat>(X0scaled);
arma::mat H = trans(X0scaledRcpp) * V * X0scaledRcpp;
//Rcpp::Rcout << H << std::endl;

arma::mat a = Rcpp::as<arma::mat>(X1scaled);
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
arma::mat Z0Rcpp = Rcpp::as<arma::mat>(Z0);
arma::mat Z1Rcpp = Rcpp::as<arma::mat>(Z1);

arma::vec lossw = trans(a - X0scaledRcpp * out2) * V * (a - X0scaledRcpp * out2);

arma::vec lossv = trans(Z1Rcpp - Z0Rcpp * out2) * (Z1Rcpp - Z0Rcpp * out2);

lossv = lossv/Z0Rcpp.n_rows;
return(wrap(lossv));
//return Rcpp::List::create(Rcpp::Named("c")=c1,
//                          Rcpp::Named("H")=H,
//                          Rcpp::Named("A")=A,
//                          Rcpp::Named("b")=b, Rcpp::Named("l")=l,Rcpp::Named("u")=u, Rcpp::Named("r")=r,
//                          Rcpp::Named("k")=k);
'

x <- cxxfunction(signature(variablesv = 'numeric', X0scaled = 'matrix', X1scaled = 'matrix',
                           Z0 = 'matrix', Z1 = 'matrix'), fn.VRcpp, "RcppArmadillo")



x(c(1, 0, 1), as.matrix(X0.scaled), as.matrix(X1.scaled), as.matrix(s.data$Z0), as.matrix(s.data$Z1))


