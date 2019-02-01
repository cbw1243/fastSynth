# Call the functions of rnorm from R.
src <- '
Function rnorm("rnorm");
double meanRcpp = Rcpp::as<double>(mean);

Rcpp::NumericVector x = rnorm(1, _["mean"] = meanRcpp, _["sd"] = 3.2 );
return wrap(x);
'

x <- cxxfunction(signature(mean = 'numeric'), src, "RcppArmadillo")
x(1)




