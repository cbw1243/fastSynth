#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
S4 ipopInRcpp(NumericMatrix c, NumericMatrix H, NumericVector A, NumericVector b,
                         NumericMatrix l, NumericMatrix u, NumericVector r){

  // Obtaining namespace of Matrix package
  Environment pkg = Environment::namespace_env("kernlab");

  // Picking up Matrix() function from Matrix package
  Function f = pkg["ipop"];

  // Executing Matrix( m, sparse = TRIE )
  return f( c, H, A, b, l, u, r,
            Named("sigf", 7),  Named("maxiter", 40), Named("margin", 0.05), Named("bound", 10), Named("verb", 0));
}
