#include <Rcpp.h>

using namespace Rcpp;

// sigmoid_rcpp
double sigmoid_rcpp(double x);
RcppExport SEXP varbvs_sigmoid_rcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmoid_rcpp(x));
    return rcpp_result_gen;
END_RCPP
}

// varbvsnormupdate_rcpp
void varbvsnormupdate_rcpp(const NumericMatrix& X, double sigma, double sa, double b0, const NumericVector& logodds, const NumericVector& xy, const NumericVector& d, NumericVector& alpha, NumericVector& mu, NumericVector& Xr, const IntegerVector& i);
RcppExport SEXP varbvs_varbvsnormupdate_rcpp(SEXP XSEXP, SEXP sigmaSEXP, SEXP saSEXP, SEXP b0SEXP, SEXP logoddsSEXP, SEXP xySEXP, SEXP dSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP XrSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type sa(saSEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type i(iSEXP);
    varbvsnormupdate_rcpp(X, sigma, sa, b0, logodds, xy, d, alpha, mu, Xr, i);
    return R_NilValue;
END_RCPP
}
