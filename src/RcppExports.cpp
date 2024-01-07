// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// invert
arma::cx_mat invert(arma::cx_mat x);
RcppExport SEXP _complex_invert(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(invert(x));
    return rcpp_result_gen;
END_RCPP
}
// polyprodcomplex
ComplexVector polyprodcomplex(ComplexVector const& poly1, ComplexVector const& poly2);
RcppExport SEXP _complex_polyprodcomplex(SEXP poly1SEXP, SEXP poly2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ComplexVector const& >::type poly1(poly1SEXP);
    Rcpp::traits::input_parameter< ComplexVector const& >::type poly2(poly2SEXP);
    rcpp_result_gen = Rcpp::wrap(polyprodcomplex(poly1, poly2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_complex_invert", (DL_FUNC) &_complex_invert, 1},
    {"_complex_polyprodcomplex", (DL_FUNC) &_complex_polyprodcomplex, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_complex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
