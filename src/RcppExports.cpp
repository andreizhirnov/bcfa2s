// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_bcfa_estim_spa
Rcpp::List rcpp_bcfa_estim_spa(Rcpp::NumericMatrix xY, Rcpp::IntegerMatrix xlam_ad, Rcpp::NumericVector xlambdabar, Rcpp::NumericMatrix xVlambda, Rcpp::NumericMatrix xL, Rcpp::NumericVector xalpha0, Rcpp::NumericVector xbeta0, int xintercepts, int xchains, int xkeep, int xstart, int xend);
RcppExport SEXP _bcfa2s_rcpp_bcfa_estim_spa(SEXP xYSEXP, SEXP xlam_adSEXP, SEXP xlambdabarSEXP, SEXP xVlambdaSEXP, SEXP xLSEXP, SEXP xalpha0SEXP, SEXP xbeta0SEXP, SEXP xinterceptsSEXP, SEXP xchainsSEXP, SEXP xkeepSEXP, SEXP xstartSEXP, SEXP xendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xY(xYSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type xlam_ad(xlam_adSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xlambdabar(xlambdabarSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xVlambda(xVlambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xL(xLSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xalpha0(xalpha0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xbeta0(xbeta0SEXP);
    Rcpp::traits::input_parameter< int >::type xintercepts(xinterceptsSEXP);
    Rcpp::traits::input_parameter< int >::type xchains(xchainsSEXP);
    Rcpp::traits::input_parameter< int >::type xkeep(xkeepSEXP);
    Rcpp::traits::input_parameter< int >::type xstart(xstartSEXP);
    Rcpp::traits::input_parameter< int >::type xend(xendSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bcfa_estim_spa(xY, xlam_ad, xlambdabar, xVlambda, xL, xalpha0, xbeta0, xintercepts, xchains, xkeep, xstart, xend));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_bcfa_simul_mu
Rcpp::List rcpp_bcfa_simul_mu(Rcpp::IntegerMatrix xordY, Rcpp::List xbounds, Rcpp::NumericMatrix xsigma, Rcpp::IntegerVector xplanwhere, Rcpp::IntegerMatrix xplanhow, int xnrep);
RcppExport SEXP _bcfa2s_rcpp_bcfa_simul_mu(SEXP xordYSEXP, SEXP xboundsSEXP, SEXP xsigmaSEXP, SEXP xplanwhereSEXP, SEXP xplanhowSEXP, SEXP xnrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type xordY(xordYSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type xbounds(xboundsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xsigma(xsigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type xplanwhere(xplanwhereSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type xplanhow(xplanhowSEXP);
    Rcpp::traits::input_parameter< int >::type xnrep(xnrepSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_bcfa_simul_mu(xordY, xbounds, xsigma, xplanwhere, xplanhow, xnrep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bcfa2s_rcpp_bcfa_estim_spa", (DL_FUNC) &_bcfa2s_rcpp_bcfa_estim_spa, 12},
    {"_bcfa2s_rcpp_bcfa_simul_mu", (DL_FUNC) &_bcfa2s_rcpp_bcfa_simul_mu, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcfa2s(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
