# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_bcfa_estim_spa <- function(xY, xlam_ad, xlambdabar, xVlambda, xL, xalpha0, xbeta0, xintercepts, xchains, xkeep, xstart, xend) {
    .Call(`_bcfa2s_rcpp_bcfa_estim_spa`, xY, xlam_ad, xlambdabar, xVlambda, xL, xalpha0, xbeta0, xintercepts, xchains, xkeep, xstart, xend)
}

rcpp_bcfa_simul_mu <- function(xordY, xbounds, xsigma, xplanwhere, xplanhow, xnrep) {
    .Call(`_bcfa2s_rcpp_bcfa_simul_mu`, xordY, xbounds, xsigma, xplanwhere, xplanhow, xnrep)
}

