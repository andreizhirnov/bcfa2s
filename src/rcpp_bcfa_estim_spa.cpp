#include "RcppArmadillo.h"
 
// [[Rcpp::depends(RcppArmadillo)]]

// draw the inverse of the variance-covariance of residuals
arma::mat azs_P(arma::vec const& alpha0, arma::vec const& b, int nmv, int nobs) {
  // Purpose: draw loadings
  // use calculated g and G (see Sik-Yum Lee's derivation)

  arma::mat P(nmv, nmv);
  P.zeros();
  for (int j=0; j<nmv; j++ ) { 
    P(j,j) = R::rgamma(1.0*nobs/2+alpha0(j), 1/b(j));
  }
  return (P);
}
 
// draw lambda
arma::vec azs_lambda(arma::mat const& G, arma::vec const& g) {
  // Purpose: draw loadings
  // use calculated g and G (see Sik-Yum Lee's derivation)
  return (arma::mvnrnd(g, G));
}

// draw eta
arma::mat azs_eta(arma::mat const& L, int const& nobs, int const& nlv,
                  arma::mat const& Y, arma::mat const& B, arma::mat const& P,
                  arma::vec const& mubar) {

  // Purpose: draw etas

  // Arguments: lambda, mu, inverse of correlation matrix (B)

  // with L: rectangular matrix with loadings

  // eta ~ MVN(mean=(L'PL + B)^(-1)(L'P Y), precision = L'PL + B)
 
  arma::mat M = L.head_cols(B.n_rows);
  arma::mat vcov = arma::inv_sympd(trans(M)*P*M + B);
  arma::mat meanpref = vcov*M.t()*P;
  arma::mat eta(nobs, nlv);
  for (int i=0; i < nobs; i++) {
     eta.row(i) = arma::mvnrnd(meanpref*(Y.row(i).t()-mubar), vcov).t();
  }
  return (eta);
}

// draw inverse of correlation matrix B
arma::mat azs_B(arma::mat const& eta, double const& k) {
  
  // with psi = (trans(eta)*eta)^{-1} and nu = N/2 + p + 3-2k
  // B|data and else ~ Wishart( psi, nu)

  int p = eta.n_cols;
  int N = eta.n_rows;

  arma::mat psi = arma::inv_sympd(eta.t()*eta);
  return(arma::wishrnd(psi, N + p + 3-2*k));
}

// [[Rcpp::export]]
Rcpp::List rcpp_bcfa_estim_spa(Rcpp::NumericMatrix xY,                              
                               Rcpp::IntegerMatrix xlam_ad, 
                               Rcpp::NumericVector xlambdabar,
                               Rcpp::NumericMatrix xVlambda,
                               Rcpp::NumericMatrix xL,
                               Rcpp::NumericVector xalpha0,                          
                               Rcpp::NumericVector xbeta0,
                               int xintercepts,
                               int xchains,
                               int xkeep,
                               int xstart,
                               int xend){

// useful constants; convert Rcpp objects into arma
  arma::Mat<int> lam_ad = Rcpp::as<arma::Mat<int>>(xlam_ad);
  arma::mat Y = Rcpp::as<arma::mat>(xY);
  arma::vec lambdabar = Rcpp::as<arma::vec>(xlambdabar);
  arma::mat L = Rcpp::as<arma::mat>(xL);
  arma::mat Vlambda = Rcpp::as<arma::mat>(xVlambda); 
  arma::mat A0 = arma::inv_sympd(Vlambda);
  
  arma::vec alpha0 = Rcpp::as<arma::vec>(xalpha0);
  arma::vec beta0 = Rcpp::as<arma::vec>(xbeta0);
  
  const bool intercepts = !!xintercepts;
  const int nmv = Y.n_cols; // number of manifest variables 
  const int nx = L.n_cols; // number of latent variables in the matrix
  const int nlv = nx - xintercepts; // number of latent variables
  const int nobs = Y.n_rows; // number of observations
  const int nlam = lam_ad.n_rows; // number of loadings
  const int keep = round(xkeep);
  const int start = round(xstart);
  const int end = round(xend);
  const int chains = round(xchains);
  
  arma::vec Yvar(nmv);
  for (int j=0; j<nmv; j++) {
    Yvar(j) = arma::dot(Y.col(j), Y.col(j));
  }
 
  arma::mat V;
  arma::mat B;
  arma::mat eta(nobs, nlv);
  arma::vec lambda(nlam);
  arma::mat P(nmv,nmv);
  arma::vec g(nlam);
  arma::mat G(nlam,nlam);
  arma::vec b(nmv);
  arma::vec gj(nlv);
  arma::mat Gj(nlv,nlv);
  arma::mat iGj(nlv,nlv);
  arma::mat A0j[nmv];
   
  arma::mat xeta(nobs, nx);
  if (intercepts) {
    xeta.tail_cols(1).ones(); 
    }
  
  arma::vec mubar(nmv);
    mubar.zeros();
    
  // conversion matrix for the loadings
  arma::cube lam2Lj(nx, nlam, nmv);
  lam2Lj.zeros();

  for (int j=0; j<nlam; j++){
    lam2Lj(lam_ad(j,1),j,lam_ad(j,0)) = 1;
    }

  for (int j=0; j<nmv; j++){
      A0j[j] = lam2Lj.slice(j)*A0*lam2Lj.slice(j).t();
    }

  Rcpp::List ret(chains);
//// objects for export
  int mkeep;
  int mkeepmax = floor(end/keep);
  arma::mat lambdadraw(mkeepmax, nlam);
  arma::mat etadraw(mkeepmax, nlv*nobs);
  arma::mat Vdraw(mkeepmax, nlv*nlv);
  arma::mat psidraw(mkeepmax, nmv);

  int tpct = ceil((end-start+1)/10);
  int ready;

/// separate by chain

for (int h=0; h<chains; h++) {
  Rcpp::Rcout << "Running chain " << h+1 << "\n";
  // initialize objects
  lambda = arma::mvnrnd(lambdabar, Vlambda);
  for (size_t i=0; i<lam_ad.n_rows; i++) {
    L(lam_ad(i,0), lam_ad(i,1)) = lambda(i); 
  }
  if (intercepts) {
    mubar = L.tail_cols(1);
  }
  
  P = azs_P(alpha0, beta0, nmv, nobs);

  eta = (Y - arma::repmat(mubar.t(), nobs, 1))*P*L.head_cols(nlv)*arma::inv_sympd(L.head_cols(nlv).t()*P*L.head_cols(nlv)); 
  xeta.head_cols(nlv) = eta;
  B = azs_B(eta, 2);

/// perform iteration
  for (int rep = start; rep < end; rep++) {

  eta = azs_eta(L, nobs, nlv, Y, B, P, mubar);
  xeta.head_cols(nlv) = eta;
  G.zeros();
  g.zeros();

for (int j=0; j<nmv; j++) {
  iGj = A0j[j] + xeta.t()*xeta;
  Gj = arma::inv_sympd(iGj);
  G += lam2Lj.slice(j).t()*Gj*lam2Lj.slice(j)/P(j,j);
  gj = Gj*(A0j[j]*lam2Lj.slice(j)*lambdabar + xeta.t()*Y.col(j));
  g += lam2Lj.slice(j).t()*gj;
  b(j) = beta0(j) + 0.5*arma::as_scalar(Yvar(j) - gj.t()*iGj*gj + trans(lam2Lj.slice(j)*lambdabar)*A0j[j]*(lam2Lj.slice(j)*lambdabar));
}
  P = azs_P(alpha0, b, nmv, nobs);
  lambda = azs_lambda(G, g);
  for (size_t i=0; i<lam_ad.n_rows; i++) {
      L(lam_ad(i,0), lam_ad(i,1)) = lambda(i);
    }
  if (intercepts) {
      mubar = L.tail_cols(1);
    }
  B = azs_B(eta, 2);

    /// parameter 2 here means that the covariance matrices with a positive ridge are
    //   slightly favored over those with high off-diagonal elements

////////////////////////////////////// keeper
  if (rep>=0 && rep%keep==0) {
    mkeep = rep/keep;
    V = arma::inv_sympd(B);
    lambdadraw.row(mkeep) = lambda.t();
    psidraw.row(mkeep) = trans(1/P.diag());
    for (int j=0; j< nlv; j++) {
      etadraw.submat(mkeep, nobs*j, mkeep, nobs*(j+1)-1) = eta.col(j).t();
      Vdraw.submat(mkeep, nlv*j, mkeep, nlv*(j+1)-1) = V.row(j);
    }
  }
  if ((rep-start)%tpct==0) {
    ready = floor(10*(rep-start)/tpct);
    Rcpp::Rcout << "..." << ready << "% complete\n";
  }
}
  ret[h] = Rcpp::List::create(
    Rcpp::Named("lambdadraw") = lambdadraw,
    Rcpp::Named("etadraw") = etadraw,
    Rcpp::Named("Vdraw") = Vdraw,
    Rcpp::Named("psidraw") = psidraw
    );
}
  return (ret);
}


