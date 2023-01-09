#include "RcppArmadillo.h"
 
// [[Rcpp::depends(RcppArmadillo)]]
double azs_tnorm(double const& mean, double const& sd, double const& left, double const& right){
  
  // Purpose: draw from a truncated normal distribution 
  double uni = arma::randu<double>();
  double lprob = R::pnorm(left, mean, sd, 1,0);
  double rprob = R::pnorm(right, mean, sd, 1,0);
  return(R::qnorm(lprob + (rprob-lprob)*uni, mean, sd, 1, 0)); 
}

// [[Rcpp::export]]
Rcpp::List rcpp_bcfa_simul_mu(Rcpp::IntegerMatrix xordY,
                               Rcpp::List xbounds,
                               Rcpp::NumericMatrix xsigma,
                               Rcpp::IntegerVector xplanwhere,
                               Rcpp::IntegerMatrix xplanhow,
                               int xnrep){

// useful constants; convert Rcpp objects into arma
  arma::Mat<int> ordY = Rcpp::as<arma::Mat<int>>(xordY);
  arma::mat sigma = Rcpp::as<arma::mat>(xsigma);
  const int nmv = ordY.n_cols; // number of manifest variables 
  const int nobs = ordY.n_rows; // number of observations  
  const int nrep = round(xnrep);
  
  arma::mat bounds[nmv];
  for (int j=0; j<nmv; j++) {
    bounds[j] = Rcpp::as<arma::mat>(xbounds[j]);
  } 
  arma::Mat<int> planhow = Rcpp::as<arma::Mat<int>>(xplanhow);
  arma::ivec whichsmethod = Rcpp::as<arma::ivec>(xplanwhere);  
  
// prepare simulation methods
  const int nsmethod = planhow.n_rows;
  arma::ivec smethodstop = planhow.col(0);
  
  arma::Mat<int> smethodseq = planhow.submat(0, 1, nsmethod-1, nmv);
  arma::mat smethodsd(nsmethod, nmv);
  arma::vec smethodmean[nsmethod][nmv];
  arma::mat smethodperm[nsmethod];
  arma::mat smethodiperm[nsmethod];

// permutation matrices
  arma::mat tempperm(nmv,nmv);
  arma::mat tempvcperm(nmv,nmv);
  arma::mat tempcov;
  arma::mat tempprec; 
  
  for (int m=0; m<nsmethod; m++) {
    tempperm.zeros();
    for (int j=0; j< nmv; j++) {
      for (int u=0; u<nmv; u++) {
        if (smethodseq(m,u)==j) {
          tempperm(u,j)=1;
          }
        }
      }
    smethodperm[m]=tempperm;
    smethodiperm[m]=arma::inv(tempperm);
    tempvcperm = tempperm*sigma*tempperm.t(); 
    
    for (int j=1; j<nmv; j++) {
      tempcov = tempvcperm.submat(j, 0, j, j-1).t();
      tempprec = arma::inv(tempvcperm.submat(0, 0, j-1, j-1));
      smethodmean[m][j] = tempprec*tempcov;
      smethodsd(m, j) = std::sqrt(1.0 - arma::as_scalar(tempcov.t()*tempprec*tempcov)); 
    }
  } 
  
//  initialize mu
    arma::cube mu(nobs, nmv, nrep);

// loop over rows
  arma::vec muhat(nmv);
  arma::ivec tseq(nmv);

  int m;
  double mean;
  double sd;

  for (int i=0; i<nobs; i++) {
 //   find out which method to use
     m = whichsmethod(i);
     tseq = smethodseq.row(m).t();
     for (int k=0; k<nrep; k++) {
      muhat.zeros();
 //   draw the first entry,
 //   apply the thresholds that correspond to the response in the first entry
      muhat(0) = azs_tnorm(0, 1, bounds[tseq(0)](ordY(i,tseq(0)),0), bounds[tseq(0)](ordY(i,tseq(0)),1));
 
//   fill other entries before the stop
     for (int j=1; j<nmv; j++) {
       mean = arma::dot(smethodmean[m][j], muhat(arma::span(0,j-1)));
       sd = smethodsd(m, j); 
       if (ordY(i,tseq(j))==-1) { 
        muhat(j) = R::rnorm(mean, sd);
        } else {
          muhat(j) = azs_tnorm(mean, sd, bounds[tseq(j)](ordY(i,tseq(j)),0), bounds[tseq(j)](ordY(i,tseq(j)),1));
          }
     }
  //   permute back and insert into the mu table
      mu.slice(k).row(i) = trans(smethodiperm[m]*muhat);
      }
  }

Rcpp::List muli(nrep);
for (int k=0; k<nrep; k++) {
  muli[k] = mu.slice(k);
  }
  return(muli);
}
