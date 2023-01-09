
## this is an adjusted function from bayesm
## converts dstar into the vector of thresholds
# dstartoc <- function(dstar) {
#   c(-100, dstar[1L], cumsum(exp(dstar[2L:length(dstar)])), 100)
# } 

dstartoc <- function(dstar) {
  c(-100, 0, cumsum(exp(dstar)), 100)
}

## this is an adjusted function from bayesm
## computes the log likelihood for the ordered probit without any covariates
lldstar0 <- function(dstar,y){
  y <- stats::na.omit(y)
  gamma <- dstartoc(dstar)
  arg <- stats::pnorm(gamma[y+2L])-stats::pnorm(gamma[y+1L])
  epsilon <- 1.0e-50
  arg <- ifelse(arg < epsilon, epsilon, arg)
  return(sum(log(arg)))
}

## computes the log likelihood for a polychoric correlation coefficient
llpoly <- function(rho, cutoffs, counts){
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  lpd <- numeric()
  for (i in seq_len(nrow(cutoffs))) { 
    lpd[i] <- log(
      mvtnorm::pmvnorm(lower = cutoffs[i,c(1L,3L)],
                       upper = cutoffs[i,c(2L,4L)],
                       corr=R))
  }
  sum(lpd*counts)
}

## computes the log likelihood for the estimations of structural parameters for CFA
negllspa <- function(theta, S, L, lam_ad, phi_ad) {
  lambda <- theta[1:nrow(lam_ad)]
  lvsd <- theta[nrow(lam_ad) + 1:ncol(L)]
  lvcor <- theta[nrow(lam_ad) + ncol(L) + 1:nrow(phi_ad)]
  var_Psi <- theta[nrow(lam_ad) + ncol(L) + nrow(phi_ad) + 1:nrow(L)]
    
  Lambda <- L
  for (i in seq_along(lambda)) {
    Lambda[lam_ad[i,1],lam_ad[i,2]] <- lambda[i]
   }
  phi <- diag(ncol(L))
  for (i in seq_along(lvcor)) {
    phi[phi_ad[i,1],phi_ad[i,2]] <- phi[phi_ad[i,2],phi_ad[i,1]] <- lvcor[i]
  }
  Psi <- diag(var_Psi) 
  Sigma <- Lambda %*% diag(lvsd) %*% phi %*% diag(lvsd) %*% t(Lambda) + Psi 
  return(determinant(Sigma, logarithm=TRUE)$modulus+sum(diag(solve(Sigma) %*%S)))
}

### mllspa(theta=qq$ini, S=qq$S, L=qq$L, lam_ad=qq$lam_ad, phi_ad=qq$phi_ad)
