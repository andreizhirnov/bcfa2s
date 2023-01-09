#' Estimate the structural parameters of the model
#' @param n.iter total number of iterations
#' @param n.thin thinning factor
#' @param n.burnin the number of iterations to burn
#' @return an \code{mcmc} object
#' @export

estim_spa <- function(y, 
                      L, lam_ad, lambda_mean, psi_mean, 
                      intercepts=TRUE,
                      n.iter=1000, n.thin=1, n.burnin=floor(n.iter/2), n.chains=3) {
## collect
  k <- list(Y=y, L=L, lam_ad=lam_ad-1L, lambdabar=lambda_mean, intercepts=as.numeric(intercepts),
            start=-n.burnin-1L, end=n.iter-n.burnin, keep=n.thin, chains = n.chains) 
  k$alpha0 <- rep(3, nrow(L))
  k$beta0 <- (k$alpha0-1)*psi_mean
  nmarkers <- sort(unique(k$lam_ad[,1L]))
  if (intercepts) {
    k$lambdabar <- c(k$lambdabar, rep(0, length(nmarkers)))
    k$L <- cbind(k$L,0)
    k$lam_ad <- rbind(k$lam_ad, cbind(nmarkers, ncol(L)))
  }
  k$Vlambda <- diag(length(k$lambdabar)) 
  names(k) <- paste0("x",names(k))
## send  
  w <- do.call("rcpp_bcfa_estim_spa", k)
## receive
  draw.nam <- character()
  for (i in 1:nrow(lam_ad)) {
    draw.nam <-c(draw.nam, paste0("lambda_",lam_ad[i,1L],"_",lam_ad[i,2L]))
  }
  if (intercepts) {
    for (j in nmarkers) {
      draw.nam <-c(draw.nam, paste0("mubar_",j))
    }
  }
  for (i in 1:ncol(L)) {
    draw.nam <-c(draw.nam, paste0("eta_",1:nrow(y),"_",i))
  }
  for (i in 1:ncol(L)) {
    draw.nam <-c(draw.nam, paste0("V_",i,"_",1:ncol(L)))
  }
  draw.nam <-c(draw.nam, paste0("psi_",1:nrow(L)))
####   draw.nam <-c(draw.nam, paste0("coef_",1:nrow(L)))  
  w <- lapply(w, function(u) {
    x <- cbind(u[["lambdadraw"]],u[["etadraw"]],u[["Vdraw"]],u[["psidraw"]])
####    x <- cbind(u[["lambdadraw"]],u[["etadraw"]],u[["Vdraw"]],u[["psidraw"]],u[["coefdraw"]])
    colnames(x) <- draw.nam
    coda::mcmc(x)
  })
  return(coda::mcmc.list(w))
}
