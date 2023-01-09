#' Estimate polychoric correlations
#' @export

estim_spa_mle <- function(model, data) {
  u <- trimws(strsplit(model[1], ";")[[1]])
  u <- u[grepl("=~", u)]
  v <- list() 
  for (g in seq_along(u)) {
    p <- trimws(strsplit(u[g], "=~")[[1]])
    v[[g]] <- unique(trimws(strsplit(p[2], "\\+")[[1]]))
    names(v)[g] <- trimws(p[1])
  }
  chmarkers <- sapply(v, "[",1)
  lvs <- names(v)
  mvs <- sort(unique(unlist(v)))
  ordY <- integer()
  ordvar <- character()
  for (j in mvs) {
    if (inherits(data[[j]], "ordered")) {
      d <- as.integer(data[[j]])
      ordY <- cbind(ordY,d)
      ordvar <- c(ordvar,j)
    }
  }
  thr <- apply(ordY, 2L, find_thresholds, simplify=FALSE)
  colnames(ordY) <- names(thr) <- NULL
  mvs <- c(ordvar, setdiff(mvs, ordvar))

  sigma <- diag(length(mvs))
  for (lefti in 2:length(mvs)) {
    for (righti in 1:(lefti-1)) {
      dt <- stats::na.omit(ordY[,c(lefti,righti)])
      agg <- stats::aggregate(list(N=rep(1,nrow(dt))), 
                       by=list(left=dt[,1], right=dt[,2]), 
                       FUN=length)
      co <- cbind(thr[[lefti]][agg$left],
                  thr[[lefti]][agg$left+1L],
                  thr[[righti]][agg$right],
                  thr[[righti]][agg$right+1L])
      e  <- stats::optimize(f=llpoly, cutoffs=co, counts=agg$N,
                     lower=-0.99, upper=0.99,
                     maximum=TRUE)
      sigma[lefti,righti] <- e$maximum
    }
  }
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  
  config <- do.call("cbind", lapply(v, function(g) as.numeric(mvs %in% g))) 
  lam_ad <- which(config>0, arr.ind=TRUE)
  markers <- integer()
  for (j in seq_along(lvs)) {
    t <- which(mvs==chmarkers[lvs[j]][1])
    markers <- c(markers, which(lam_ad[,"row"]==t & lam_ad[,"col"]==j))
  }
  lam_ad <- lam_ad[-markers,]
  lam_ad <- lam_ad[order(lam_ad[,"row"], lam_ad[,"col"]),]
  phi_ad <- do.call("rbind", lapply(2:length(lvs), function(u) cbind(rep(u, u-1L),1:(u-1L)))) 

### find MLE estimates for the structural parameters
  ini <- c(rep(1, nrow(lam_ad)+length(lvs)), rep(0, (length(lvs)*(length(lvs)-1)/2)), rep(1, length(mvs)))
  lb <- c(rep(-Inf, nrow(lam_ad)), rep(.Machine$double.eps, length(lvs)), rep(-0.99, nrow(phi_ad)), rep(.Machine$double.eps, length(mvs)))
  ub <- c(rep(Inf, nrow(lam_ad)), rep(Inf, length(lvs)), rep(0.99, nrow(phi_ad)), rep(Inf, length(mvs)))
  res <- stats::optim(par=ini,
               fn=negllspa, 
               S=sigma, L=config, lam_ad=lam_ad, phi_ad=phi_ad,
               lower=lb, upper=ub, method="L-BFGS-B")
  if (res$convergence>0) {
    stop("The algorithm did not converge")
    return(NULL)
  }
  
  lambda.hat <- res$par[1:nrow(lam_ad)]
  lvsd.hat <-  res$par[nrow(lam_ad) + seq_along(lvs)]
  lvcor <- res$par[nrow(lam_ad) + length(lvs) + 1:nrow(phi_ad)]
  lvcor.hat <-  diag(length(lvs))
  for (i in seq_along(lvcor)) {
    lvcor.hat[phi_ad[i,1],phi_ad[i,2]] <- lvcor.hat[phi_ad[i,2],phi_ad[i,1]] <- lvcor[i]
  }
  psi.hat <-  res$par[nrow(lam_ad) + length(lvs) + nrow(phi_ad) + seq_along(mvs)]
  
  return(
    list(
      mvs=mvs,
      lvs=lvs,
      L=config,
      lam_ad=lam_ad,
      phi_ad=phi_ad,
      S=sigma,
      thresholds=thr,
      lambda.hat = lambda.hat,
      lvsd.hat = lvsd.hat,
      lvcor.hat = lvcor.hat,
      psi.hat = psi.hat
    )
  )
}


