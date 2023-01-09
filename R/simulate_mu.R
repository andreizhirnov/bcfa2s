#' Create replicates of the dataset with continuous 'latent' variables consistent with the ordinal data

simulate_mu <- function(data, mvs, thresholds, sigma, nrep=20L) {
  ordY <- integer() 
  for (j in mvs) {
    if (inherits(data[[j]], "ordered")) {
      d <- as.integer(data[[j]])
      ordY <- cbind(ordY,d)
    }
  } 
  colnames(ordY) <- NULL 
  bounds <- lapply(thresholds, function(u) cbind(u[1:(length(u)-1)], u[2:length(u)]))
  wtr <- apply(ordY, 1L, function(u) {
    v <- which(!is.na(u))
    return(c(length(v),v-1L,setdiff(seq_along(mvs), v)-1L))
    }, simplify=FALSE)
  uwtr <- unique(wtr)
  planbase <- do.call("rbind", wtr)
  planhow <- do.call("rbind", uwtr)
  
  planwhere <- rep(NA, nrow(planbase))
  for (j in 1:nrow(planhow)) { 
    planwhere[which(apply(planbase, 1L, function(x) return(all(x == planhow[j,]))))] <- j
  }
  k <- list(
    ordY=replace(ordY-1L, is.na(ordY), -1L),
    bounds=bounds,
    sigma=sigma,
    nrep=nrep,
    planhow=as.matrix(planhow),
    planwhere=planwhere-1L
    )
  names(k) <- paste0("x",names(k))
  w <- do.call("rcpp_bcfa_simul_mu", k)
  return(w)
}