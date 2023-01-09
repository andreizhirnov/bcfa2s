#' Estimate polychoric correlations
#' @export

polychor_mle <- function(model, data) {
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
  ordk <- integer()
  for (j in mvs) {
    if (inherits(data[[j]], "ordered")) {
      d <- as.integer(data[[j]])
      ordY <- cbind(ordY,d)
      ordvar <- c(ordvar,j)
      ordk <- c(ordk, max(d, na.rm=TRUE)-2L)
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
  return(sigma)
}


