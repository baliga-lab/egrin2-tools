#
# resample.R - Resampling utility functions
#

cvar <- function(genes, conditions, ratios, mode=c("median", "mean", "none")[1]) {

  m <- abs(apply(ratios[genes, conditions, drop=F], 2, mean, na.rm=T))
  m[which(m == 0)] = 1e-6  # replace the 0 places with 1e-6

  if (mode=="mean") { 
    var.m <- mean(apply(ratios[genes, conditions, drop=F], 2, sd, na.rm=T) / m, na.rm=T)
  } else if (mode=="median"){
    var.m <- median(apply(ratios[genes, conditions, drop=F], 2, sd, na.rm=T) / m, na.rm=T)
  } else {
    var.m <- apply(ratios[genes, conditions, drop=F], 2, sd, na.rm=T) / m
  }
  return(var.m)
}

make.to.r.item <- function(ratios, geneSetSize, resamples, mode, genePool) {

  message(paste(c("gene set size=", geneSetSize), collapse=""))
  i <- do.call(rbind, lapply(seq(1:resamples), function(j) {
    # returns the first to.r that is non completely NA
    repeat {
      to.r <- cvar(genes=sample(genePool, geneSetSize),
                   conditions=colnames(ratios),
                   ratios=ratios,
                   mode=mode)
      if (sum(is.na(to.r)) != length(to.r)) break
    }
    return(to.r)
  }))

  i.2 <- do.call(cbind, lapply(seq(1:dim(i)[2]), function(j) {          
    j <- i[,j]
    j.ecdf <- ecdf(j)  # ecdf is R function

    # make everything above lowest 7.5% quantile = 0
    j[which(j > quantile(j.ecdf, probs=seq(0, 1, .075))[2])] = 0
    return(j)
  }))

  colnames(i.2) <- colnames(i)
  i.2 <- as(i.2, "sparseMatrix")
}

resampleRandomConditions <- function(ratios, geneSetSizes=seq(3, 200, 1),
                                     resamples=20000, mode="none") {
  require('multicore')
  genePool = rownames(ratios)
  to.r <- mclapply(seq(1:length(geneSetSizes)), function(i) {
    make.to.r.item(ratios, geneSetSizes[i], resamples, mode, genePool)
  })
  names(to.r) <- as.character(geneSetSizes)
  invisible(to.r)
}
