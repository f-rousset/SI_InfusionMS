regAbcrf.formula2 <- function (
    formula, data, ntree = 500, mtry = max(floor((dim(data)[2] - 
                                                    1)/3), 1), sampsize = min(1e+05, nrow(data)), 
    paral = FALSE,
    ncores = if (paral) max(mydetectCores() - 1L, 1) else 1,  # when run on a shared cluster, it is important to control this explicitly
    # Otherwise ranger's 'Default is number of CPUs available.' 
    importance = "impurity", ...) 
{
  if (!inherits(formula, "formula")) 
    stop("regAbcrf.formula is only for formula objects")
  if (!inherits(data, "data.frame")) 
    stop("data needs to be a data.frame object")
  if (nrow(data) == 0L || is.null(nrow(data))) 
    stop("no simulation in the reference table (resp, sumstat)")
  if (sampsize > nrow(data)) 
    stop("sampsize too large")
  if ((!is.logical(paral)) || (length(paral) != 1L)) 
    stop("paral should be TRUE or FALSE")
  if (is.na(ncores)) {
    warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
    ncores <- 1
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf))
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (!is.numeric(model.response(mf))) 
    stop("response variable should be numeric")
  model.rf <- ranger::ranger(formula, data = data, num.trees = ntree, 
                     mtry = mtry, sample.fraction = sampsize/nrow(data), 
                     num.threads = ncores, # when run on a shared cluster, it is important to control this explicitly
                     # Otherwise ranger's 'Default is number of CPUs available.' 
                     keep.inbag = TRUE, importance = importance, ...)
  model.rf$NMAE <- mean(abs((model.response(mf) - model.rf$predictions)/model.response(mf)))
  cl <- match.call()
  cl[[1]] <- as.name("regAbcrf")
  x <- list(call = cl, formula = formula, model.rf = model.rf)
  class(x) <- "regAbcrf"
  x
}


build_abcrf_projections <- function(
    reftable_abcrf, parNames, Sobs, ntree, min.node.size, 
    ncores, # when run on a shared cluster, it is important to control this explicitly
    # Otherwise ranger's 'Default is number of CPUs available.
    importance="none") {
  resu <- vector("list", length(parNames))
  names(resu) <- parNames
  for (parm in parNames) {
    cat(crayon::yellow(parm," "))
    resu[[parm]] <- regAbcrf.formula2(as.formula(paste(parm,"~.")),
                             reftable_abcrf[,c(parm, names(Sobs))], ntree=ntree, 
                             min.node.size = min.node.size, paral = TRUE, ncores = ncores,
                             importance=importance)
  } 
  resu
}

.findweights <- get("findweights", asNamespace("abcrf"))

HPDint <- function (object, obs, training, level=0.95, densityEst="stats::density", ncores=1L, ...) {
  if (!inherits(object, "regAbcrf")) 
    stop("object not of class regAbcrf")
  if (!inherits(training, "data.frame")) 
    stop("training needs to be a data.frame object")
  if (!inherits(obs, "data.frame")) 
    stop("obs needs to be a data.frame object")
  if (nrow(obs) == 0L || is.null(nrow(obs))) 
    stop("no data in obs")
  if (nrow(training) == 0L || is.null(nrow(training))) 
    stop("no simulation in the training reference table (response, sumstat)")
  x <- obs
  if (!is.null(x)) {
    if (is.vector(x)) {
      x <- matrix(x, ncol = 1)
    }
    if (nrow(x) == 0) 
      stop("obs has 0 rows")
    if (any(is.na(x))) 
      stop("missing values in obs")
  }
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  mf$data <- training
  training <- mf$data
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  resp <- model.response(mf)
  obj <- object$model.rf
  inbag <- matrix(unlist(obj$inbag.counts, use.names = FALSE), 
                  ncol = obj$num.trees, byrow = FALSE)
  obj[["origNodes"]] <- predict(obj, training, predict.all = TRUE, 
                                num.threads = ncores)$predictions
  obj[["origObs"]] <- model.response(mf)
  origObs <- obj$origObs
  origNodes <- obj$origNodes
  nodes <- predict(obj, x, predict.all = TRUE, num.threads = ncores)$predictions
  if (is.null(dim(nodes))) 
    nodes <- matrix(nodes, nrow = 1)
  ntree <- obj$num.trees
  nobs <- object$model.rf$num.samples
  nnew <- nrow(x)
  if (nnew>1L) stop("density estimation for multiple obs not yet implemented") 
  weights <- .findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree))
  weights.std <- weights/ntree
  if (densityEst=="stats::density") {
    postDensity <- density(resp, weights = weights.std[, 1], warnWbw=FALSE,  ...)
    dens <- postDensity$y
    xin <- postDensity$x
  } else if (densityEst=="locfit") {
    requireNamespace("locfit",quietly = TRUE) # lp()
    postDensity <- locfit::locfit( ~ locfit::lp(resp), weights=weights.std[, 1], ...)
    xrange <- range(resp)
    xin <- seq(xrange[1], xrange[2], length.out=10000L)
    dens <- predict(postDensity, newdata = xin)
  } else if (densityEst=="hdrcde") {
    requireNamespace("hdrcde",quietly = TRUE) # lp()
    ##postDensity <- hdrcde::hdr.2d(x = pod.weighted.sample.rf[[1]][[i.pod]]$Ne, y = pod.weighted.sample.rf[[2]][[i.pod]]$µ, 
    ##                                            prob = c(0.95)) 
    xrange <- range(resp)
    xin <- seq(xrange[1], xrange[2], length.out=10000L)
    dens <- predict(postDensity, newdata = xin)
  }
  ord <- order(dens, decreasing = TRUE)
  cumw <- cumsum(dens[ord])
  cumw <- cumw/cumw[length(cumw)]
  mwc <- which.max(cumw>level) # ~ min(which(cumw>level)) as which.max returns the first value. Measurably faster 
  #(more that in comparisons here :https://stackoverflow.com/questions/29388334/find-position-of-first-value-greater-than-x-in-a-vector)
  range(xin[ord[1:mwc]])
}

wrap_predict.regAbcrf <- function(object, parm, obs, reftable, ret_HPDint, level=0.95, 
                                  ncores, # names in predict.regAbcrf
                                  ...) {
  resu_abcrf <- predict(object,  obs=obs, training=reftable[,c(parm, names(obs))],...)
  if (ret_HPDint) {
    if (FALSE) class(resu_abcrf) <- "hack_FR" # the returned resu_abcrf object does NOT inherit from class "list"
    # # its class has specific methods that interfere with str(), 
    # # even causing a form of error when elements are added to this object (as below)
    # # changing the class avoids the error. (But the error is non-stopping so we keep this as 'comment' only)
    resu_abcrf$HPDint_d <- HPDint(object,  obs=obs, training=reftable[,c(parm, names(obs))], level=level,
                                  ncores=ncores, 
                                  densityEst="stats::density")
    if (ret_HPDint>1.5) { # this tends to generate errors
      resu_abcrf$HPDint_l <- HPDint(object,  obs=obs, training=reftable[,c(parm, names(obs))], level=level,
                                    ncores=ncores,
                                    densityEst="locfit")
    }
  }
  resu_abcrf
} 

abcrf_replicate <- function(Sobs, abcrf_proj, reftable_abcrf,
                            parNames=names(abcrf_proj),
                            plot.=interactive(),
                            ncores_abcrf=1L,
                            ret_HPDint=TRUE
                            ) { 
  Sobs.df <- as.data.frame(t(Sobs))
  resu <- vector("list", length(parNames))
  names(resu) <- parNames
  for (parm in parNames) {
    resu[[parm]] <- wrap_predict.regAbcrf(abcrf_proj[[parm]], parm=parm, obs=Sobs.df, reftable=reftable_abcrf,
                                          ret_HPDint=ret_HPDint, ncores=ncores_abcrf, paral=(ncores_abcrf>1L))
    if (plot.) abcrf::densityPlot(abcrf_proj[[parm]], obs=Sobs.df, training=reftable_abcrf, warnWbw=FALSE,
                                  main=paste("Posterior density for",parm))
  } 
  resu
}

