## Minor variations and comments for likelihoodAsy::rstar()
my_rstar <- function(data, thetainit, floglik, fscore = NULL, fpsi, psival, 
                     datagen, R = 1000, seed = NULL, trace = TRUE, ronly = FALSE, 
                     psidesc = NULL, constr.opt = "alabama") # constr.opt ignored
{
  if (!is.list(data)) {
    warning("data should be provided as a list\n")
    data <- as.list(data)
  }
  if (!is.numeric(thetainit)) 
    stop("a starting point for the parameter theta is required \n")
  if (!is.numeric(psival)) 
    stop("a value for the parameter of interest is required \n")
  f0 <- floglik(thetainit, data)
  if (!is.numeric(f0)) 
    stop("problems in the loglikelihood function \n")
  if (!is.null(fscore)) {
    g0 <- fscore(thetainit, data)
    if (length(g0) != length(thetainit)) 
      stop("size of starting point different from the size of score function \n")
  }
  if (!ronly) {
    data0 <- datagen(thetainit, data)
    f0 <- floglik(thetainit, data0)
    if (!is.numeric(f0)) 
      stop("problems in the function to simulate data \n")
  }
  if ((constr.opt != "solnp") & (constr.opt != "alabama")) 
    stop("constrained optimizer must be either 'solnp' or 'alabama'")
  p <- length(thetainit)
  if (trace) 
    cat("get mle ....", "\t")
  min.floglik <- function(theta, data) (-1) * floglik(theta, data)
  min.fscore <- if (is.null(fscore)) {
    NULL
  } else function(theta, data) (-1) * fscore(theta, data)
  if (FALSE) { # rstar original code
    obj.hat <- nlminb(thetainit, min.floglik, min.fscore, data = data)
    theta.hat <- obj.hat$par
    el.hat <- floglik(theta.hat, data)
  } else {
    newMSL <- summLik(data$slik,parm = NULL, data=data$newobs)
    theta.hat <- attr(newMSL,"profpt")
    el.hat <- newMSL[[1]]
  }
  if (!ronly) {
    j.hat <- if (is.null(fscore)) 
      -pracma::hessian(floglik, theta.hat, data = data)
    else -pracma::jacobian(fscore, theta.hat, data = data)
    score.hat.data <- if (is.null(fscore)) 
      pracma::grad(floglik, theta.hat, data = data)
    else fscore(theta.hat, data = data)
  }
  if (length(thetainit) > 1) {
    if (!ronly) {
      var.theta.hat <- try(solve(j.hat))
      if (inherits(var.theta.hat,"try-error")) {
        j.hat <- spaMM::regularize(j.hat) # quick fix for inherently defiective procedure
        var.theta.hat <- solve(j.hat)
      }
      se.theta.hat <- sqrt(diag(var.theta.hat))
    }
    if (trace) 
      cat("get mle under the null....", "\n") # the profile lik for h0, for the standard LRT anyway
    if (FALSE) {
      psifcn.mod <- if (constr.opt == "solnp") 
        function(theta, data) fpsi(theta)
      else function(theta, data) fpsi(theta) - psival
      objHyp <- if (constr.opt == "solnp") 
        solnp(theta.hat, fun = min.floglik, eqfun = psifcn.mod, 
              eqB = psival, control = list(trace = 0), data = data)
      else constrOptim.nl(theta.hat, fn = min.floglik, heq = psifcn.mod, 
                          gr = min.fscore, control.outer = list(trace = FALSE), 
                          data = data)
      theta.til <- objHyp$par
      el.til <- floglik(theta.til, data)
    } else {
      newMSL <- summLik(data$slik,parm = psival, data=data$newobs)
      theta.til <- attr(newMSL,"profpt")
      el.til <- newMSL[[1]]
    }
    psi.hat <- fpsi(theta.hat)
    if (!ronly) {
      j.til <- if (is.null(fscore)) 
        -pracma::hessian(floglik, theta.til, data = data)
      else -pracma::jacobian(fscore, theta.til, data = data)
      score.til.data <- if (is.null(fscore)) 
        pracma::grad(floglik, theta.til, data = data)
      else fscore(theta.til, data)
      dpsi.dtheta <- pracma::grad(fpsi, theta.hat)
      var.psi.hat <- dpsi.dtheta %*% var.theta.hat %*% 
        dpsi.dtheta
      se.psi.hat <- sqrt(var.psi.hat)
    }
    dLR <- el.hat - el.til
    if (dLR<0) {
      newmax <- MSL(data$slik, init=theta.til, eval_RMSEs=FALSE, CIs=FALSE)
      theta.hat <- newmax$MSL$MSLE
      el.hat <- newmax$MSL$maxlogL
      dLR <- el.hat - el.til
      if (dLR<0) dLR <- 0 # (___F I X M E___ quick patch but minor issue relative to broader ones with rstar)
    }
    r <- sqrt(2 * (dLR)) * sign(psi.hat - psival)
    if (!ronly) {
      C.hat <- pracma::grad(fpsi, theta.hat)
      C.til <- pracma::grad(fpsi, theta.til)
      k <- which(C.hat != 0)[1]
      obj.info <- likelihoodAsy:::.newinfo(p, k, C.hat, C.til, j.hat, j.til, 
                                           score.hat.data, score.til.data, theta.hat, theta.til, 
                                           fpsi)
      j.hat.new <- obj.info$j.hat.new
      j.til.new <- obj.info$j.til.new
    }
    if (ronly) 
      out <- list(r = r, theta.hat = theta.hat, psi.hat = psi.hat, 
                  theta.hyp = theta.til, psi.hyp = psival)
    if (!ronly) {
      if (trace) 
        cat("start Monte Carlo computation", "\n")
      meanAll <- rep(0, 2 * p + 1)
      prodAll <- matrix(0, 2 * p + 1, 2 * p + 1)
      dataSim <- as.list(data, all.names = TRUE)
      seed.in <- if (is.null(seed)) 
        sample.int(2^30, 1)
      else seed
      set.seed(seed.in)
      if (trace) 
        pb <- txtProgressBar(style = 3)
      for (i in 1:R) {
        if ((i%%(R/10) == 0) & trace) 
          setTxtProgressBar(pb, i/R)
        dataSim <- datagen(theta.hat, data = data)
        l1 <- floglik(theta.hat, dataSim) # for full theta vectors, hence no profiling
        l0 <- floglik(theta.til, dataSim)
        score.hat <- if (is.null(fscore)) 
          pracma::grad(f = floglik, x0 = theta.hat, data = dataSim) # or numDeriv::grad(f = floglik, x = theta.hat, data = dataSim)
        else fscore(theta.hat, dataSim)
        score.til <- if (is.null(fscore)) 
          pracma::grad(f = floglik, x0 = theta.til, data = dataSim) # or numDeriv::grad(f = floglik, x = theta.til, data = dataSim)
        else fscore(theta.til, dataSim)
        obj.score <- likelihoodAsy:::.newscores(p, k, C.hat, C.til, score.hat, 
                                                score.til)
        uhh <- c(obj.score$score.new.hat, obj.score$score.new.til, 
                 (l1 - l0))
        meanAll <- meanAll + uhh/R
        prodAll <- prodAll + tcrossprod(uhh)/R
      }
      if (trace) 
        close(pb)
      covAll <- prodAll * R/(R - 1) - tcrossprod(meanAll) * 
        R/(R - 1)
      S <- covAll[1:p, (p + 1):(2 * p)]
      i.hat <- covAll[1:p, 1:p]
      i.hatInv <- qr.solve(i.hat, tol = 10^-20)
      q <- covAll[2 * p + 1, 1:p]
      SS2 <- t(S) %*% i.hatInv %*% j.hat.new
      SS1 <- q %*% i.hatInv %*% j.hat.new
      indpsi <- k
      locSS2 <- t(SS2[-indpsi, -indpsi])
      locinv <- try(solve(locSS2))
      if (inherits(locinv,"try-error")) {
        locSS2 <- spaMM::regularize(locSS2) # quick fix for inherently defiective procedure
        locinv <- solve(locSS2)
      }
      
      numU <- SS1[indpsi] - SS2[-indpsi, indpsi] %*% locinv %*% SS1[-indpsi]
      j.hatInv.new <- try(solve(j.hat.new))
      if (inherits(j.hatInv.new,"try-error")) {
        j.hatInv.new <- spaMM::regularize(j.hatInv.new) # quick fix for inherently defiective procedure
        j.hatInv.new <- solve(j.hat.new)
      }
      jProf.new <- 1/j.hatInv.new[indpsi, indpsi]
      u <- numU/sqrt(jProf.new)
      CPsi <- det(as.matrix(SS2[-indpsi, -indpsi]))/sqrt(det(as.matrix(j.til.new[-indpsi, 
                                                                                 -indpsi])) * det(as.matrix(j.hat.new[-indpsi, 
                                                                                                                      -indpsi])))
      # cf eq 12 in PierceB17:
      NP <- (1/r) * log(CPsi)
      INF <- (1/r) * log(u/r)
      rs <- r + NP + INF
      
      out <- list(r = r, NP = drop(NP), INF = drop(INF), 
                  rs = drop(rs), theta.hat = theta.hat, info.hat = j.hat, 
                  se.theta.hat = se.theta.hat, psi.hat = psi.hat, 
                  se.psi.hat = drop(se.psi.hat), theta.hyp = theta.til, 
                  psi.hyp = psival, seed = seed.in,
                  utot=CPsi*u/r) # added 'utot' for experiments
    }
  }
  if (length(thetainit) == 1) { # no nuisance param -- I don't really use this case
    if (!ronly) {
      j.hat <- as.numeric(j.hat)
      var.theta.hat <- 1/(j.hat)
      se.theta.hat <- sqrt(var.theta.hat)
    }
    psifcn.mod <- function(theta, data) fpsi(theta) - psival
    objHyp <- nleqslv::nleqslv(theta.hat, psifcn.mod, data = data)
    theta.til <- objHyp$x
    el.til <- floglik(theta.til, data)
    psi.hat <- fpsi(theta.hat)
    if (!ronly) {
      j.til <- if (is.null(fscore)) 
        -pracma::hessian(floglik, theta.til, data = data)
      else -pracma::jacobian(fscore, theta.til, data = data)
      score.til.data <- if (is.null(fscore)) 
        pracma::grad(floglik, theta.til, data = data)
      else fscore(theta.til, data)
      dpsi.dtheta <- pracma::grad(fpsi, theta.hat)
      var.psi.hat <- dpsi.dtheta^2 * var.theta.hat
      se.psi.hat <- sqrt(var.psi.hat)
      j.hat.new <- j.hat/dpsi.dtheta^2
    }
    r <- sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
    if (ronly) 
      out <- list(r = r, theta.hat = theta.hat, psi.hat = psi.hat, 
                  theta.hyp = theta.til, psi.hyp = psival)
    if (!ronly) {
      if (trace) 
        cat("start Monte Carlo computation", "\n")
      meanAll <- rep(0, 2)
      prodAll <- matrix(0, 2, 2)
      dataSim <- as.list(data, all.names = TRUE)
      seed.in <- if (is.null(seed)) 
        sample.int(2^30, 1)
      else seed
      set.seed(seed.in)
      if (trace) 
        pb <- txtProgressBar(style = 3)
      for (i in 1:R) {
        if ((i%%(R/10) == 0) & trace) 
          setTxtProgressBar(pb, i/R)
        dataSim <- datagen(theta.hat, data = data)
        l1 <- floglik(theta.hat, dataSim)
        l0 <- floglik(theta.til, dataSim)
        score.hat <- if (is.null(fscore)) 
          pracma::grad(f = floglik, x0 = theta.hat, data = dataSim)
        else fscore(theta.hat, dataSim)
        obj.score.new.hat = score.hat/dpsi.dtheta
        uhh <- c(obj.score.new.hat, (l1 - l0))
        meanAll <- meanAll + uhh/R
        prodAll <- prodAll + tcrossprod(uhh)/R
      }
      if (trace) 
        close(pb)
      covAll <- prodAll * R/(R - 1) - tcrossprod(meanAll) * 
        R/(R - 1)
      i.hat <- covAll[1, 1]
      q <- covAll[2, 1]
      u <- q * sqrt(j.hat.new)/i.hat
      INF <- (1/r) * log(u/r)
      rs <- r + INF
      out <- list(r = r, INF = INF, rs = rs, theta.hat = theta.hat, 
                  se.theta.hat = se.theta.hat, info.hat = j.hat, 
                  psi.hat = psi.hat, se.psi.hat = se.psi.hat, theta.hyp = theta.til, 
                  psi.hyp = psival, seed = seed.in)
    }
  }
  out$psidesc <- psidesc
  out$R <- R
  if ((!ronly) & (abs(r) < 0.1)) {
    cat("Value under testing close to the MLE - there might be a singularity in r*\n")
    warning("Value under testing close to the MLE - there might be a singularity in r*\n")
  }
  return(structure(out, class = "rstar"))
}

wrap_my_rstar <- function(slik, h0, R=1000L, which="safe") {
  
  floglik <- function(theta, 
                      data # here data is the result of datagen()
  ) { summLik(object=data$slik, data=data$newobs, parm=theta, which=which) }
  
  datagen <- function(theta, data) {
    data$newobs <- simulate(object = data$slik, nsim=1L,theta = theta)
    data
  }
  
  
  fpsi <- function(theta) theta[names(h0)]
  rstar_resu <- my_rstar(data=list(slik=slik,  ## cf also diy2inf::my_rstar() for equivalent results
                                               newobs=t(get_from(slik,"proj_data"))),
                                     thetainit = slik$MSL$MSLE,
                                     floglik=floglik,
                                     fpsi=fpsi,
                                     # constr.opt = "alabama",
                                     datagen=datagen,
                                     psival=h0,
                                     R=R)
  chi2_LR <- rstar_resu$rs^2
  basic_LR <- rstar_resu$r^2
  
  if (FALSE) { # unidirectional tail  by Lugannani-Rice-like expression 
    # but this fails miserably, so nothing included in return value. 
    r <- rstar_resu$r
    Lapl_tail <- 1-pnorm(r)+(r/rstar_resu$utot-1) *pnorm(r)/r
  }
  
  
  data.frame(chi2_LR=chi2_LR, df=1,  p_value=1-pchisq(chi2_LR, df=1), # Lapl_tail=2*Lapl_tail,
             # keep info about basic LRT, which should match the "basicLRT" computation:
             basic=data.frame(chi2_LR=basic_LR, df=1,  p_value=1-pchisq(basic_LR, df=1)))
}



# floglik <- function(theta, 
#                     data # here data is the result of datagen()
# ) { 
#   summLik(object=data$slik, data=data$newobs, parm=theta)
# }
# 
# fpsi <- function(theta) theta["log1p.t23."]
# 
# datagen <- function(theta, data) {
#   statdens_h0 <- Infusion:::.conditional_Rmixmod(data$slik$jointdens, given=theta, expansion=1) # stat dens|ML parameter estimates
#   data$newobs <- Infusion:::.simulate.MixmodResults(statdens_h0, nsim=1L, size=1L, drop=TRUE,
#                                                     norm_or_t=Infusion:::.wrap_rmvnorm) # directly in projected space
#   data
# }
# 
