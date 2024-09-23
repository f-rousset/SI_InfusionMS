wrap_rstar <- function(slik, h0, R=1000L, which="safe") {
  
  floglik <- function(theta, 
                      data # here data is the result of datagen()
  ) { summLik(object=data$slik, data=data$newobs, parm=theta, which=which) }
  
  datagen <- function(theta, data) {
    data$newobs <- simulate(object = data$slik, nsim=1L,theta = theta)
    data
  }
  
  
  fpsi <- function(theta) theta[names(h0)]
  rstar_resu <- likelihoodAsy::rstar(data=list(slik=slik,  ## cf also diy2inf::my_rstar() for equivalent results
                                       newobs=t(get_from(slik,"proj_data"))),
                             thetainit = slik$MSL$MSLE,
                             floglik=floglik,
                             fpsi=fpsi,
                             constr.opt = "alabama",
                             datagen=datagen,
                             psival=h0,
                             R=R)
  chi2_LR <- rstar_resu$rs^2
  basic_LR <- rstar_resu$r^2
  data.frame(chi2_LR=chi2_LR, df=1,  p_value=1-pchisq(chi2_LR, df=1),
             # keep info about basic LRT, which should match the "basicLRT" computation:
             basic=data.frame(chi2_LR=basic_LR, df=1,  p_value=1-pchisq(basic_LR, df=1)))
}

get_performance <- function(ecdfs_df, # ! this is NOT do.call(rbind, ecdfs): see ecdfs_5K <- ...
                            target, plot.=TRUE, mfrow=c(2L,1L),
                            relative=TRUE) {
  if (is.matrix(target)) {
    parm <- colnames(target)
  } else parm <- names(target)
  if ("MSL" %in% colnames(ecdfs_df)) {
    MLSEv <- do.call(rbind, do.call(rbind, ecdfs_df[,"MSL"])[,"MSLE"])
    mslev <- MLSEv[,parm]
    # hessians <- do.call(rbind,ecdfs_df[,"MSL"])[,"hessian"]
    parNames <- colnames(MLSEv) 
    parmpos <- which(parNames==parm)
    # sing <- sapply(hessians, function(v) diag(solve(v))[parmpos]<0)
    loc_df <- data.frame(MSLE=mslev) 
  } else loc_df <- as.data.frame(matrix(nrow=nrow(ecdfs_df),ncol=0L))
  if ("CIs" %in% colnames(ecdfs_df)) {
    CIs_df <-  do.call(rbind, do.call(rbind, ecdfs_df[,"CIs"])[,parm])
    if ("0.95" %in% colnames(CIs_df)) # weibull...
      CIs_df <- do.call(rbind, CIs_df[,"0.95"])
    mslCI <- do.call(rbind,CIs_df[,"interval"])
    coverage <- c(mslCI=sum((is.na(mslCI[,1]) | mslCI[,1] < target) & (is.na(mslCI[,2]) | mslCI[,2] > target))/nrow(mslCI))
    if ("bootCI" %in% colnames(CIs_df)) {
      # ad hoc handling of case where boot CIs were computed for some replicates, but not all:
      # bootCIs <- mslCI[,"bootCI"]
      # whichboot <- unlist(lapply(bootCIs, inherits, what="bootci"))
      #percent <- do.call(rbind, do.call(rbind, bootCIs[whichboot])[,"percent"])[,4:5]
      percent <- do.call(rbind, do.call(rbind, CIs_df[,"bootCI"])[,"percent"])[,4:5]
      coverage <- c(coverage,
                    perc=sum(percent[,1] < target & percent[,2] > target)/nrow(mslCI))
    }
  } else coverage <- c()
  if ("abcrf" %in% colnames(ecdfs_df)) {
    abcrf_resu <- do.call(rbind,do.call(rbind, ecdfs_df[,"abcrf"])[,parm])
    loc_df$postEv <- unlist(abcrf_resu[,"expectation"])
    loc_df$postMed <- unlist(abcrf_resu[,"med"])
    postCI <- do.call(rbind, abcrf_resu[,"quantiles"])
    coverage <- c(coverage,
                  postCI=sum(postCI[,1] < target & postCI[,2] > target)/nrow(postCI))
    if ("HPDint_d" %in% colnames(abcrf_resu)) {
      HPD_int <- do.call(rbind, abcrf_resu[,"HPDint_d"])
      coverage <- c(coverage,
                    HPDint_d=sum(HPD_int[,1] < target & HPD_int[,2] > target)/nrow(HPD_int))
    }
    if ("HPDint_l" %in% colnames(abcrf_resu)) {
      HPD_int <- do.call(rbind, abcrf_resu[,"HPDint_l"])
      coverage <- c(coverage,
                    HPDint_l=sum(HPD_int[,1] < target & HPD_int[,2] > target)/nrow(HPD_int))
    }
  }
  if ("SLRTs" %in% colnames(ecdfs_df)) {
    LRTs <- do.call(rbind,do.call(rbind,ecdfs_df[,"SLRTs"])[,parm])
    basicLRTs <- do.call(rbind, LRTs[,"basicLRT"])
    coverage <- c(coverage, basicP=1-sum(basicLRTs$p_value <=0.05)/nrow(basicLRTs)) 
    if ("rawBootLRT" %in% colnames(LRTs)) {
      rawBootLRT <- do.call(rbind, LRTs[,"rawBootLRT"])
      coverage <- c(coverage, rawbootP=1-sum(rawBootLRT$p_value <=0.05)/nrow(rawBootLRT)) 
    }
    if ("BartBootLRT" %in% colnames(LRTs)) {
      BartBootLRT <- do.call(rbind, LRTs[,"BartBootLRT"])
      coverage <- c(coverage, bartlettP=1-sum(BartBootLRT$p_value <=0.05)/nrow(BartBootLRT)) 
    }
    if ( length(pstarname <- intersect(c("safeBartBootLRT", "pstarBoot"), # I very transiently used safeBartBootLRT to store pstarBoot
                                       colnames(LRTs)))) {
      pstarBoot <- do.call(rbind, LRTs[,pstarname])
      coverage <- c(coverage, pstarBoot=1-sum(pstarBoot$p_value <=0.05)/nrow(pstarBoot)) 
    }
    
    rstar_name <- intersect(c("LRT_rstar","p_rstar"), colnames(LRTs))
    if (length(rstar_name)) {
      LRT_rstar <- do.call(rbind, LRTs[,rstar_name])
      LRT_rstar <- na.omit(LRT_rstar)
      if (length(naa <- attr(LRT_rstar,"na.action"))) {
        message(paste0("NA/NaN in ",rstar_name," for ",parm," for replicate(s) ",paste(naa, collapse=","),"."))
      }
      coverage <- c(coverage, rstarP=1-sum(LRT_rstar$p_value <=0.05)/nrow(LRT_rstar)) 
    }
  }
  if (plot. && "MSL" %in% colnames(ecdfs_df) ) {
    
    oldpar <- par(mfrow=mfrow)
    xylim <- range(c(loc_df,(LOWER[parm]+UPPER[parm])/2,target))
    if (ncol(loc_df)==1L) { # loc_df has one or two column dependening on abcrf presence
      p <- ggplot(loc_df, aes(x="",y=MSLE)) + 
        geom_violin() + 
        ylim(LOWER[parm],UPPER[parm]) +
        geom_point(x = 1, y=target)
      print(p)
    } else {
      if (ncol(loc_df)==2L) {
        plot(loc_df, main=parm, xlim=xylim, ylim=xylim); abline(0,1) 
      } else {
        plot(loc_df[,c("MSLE","postMed")], main=parm, xlim=xylim, ylim=xylim); abline(0,1) # loc_df has one or two column dependening on abcrf presence
        points(loc_df[,c("MSLE","postEv")], col="grey70")
      }
      points(target, target, col="red", pch=19)
      abline((LOWER[parm]+UPPER[parm])/2, 0, col="blue") 
    }
    
    inplot <- "basicLRT"
    colors. <- c("basicLRT"="black", "rawBootLRT"="red", "BartBootLRT"="blue",
                 # "safeBartBootLRT"="red",
                 "pstarBoot"="cyan", "LRT_rstar"="green")
    pchs <- list("basicLRT"=20, "rawBootLRT"=".", "BartBootLRT"=".",
                 #"safeBartBootLRT"=20,
              "pstarBoot"=".", "LRT_rstar"=".")
    plot(ecdf(basicLRTs$p_value), xlim=c(0,1), ylim=c(0,1), main=paste("CDF of p-values for",parm),
         pch=pchs[["basicLRT"]]); abline(0,1) 
    for (st in c("LRT_rstar", "pstarBoot","BartBootLRT", "rawBootLRT" #, "safeBartBootLRT" 
                 )) {
      if (exists(st, inherits = FALSE)) {
        plot(ecdf(get(st)$p_value), add=TRUE, col=colors.[[st]], pch=pchs[[st]])
        inplot <- c(inplot,st)
      }
    }
    inplot <- intersect(names(colors.), inplot) # reordering
    if (length(inplot)>1L) {
      legend("bottomright", 
             c("basicLRT"="basic", "rawBootLRT"="raw bt.", "BartBootLRT"="Bartlett bt.",
               # "safeBartBootLRT"="safe Bart.",
               "pstarBoot"="p* bt.", "LRT_rstar"="r*")[inplot], 
             col=colors.[inplot], pch=20)
    }
    par(oldpar)
  }
  err <- (loc_df-target)
  if (relative) err <- err/(UPPER[parm]-LOWER[parm])
  if (length(loc_df$"postEv")) {
    priorMean <- (UPPER[parm]+LOWER[parm])/2
    sdpostEv <- (mean(loc_df$"postEv"^2) - mean(loc_df$"postEv")^2)/(length(loc_df$"postEv")-1L)
    var_ratio <- var(loc_df$"postEv")/var(loc_df$"MSLE")
    
    info <- c(corr=cor(loc_df$"MSLE",loc_df$"postEv"),
              var_ratio=var_ratio,
              pos2Prior=((mean(loc_df$"postEv")-priorMean)/
                             (target-priorMean+c(-1,1)*2*sdpostEv))[1:2]) # with fix for target-priorMean=0
  } else info <- NULL
  
  list(bias=colMeans(err), RMSE=sqrt(colMeans(err^2)), cvrg=coverage, info=info)
}

ad_hoc_fn <- function(ecdfs_df, # a data frame with a complex structure
                      DGparlist, parNames=NULL, 
                      par.grid=NULL, 
                      ...) {
  if (is.null(par.grid)) {
    if (is.null(parNames)) parNames <- names(DGparlist)
    resu <- vector("list", length(parNames))
    names(resu) <- parNames
    for (parm in parNames) {
      resu[[parm]] <- 
        unlist(get_performance(ecdfs_df, 
                             target=unlist(DGparlist[parm]), ...))
    }
  } else {
    if (is.null(parNames)) parNames <- colnames(par.grid)
    resu <- vector("list", length(parNames))
    names(resu) <- parNames
    for (parm in parNames) {
      resu[[parm]] <- 
        unlist(get_performance(ecdfs_df, 
                               target=as.matrix(par.grid[,parm,drop=FALSE]), ...))
    }
  }
  resu
}

diagnose <- function(ecdfs_df, parm) {
  mslev <- do.call(rbind, do.call(rbind, ecdfs_df[,"MSL"])[,"MSLE"])[,parm]
  LRTs <- (do.call(rbind,do.call(rbind,do.call(rbind,ecdfs_df[,"SLRTs"])[,parm])))
  data.frame(MSLE=mslev, chi2_LR=LRTs[,"chi2_LR"])
}

summary_inference_K <- function(slik_siz=upsliks[[length(upsliks)]], upsliks, # either upsliks or slik_siz non-NULL
                                S_obs=get_from(slik_siz,"raw_data"),  
                                nsim, nb_cores_LRboot=1L, ncores_abcrf=mydetectCores(), abcrf_proj_K=NULL, 
                                h0s, reftable_abcrf, ret_HPDint, 
                                verbose=FALSE, 
                                rstar_R=1000L,
                                #profile_first=FALSE, # ideally not needed too
                                ii,
                                ...) {
  ### hack to handle interrupted jobs, for out-of memory errors in abcrf run for example
  tmpfile <- paste0("interrupted_summ_",ii,".rda")
  chk <- suppressWarnings(try(load(tmpfile), silent=TRUE))
  if (inherits(chk,"try-error")) {
    resu <- NULL
  } else print(paste0(tmpfile," loaded.")) # we have a partial resu list
  ### 
  # if (profile_first) {
  #   if (verbose) cat(crayon::yellow("\nComputing profiles for more secure maximization... "))
  #   plot1Dprof(slik_siz, do_plot=interactive()) 
  # } 
  # For old objects:
  if (is.null(resu)) {
    if (is.null(slik_siz$CIobject)) slik_siz$CIobject <- list2env(list(CIs=NULL),
                                                                  parent=emptyenv())
    if (verbose) cat(crayon::yellow("\nComputing SLik CIs... "))
    CIs <- allCIs(slik_siz, verbose=verbose, nsim=nsim, nb_cores=nb_cores_LRboot)$CIs
    resu <- list(MSL=as.list(slik_siz$MSL), CIs=CIs, S_obs=S_obs)
    save(resu, file=tmpfile)
  }
  if ( ! inherits(slik_siz$jointdens,"MAF")) resu$nbCluster <- slik_siz$jointdens@nbCluster
  #the boostrap is the slow step.
  
  if ( ! is.null(h0s) && is.null(resu$SLRTs)) {
    if (verbose) cat(crayon::yellow("\nComputing SLRTs... "))
    if (nsim>1L) cat(crayon::yellow(paste("(!) with bootstrap,",nsim,"simulations:")))
    SLRTs <- vector("list", length(h0s))
    names(SLRTs) <- names(h0s)
    latest_SLRTs <- SLRTs # never reset to a list of NULLs => stores the latest 'best' results &
                          # is returned in case of manual interrupt:
    tryCatch(
      while (is.null(SLRTs[[1]])) {
        for (it in seq_along(SLRTs)) {
          if (is.null(SLRTs[[it]])) {
            cat(" ",crayon::yellow(names(h0s)[it]))
            h0 <- h0s[[it]]
            h0names <- names(h0)
            prevmaxlogL <- slik_siz$MSL$maxlogL
            slrt <- SLRT(slik_siz, h0=h0, nsim=nsim, 
                         nb_cores=nb_cores_LRboot, #passed through the ... all the way down to .ecdf_2lr() 
                         ...)
            is_h0_at_bound <-  any(c(h0-slik_siz$lower[h0names],slik_siz$upper[h0names]-h0)<1e-14)
            ## If h0 is at bound (or beyond), bobyqa will be used with a local init pushed inside
            ## and the new logL may be below that of init (notably in the 'beyond bounds' case)
            ## So that an infinite loop would result. 
            ## In general we should avoid setting h0 at bound in the performance checks, but this has
            ## occurred by oversight, causing an infinite loop.
            ## Hence the check here, and an added warning in SLRT().
            ## But...
            
            if ( ! is_h0_at_bound) {
              while (attr(slrt,"MSL_updated")) {
                newmaxlogL <- slik_siz$MSL$maxlogL
                if (verbose) print(paste0("new maxlogL: ", newmaxlogL))
                if (newmaxlogL <= prevmaxlogL) {
                  # Apparently I still have cases where the profile found a new max,
                  # so I re-optimize but the 'new new max' is not above the old max.
                  # The solution is at a boundary and bobyqa is bobyqa...
                  # This can make a long loop (synthetic/ evaluations, #72)
                  break()
                } else prevmaxlogL <- newmaxlogL
                SLRTs[1:length(h0s)] <- list(NULL) # all: consider the case where 1st h0 induces updating AFTER 2nd one did so.
                slrt <- SLRT(slik_siz, h0=h0s[[it]], nsim=nsim, nb_cores=nb_cores_LRboot, ...)
              }
            }
            SLRTs[[it]] <- slrt
            latest_SLRTs[[it]] <- slrt
          }
          if (is.null(SLRTs[[1]])) break
        }
      },
      interrupt = function(e){ 
        message("Manual interrupt of loop!")
        SLRTs <<- latest_SLRTs
      }
    )
    if (requireNamespace("likelihoodAsy") && rstar_R>0L) for (it in seq_along(SLRTs)) {
      SLRTs[[it]]$LRT_rstar <- wrap_my_rstar(slik=slik_siz, h0=h0s[[it]]) 
    }
    resu$SLRTs <- SLRTs
    resu$MSLpost=as.list(slik_siz$MSL) # suitable for diagnosing SLRTs
    save(resu, file=tmpfile)
  }
  if (is.null(resu$abcrf) &&
      ! (is.null(abcrf_proj_K) || is.null(reftable_abcrf))
    ) {
    if (verbose) cat(crayon::yellow("\nPerforming abcrf analyses... "))
    resu$abcrf <- abcrf_replicate(Sobs=S_obs, # fn def'd in abcrf_addl_defs.R
                                  abcrf_proj=abcrf_proj_K,
                                  reftable_abcrf=reftable_abcrf,
                                  ret_HPDint = ret_HPDint,
                                  ncores_abcrf=ncores_abcrf)
  } else if (verbose) {
    cat(crayon::yellow(
      "\nNot performing abcrf analyses because ",    
      paste(c(" ! is.null(resu$abcrf)",
              "is.null(abcrf_proj_K)", 
              "is.null(reftable_abcrf)")[c(! is.null(resu$abcrf),
                                           is.null(abcrf_proj_K),
                                           is.null(reftable_abcrf))], collapse=", ")
    ))
  }
  unlink(tmpfile)
  resu
}


compare_intervals <- function(ecdfs_df, parm) {
  mslCI <- do.call(rbind, do.call(rbind, do.call(rbind, ecdfs_df[,"CIs"])[,parm])[,"interval"])
  resu <- as.data.frame(mslCI)
  if ("abcrf" %in% colnames(ecdfs_df)) {
    abcrf_resu <- do.call(rbind,do.call(rbind, ecdfs_df[,"abcrf"])[,parm])
    if ("HPDint_d" %in% colnames(abcrf_resu)) {
      HPDint_d <- do.call(rbind, abcrf_resu[,"HPDint_d"])
      resu <- cbind(resu, HPDint_d=as.data.frame(HPDint_d))
    }
    if ("HPDint_l" %in% colnames(abcrf_resu)) {
      HPDint_l <- do.call(rbind, abcrf_resu[,"HPDint_l"])
      resu <- cbind(resu, HPDint_l=as.data.frame(HPDint_l))
    }
  }
  plot(resu[,1],ylim=range(resu), main=paste0("Bounds of intervals for '",parm,"'"))
  points(resu[,2])
  whch <- c(TRUE, "HPDint_d" %in% colnames(abcrf_resu), "HPDint_l" %in% colnames(abcrf_resu))
  if (whch[2]) {
    points(HPDint_d[,1], col="red")
    points(HPDint_d[,2], col="red")
  }  
  if (whch[3]) {
    points(HPDint_l[,1], col="blue")
    points(HPDint_l[,2], col="blue")
  }  
  resu
  legend("right", c("summ.lik", "HPD default", "HPD locfit")[whch], col=c("black","red","blue")[whch], pch=1)
}


