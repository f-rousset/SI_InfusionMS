
calc_sliks <- function(ii, S_obs,reftable_sizes, control.Simulate=NULL,
                       LOWER,UPPER,
                       upsliks=NULL, # possible incomplete one
                       save_upsliks=TRUE,
                       verbose=interactive(),
                       cluster_args=list(project=list(num.threads=mydetectCores()-1L)), # -> refine -> ranger parallelisation
                       methodArgs=list(),
                       constr_crits=NULL,
                       scaDGparlist,
                       init_indiv_reft_size=min(reftable_sizes[1]/5L,1000L), 
                       update_projectors=TRUE,
                       reduce_projectors=TRUE,
                       using="Rmixmod",
                       workflow_design=get_workflow_design(length(LOWER)),
                       Simulate="DIYABC_reftable_wrapper",
                       ...) {
  on_genotoul <- (user <- Sys.info()["user"]=="frousset") # control of gcinfo()
  if (is.null(names(reftable_sizes))) stop("'reftable_sizes' must be a *named* vector.")
  if (is.logical(verbose)) verbose <- list(most=verbose, rparam=FALSE)
  blocksizes <- diff(c(0L,reftable_sizes))
  subblock_nbr <- workflow_design$subblock_nbr
  subblock_size <- workflow_design$subblock_size
  # the total number of sub-blocks over all blocks is O(npar) AND blocksize is a multiple of subblocksize 
  
  if ( ! length(upsliks)) {
    if (is.null(control.Simulate)) stop("'control.Simulate' must be provided")
    parNames <- names(LOWER)
    parsp <- init_reftable(lower=LOWER, upper=UPPER, nUnique = init_indiv_reft_size,
                           constr_crits = constr_crits)
    # control.Simulate$fix_parsTable and $fix_canon.grid are called in DIYABC_reftable_wrapper()

    # To check the ranges:
    # apply(apply(parsp,1L, canonize),1,range)

    reftable_raw <- add_reftable(    Simulate=Simulate,
                                     #
                                     par.grid=parsp, # canonized internally ( not externally as we want the reftable in composite space)
                                     #
                                     control.Simulate=control.Simulate
    )
    dprojectors <- vector("list", length(parNames))
    names(dprojectors) <- paste0("p",parNames)
    for (st in parNames) {
      dprojectors[[paste0("p",st)]] <- project(st,stats=names(S_obs), data=reftable_raw,verbose=interactive(),
                                               methodArgs=methodArgs, keep_data=FALSE)
    }
    
    dprojSimuls <- project(reftable_raw,projectors=dprojectors,verbose=FALSE,is_trainset = TRUE, use_oob = TRUE)
    dprojSobs <- project(S_obs,projectors=dprojectors, is_trainset = FALSE, ext_projdata=reftable_raw)
    
    slik_siz <- infer_SLik_joint(dprojSimuls, stat.obs=dprojSobs, verbose=verbose,
                                 constr_crits=constr_crits, using=using)
    slik_siz$projdata <- reftable_raw
    
    slik_siz <- MSL(slik_siz, CIs=FALSE, eval_RMSEs = FALSE)
    
    if (FALSE) { 
      plot1Dprof(slik_siz,
                 decorations=function(par) {
                   # yrange ref working with Infusion version > 2.1.10:
                   if ( ! is.null(scaDGparlist)) { # not for 'actual_obs'
                     yrange <- get("yrange", envir=parent.frame()) # parent frame is the plot1Dprof() evaluation frame.
                     points(y=yrange[1],x=scaDGparlist[[par]],pch=20,cex=2,col="red")
                     points(y=yrange[2],x=slik_siz$MSL$MSLE[par],pch=20,cex=2,col="blue")
                   }
                 })
    }
    ntot_this <- min(2L*subblock_size, 
                     reftable_sizes[1]-nrow(slik_siz$logLs)) # where 2L*subblock_size is >= reftable_sizes[1]
    if (on_genotoul) options(error=quote({dump.frames(to.file=TRUE, include.GlobalEnv = TRUE)}))
    if (ntot_this>0L) { # First refine
      cluster_args_ini <- cluster_args
      cluster_args_ini$ranger <- NULL # no ranger parallelisation for first refines
      slik_siz <-  refine(slik_siz, # maxit=maxit,
                          workflow_design=workflow_design,
                          ntot=ntot_this,
                          CIs=FALSE, eval_RMSEs=FALSE, # update_projectors=TRUE,
                          map.asp=1, # tempo fix spaMM <= 4.4.0
                          cluster_args=cluster_args_ini, 
                          verbose=verbose,
                          precision=0) #
    }
    # next, we finish the first block after the first refine:
    ntot_this <- reftable_sizes[1]-nrow(slik_siz$logLs) #Using actual reftable size rather than previous target size 'subblock_size'
    if (ntot_this>0L) {
      maxit <- max(1L, as.integer(ntot_this/subblock_size))
      slik_siz <-  refine(slik_siz, maxit=maxit, #loop until reftable_sizes[1] is reached.
                          workflow_design=workflow_design,
                          ntot=ntot_this, 
                          map.asp=1, # tempo fix spaMM <= 4.4.0
                          CIs=FALSE, eval_RMSEs=FALSE, # update_projectors=TRUE,  
                          cluster_args=cluster_args_ini, # no ranger parallelisation for first refines 
                          verbose=verbose,
                          precision=0) #
    }
    
    upsliks <- vector("list", length(reftable_sizes))
    names(upsliks) <- names(reftable_sizes)
    upsliks[[1]] <- slik_siz
    newsizes <- names(reftable_sizes)[-1]

  } else {
    # refine pre-existing upsliks
    slik_siz <- upsliks[[length(upsliks)]]
    if (is.null(control.Simulate)) {
      control.Simulate <- attr(slik_siz$logLs,"control.Simulate")
    } # non-NULL value allows cross-platform workflows.
    newsizes <- setdiff(names(reftable_sizes), names(upsliks))
    supp <- vector("list", length(newsizes))
    names(supp) <- newsizes
    upsliks <- c(upsliks, supp) # add empty elements to be filled by next loop
  }
  
  
  
  it <- 0L
  for (newsize in newsizes) { # loop to get fits beyond the reftable_sizes[1] one
    if (on_genotoul) {
      gcinfo(TRUE)
      gc()
      gcinfo(FALSE)
    }
    it <- it+1L
    maxit <- as.integer(ceiling(blocksizes[newsize]/subblock_size))
    
    # if (reftable_sizes[newsize]>60000L) {
    #   message("**update_projectors_it set to FALSE because of likely memory limitation**\n")
    #   update_projectors_it <- FALSE 
    # } else if (length(update_projectors)==1L && is.logical(update_projectors)) {
    #   if (update_projectors) {
    #     if (reftable_sizes[newsize]>20000L) {
    #       # Convert TRUE into numeric value for size of reftable at end of refine => induces 'final' projection:
    #       update_projectors_it <- blocksizes[newsize]
    #       # ...BUT not for too large ones bc I limit the memory:
    #     } else update_projectors_it <- TRUE
    #   }
    # } else update_projectors_it <- update_projectors
    
    # update_projectors_it <- TRUE   
    
    locverbose <- verbose
    locverbose$proj <- locverbose$movie <- FALSE
    message(paste0("refine for new size ", newsize,":\n"))
    slik_siz <- refine(slik_siz, maxit=maxit, ntot=blocksizes[newsize], 
                       workflow_design=workflow_design,
                       control.Simulate=control.Simulate,
                       verbose=locverbose, 
                       map.asp=1, # tempo fix spaMM <= 4.4.0
                       CIs=FALSE, eval_RMSEs=FALSE, # update_projectors=update_projectors_it,  
                       cluster_args=cluster_args, precision=0) #
    if (newsize %in% names(upsliks)) upsliks[[newsize]] <- slik_siz
  }
  # This reduces the projectors in the *same* environment found in *all* elements of upsliks:
  if (reduce_projectors) deforest_projectors(slik_siz) # otherwise saving may be *slow* and *bulky* 
  if (save_upsliks) save(upsliks, file=paste0("saved_fits/upsliks_",ii,".rda"))
  upsliks
}

# I need a perdata_perf_workflow to get CIs etc from the SLiks

get_upsliks <- function(ii, reftable_sizes, more_refines, workflow_design,
                        ...) { # ... passing S_obs,
                                # 'LOWER', 'UPPER', 'control.Simulate',
                                # 'constr_crits', 'constr_tuning', 'verbose', 'workflow_design'
  if (exists("upsliks",envir = .GlobalEnv)) {
    cat(crayon::yellow(paste0(" 'upsliks' already exists in .GlobalEnv")))
  } else {
    if (on_genotoul <- (user=="frousset")) {
      (chk <- safe_load(paste0("./../saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
    } else (chk <- safe_load(paste0("./saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
    
    if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
      cat(crayon::yellow(paste0("generating 'upsliks'...")))
      upsliks <- calc_sliks(ii, reftable_sizes=reftable_sizes, 
                            workflow_design=workflow_design, ...) 
    } else cat(crayon::yellow(paste0(" 'upsliks' loaded from disk")))
  }
  # an upsliks is now available
  
  if (more_refines) {
    cat(crayon::yellow(paste0(" ...augmenting 'upsliks' for 'more_refines'= ",
                              more_refines,"...\n")))
    more_reftable_sizes <- tail(reftable_sizes,1)+ seq_len(more_refines)*workflow_design$refine_blocksize
    names(more_reftable_sizes) <- paste0(more_reftable_sizes/1000L,"K")
    reftable_sizes <- c(reftable_sizes, more_reftable_sizes)
  }
  newsizes <- setdiff(names(reftable_sizes), names(upsliks))
  if (length(newsizes)) {
    if (length(newsizes)) cat(crayon::yellow(paste0("...augmenting 'upsliks' for sizes ",
                                                    paste(newsizes, collapse=" "),"...\n")))
    upsliks <- calc_sliks(ii, reftable_sizes=reftable_sizes, upsliks=upsliks, 
                          workflow_design=workflow_design, ...) 
  }

  upsliks
}
