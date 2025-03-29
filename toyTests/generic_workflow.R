# NOT the same generic_workflow as for diy2inf, but same logic

{  # Run this in all cases
  
  ## Keep it handy; but don't run the risk of saturating disks...
  # if ( ! interactive()) options(error=quote({dump.frames(to.file = TRUE, include.GlobalEnv=TRUE)}))
  {
    library("Infusion")
    library("mvtnorm")
    
    user <- Sys.info()["user"]
    
    if (on_genotoul <- (user=="frousset")) {
      on_old_genotoul <- length(grep("genologin|node",Sys.info()["nodename"]))
      if (on_old_genotoul) {
        dir_generic <- "/work/frousset/Infusionplus/diy2inf/generic/"
      } else dir_generic <- "/work/user/frousset/Infusionplus/diy2inf/package/R/"
    } else if (.Platform$OS.type=="windows") {
      # user=="francois.rousset"
      dir_generic <- "D:/home/francois/travail/stats/Infusionplus/diy2inf/package/R/"
    } else if (user=="francois") {
      dir_generic <- "/mnt/d/home/francois/travail/stats/Infusionplus/diy2inf/package/R/"
      Infusion.options(is_devel_session=TRUE)
    } else {
      print(Sys.info())
      stop("You need to adapt this generic workflow to your machine.")
    }
    
    {
      # library(diy2inf)
      source(paste0(dir_generic,"safe_simul_utils.R")) 
      source(paste0(dir_generic,"Infusion_addl_defs.R"))
      source(paste0(dir_generic,"performance_fns.R")) # these one at least are useful for ricker
      # source(paste0(dir_generic,"my_rstar.R")) #
      source(paste0(dir_generic,"output_files_manips.R")) # df_extract...
    }
    
    if (on_genotoul) {
      (thisfilepath <- paste0(getwd(),"/"))
      cmdl_args <- commandArgs(trailingOnly = TRUE)
      # bash script has '--args jobNbr=$job threadNbr=$threadNbr 
      #      and perhaps refine=TRUE or =2...
      cmdl_args <- eval(parse(text=paste0("list(", paste(cmdl_args,collapse=", "),")")))
      if (is.null(cmdl_refine <- cmdl_args$refine)) cmdl_refine <- FALSE
      if (is.null(cmdl_resummarize <- cmdl_args$resummarize)) cmdl_resummarize <- FALSE
      if (is.null(cmdl_fast_reClu <- cmdl_args$fast)) cmdl_fast_reClu <- FALSE 
      if (is.null(cmdl_reproject <- cmdl_args$reproject)) cmdl_reproject <- FALSE
      if (is.null(cmdl_boot_nsim <- cmdl_args$boot_nsim)) cmdl_boot_nsim <- 199L
      if (cmdl_reproject || cmdl_resummarize) {
        if ( ! is.null(procs_per_job <- cmdl_args$threadNbr) &&
             procs_per_job>1L ) warning("cmdl_args$threadNbr>1L with cmdl_reproject || cmdl_resummarize")
        # => in principle I prevent this at the level of the sbatch script.
        # procs_per_job <- 1L
      } else if (is.null(procs_per_job <- cmdl_args$threadNbr)) procs_per_job <- 1L # which is also mydetectCores()'s ad-hoc return value for genotoul (a safety fix)  
      if (is.null(cmdl_jobNbr <- cmdl_args$jobNbr)) cat(crayon::underline("'jobNbr' missing from commandArgs on genotoul."))
      if ( is.null(cmdl_jobNbr)) {
        jobsubdir <- "./"
      } else jobsubdir <- paste0("job_",cmdl_jobNbr,"/")
    } else {
      if (interactive()) {
        if ( ! exists("cmdl_refine")) cmdl_refine <- FALSE
        if ( ! exists("cmdl_resummarize")) cmdl_resummarize <- FALSE
        if ( ! exists("cmdl_reproject")) cmdl_reproject <- FALSE
        if ( ! exists("cmdl_jobNbr")) cmdl_jobNbr <- NULL
        if ( ! exists("cmdl_fast_reClu")) cmdl_fast_reClu <- FALSE
      } else {
        cmdl_refine <- cmdl_resummarize <- cmdl_fast_reClu <- FALSE
        cmdl_jobNbr <- NULL
      }
      if (Infusion:::.inRstudio()) {
        (thisfilepath <-  paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/") )
        setwd(thisfilepath)
      } else {
        "make sure the work dir to be the source directory. Then:"
        (thisfilepath <- paste0(getwd(),"/"))
      }
      procs_per_job <- mydetectCores() # max(1L,mydetectCores()-1L)
    }
    
    modelpath_nickname <- tail(strsplit(thisfilepath,"/")[[1]],1) 
    final_MAF <- length(grep("_final_MAF", modelpath_nickname)) > 0L
    
    (savedir <- thisfilepath) # this makes it the working directory...     with final '/'
    
    cat(crayon::yellow("Working directory is currently "),getwd(),"\n")
    
  }
  
  {
    if (length(grep("MVN6",thisfilepath))) {
      DIM <- 3L    
    } else DIM <- 5L
    npar <- DIM*(DIM+1L)/2L
    
    (statnames <-  paste0(
      "S",
      unlist(lapply((seq(DIM)), function(v) seq(v))),
      unlist(lapply(seq(DIM), function(v) rep(v,v)))
    ))
    
    if (length(grep("identif",modelpath_nickname))) {
      sample1 <- function(cholTr, verbose=FALSE) {
        Sig <- diag(DIM)
        Sig[upper.tri(Sig,TRUE)] <- cholTr
        d <- diag(Sig)   
        Sig <- diag(x=d*d)    #     so that only diagonal elements will be identifiable
        covm <- cov(mvtnorm::rmvnorm(50,mean=rep(0,DIM),sigma=Sig))
        if (verbose) {
          print(Sig)
          print(covm)
        }
        covm[upper.tri(covm,diag=TRUE)]
      }
    } else {
      sample1 <- function(cholTr, verbose=FALSE) {
        Sig <- diag(DIM)
        Sig[upper.tri(Sig,TRUE)] <- cholTr
        Sig <- crossprod(Sig)
        covm <- cov(mvtnorm::rmvnorm(50,mean=rep(0,DIM),sigma=Sig))
        if (verbose) {
          print(Sig)
          print(covm)
        }
        covm[upper.tri(covm,diag=TRUE)]
      }
    }
    
    sample_fn <- function(parsTable, verbose=FALSE) {
      resu <- apply(parsTable, 1L, sample1, verbose=verbose)
      resu <- data.frame(t(resu))
      colnames(resu) <- statnames
      resu
    }
    
    { # variables always needed
      # Data-generating parameters cholDGPv
      set.seed(123)
      DGcovmat <- rWishart(1,df=10,Sigma=diag(5L))[,,1]
      cholDGP <- chol(DGcovmat)
      cholDGP <- cholDGP[1:DIM,1:DIM]
      cholDGPv <- cholDGP[upper.tri(cholDGP,TRUE)] # 
      (parnames <-  paste0(
        "V",
        unlist(lapply((seq(DIM)), function(v) seq(v))),
        unlist(lapply(seq(DIM), function(v) rep(v,v)))
      ))
      names(cholDGPv) <- parnames
      h0s <- as.list(cholDGPv)
      for (it in seq_along(h0s)) {names(h0s[[it]]) <- names(h0s)[it]}
    }
    
    chk <- try(load("S_obs_table.rda"))
    if (inherits(chk,"try-error")) {
      set.seed(234)
      S_obs_table <-  t(replicate(1000L, {sample1(cholDGPv)}))
      colnames(S_obs_table) <- statnames
      save(S_obs_table, file="S_obs_table.rda")
    }
    
    chk <- try(load("saves4Rmd.rda"))
    if (inherits(chk,"try-error")) { # variables always needed: recompute or load
      if (FALSE) {
        # a quick look at realistic parameter ranges:
        apply(t(replicate(1000, {
          rcovmat <- rWishart(1,df=10,Sigma=diag(5L))[,,1]
          cholr <- chol(rcovmat)
          cholr[upper.tri(cholr,TRUE)] #
        })),2L,range)
      }
      # => diagonal terms in (roughly) 0.5, 6 and non-diag ones in -4, 4  
      LOWER <- UPPER <- cholDGP
      LOWER[upper.tri(LOWER,TRUE)] <- 0.5
      LOWER[upper.tri(LOWER,FALSE)] <- -4
      UPPER[upper.tri(LOWER,TRUE)] <- 6
      UPPER[upper.tri(LOWER,FALSE)] <- 4
      LOWER <- LOWER[upper.tri(LOWER,TRUE)]
      UPPER <- UPPER[upper.tri(UPPER,TRUE)]
      names(LOWER) <- names(UPPER) <- parnames
      save(LOWER,UPPER,cholDGPv, file="saves4Rmd.rda")
    }
    
    
    {
      
      workflow_design <- get_workflow_design(npar)
      (default_indiv_reft_size <- workflow_design$final_reft_size)
      (refine_blocksize <- workflow_design$refine_blocksize) # eg 7 param -> blocksize =4000
      # Save 3 or 4 sizes:
      ## these lines must match those in calc_sliks():
      ##
      (reftable_sizes <- workflow_design$reftable_sizes) # ... 4000 8000 12000...
      (init_indiv_reft_size <- workflow_design$init_reft_size) # ...4000/5 => 800 ...
      # (reftable_sizes <- reftable_sizes[reftable_sizes>=2L*init_indiv_reft_size]) # ... 4000 8000 12000... 
      (reftable_sizes <- unique(c(reftable_sizes,default_indiv_reft_size))) # ... 4000 8000 12000 14000...
      (names(reftable_sizes) <- paste0(reftable_sizes/1000L,"K")) # "4K" "8K" "12K" "14K"
      # max_table_size NOT set here bc the relevant one may depend on the refines of previous sessions,
      # not only from the above info in the current session. 
      
    }
    
  }
  
  
}

parnames <- names(LOWER)

NREPL <- nrow(S_obs_table)
if (is.null(cmdl_jobNbr)) { # 
  upto <- NREPL # min(1L, NREPL)
  ii <- prev_latest_one <- 0L # virgin directory, except if:
  # ii <- prev_latest_one <- 175L # continuing simulation: the index of the last simulated 
  cat(crayon::yellow("Starting from the ",prev_latest_one+1L,"th replicated sample.\n"))
} else {
  upto <- min(cmdl_jobNbr, NREPL)
  ii <- prev_latest_one <- cmdl_jobNbr-1L 
} 


# ! Next block contains inferences !
if (ii < NREPL) { # so that setting prev_latest_one to NREPL prevents the creation of a new summaries_list
  Infusion.options(is_devel_session= user %in% c("francois.rousset","francois")) # position handling package reinstallations
  verboses <- list(most=interactive() || on_genotoul,
                   rparam=interactive() || on_genotoul,
                   cloud_parm=NULL)
  cluster_args <- list(project=list(num.threads=procs_per_job)) # but there is no projection in this specific script
  
  if (for_timings <- FALSE) { # Beware of the else !
    Rprof(filename="Rprof.fit_real_data.out", 
          event="elapsed", # non-default under linux, necessary to count time spect in system() calls 
          interval=0.5) 
    ii <- 1L
    cat(crayon::bgYellow$black(ii)," for timings\n"); system("sleep 1")
    S_obs <- S_obs_table[ii,]
    set.seed(234L+ii)
    Infusion.options(is_devel_session= user %in% c("francois.rousset","francois")) # position handling package reinstallations
    t_begin <- Sys.time()
    ## Construct initial reference table:
    parsp_j <- init_reftable(lower=LOWER,upper=UPPER,
                             nUnique =  workflow_design$init_reft_size)
    
    set.seed(456L+ii)
    dsimuls <- add_reftable(,Simulate="sample_fn",par.grid=parsp_j,verbose=FALSE)
    # head(dsimuls)
    
    verboses <- list(most=interactive() || on_genotoul,
                     rparam=interactive() || on_genotoul)
    
    slik_j <- infer_SLik_joint(dsimuls,stat.obs=S_obs,verbose=verboses)
    slik_j <- MSL(slik_j, eval_RMSEs = FALSE, CIs=FALSE)
    for (siz in names(reftable_sizes)) { # ad hoc as sizes does not contain 1K
      slik_j <- refine(slik_j, eval_RMSEs = FALSE, CIs=FALSE,
                       verbose=verboses)
    } 
    t_end <- Sys.time()
    cat(paste0("workflow took ", signif(difftime(t_end,t_begin),4),"s\n"))
    Rprof(NULL)
    cat(crayon::yellow("browser() at end of block for timings."))
    browser()
  } else {
    TMP <- summaries_list <- replicate(
      upto-prev_latest_one, {
        ii <<- ii+1L
        cat(crayon::bgYellow$black(ii),"\n")
        S_obs <- S_obs_table[ii,]
        set.seed(234L+ii)
        Infusion.options(is_devel_session= user %in% c("francois.rousset","francois")) # position handling package reinstallations
        if (length(grep("pow29",thisfilepath))) { # e.g. "[...]/nbClu/pow29/N_7from17"
          Infusion.options(nbClu_pow_rule_fn= function(nr, ...) ceiling(nr^0.29))
          message("***using alternative pow rule nr^0.29***")
        }
        if (length(grep("pow33",thisfilepath))) { # e.g. "[...]/nbClu/pow33/N_7from17"
          Infusion.options(nbClu_pow_rule_fn= function(nr, ...) ceiling(nr^0.33))
          message("***using alternative pow rule nr^0.33***")
        }
        if (length(grep("nr8",thisfilepath))) { # e.g. "[...]/nbClu/nr8/N_7from17"
          Infusion.options(maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) { 
            (nr+8L)%/%(nc*(nc+3L)*4L+8L) 
          })
          message("***using alternative maxnbCluster ~ P/8 ***")
        }
        if (length(grep("nr3",thisfilepath))) { # e.g. "[...]/nbClu/nr3/N_7from17"
          Infusion.options(maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) { 
            (nr+3L)%/%(nc*(nc+3L)*3L/2L+3L) 
          })
          message("***using alternative maxnbCluster ~ P/3 ***")
        }
        
        if (cmdl_resummarize || final_MAF || cmdl_fast_reClu) {
          if (exists("upsliks",envir = .GlobalEnv)) {
            cat(crayon::yellow(paste0(" 'upsliks' already exists in .GlobalEnv")))
          } else {
            if (on_genotoul <- (user=="frousset")) {
              (chk <- safe_load(paste0("./saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
            } else (chk <- safe_load(paste0("./saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
            
            if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
              # there should be an error message here
            } else cat(crayon::yellow(paste0(" 'upsliks' loaded from disk.\n")))
            
          }
          max_table_size <- tail(names(upsliks),1L)
          if (cmdl_fast_reClu) { # for reclustering without reprojecting
            message("***** cmdl_fast_reClu => nbCluster set to ", cmdl_fast_reClu," ***** ")
            upsliks[[max_table_size]] <- recluster(upsliks[[max_table_size]],
                                                   nbCluster=cmdl_fast_reClu) 
          }
          if (length(grep("_MAF", modelpath_nickname))) { # not used for MVN toy tests
            config_mafR(torch_device = "cuda")
            if (final_MAF && ! inherits(upsliks[[max_table_size]]$jointdens,"MAF")) {
              upsliks[[max_table_size]] <- recluster(upsliks[[max_table_size]], using="c.mafR") 
              upsliks[[max_table_size]] <- 
                save_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", 
                          ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))  
              upsliks  <- upsliks[max_table_size]
              # save(upsliks, file=paste0("./saved_fits/upsliks_",ii,".rda"))
              # : only after deforest_projectors(...)
              message("final_MAF: MSL() called before computation of summaries.")
            } else upsliks[[max_table_size]] <- 
                load_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", 
                          ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))
            
          } 
          if (cmdl_resummarize) {
            slik_j <- upsliks[[max_table_size]]
          }
          # if (length(grep("_MAF",modelpath_nickname))) {
          #   upsliks[[max_table_size]] <- save_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))
          # }
        } else {
          # message("Projecting with num.trees=1000")
          # Infusion.options(proj_methodArgs=list(num.trees=1000L))
          ## in N_7from17 replicate 1 is 'too nice' and replicate 2 is 'difficult'
          if (length(grep("_MAF",modelpath_nickname))) {
            using <- "c.mafR"
            config_mafR(torch_device = "cuda")
            message("************** 'c.mafR' is being used. ***************** ")
            # don't forget the export in the bash script ...
          } else using <- "Rmixmod"
          
          S_obs <- S_obs_table[ii,]
          
          ################################################################
          verboses <- list(most=interactive() || on_genotoul,
                           rparam=interactive() || on_genotoul)
          
          { # code trying to automate refine jobs.
            if (exists("upsliks",envir = .GlobalEnv)) {
              cat(crayon::yellow(paste0(" 'upsliks' already exists in .GlobalEnv")))
            } else {
              ## no job subdir for the MVN toy tests:
              #if (on_genotoul <- (user=="frousset")) {
              #  (chk <- safe_load(paste0("./../saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
              #} else 
              (chk <- safe_load(paste0("./saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
              
              if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
                cat(crayon::yellow(paste0("generating 'upsliks'...")))
                # original code for the workflow
                upsliks <- vector("list", length(reftable_sizes))
                names(upsliks) <- names(reftable_sizes)
                if (DIM==5L) upsliks <- upsliks[c("28K","44K")] # ad hoc
                
                ## Construct initial reference table:
                parsp_j <- init_reftable(lower=LOWER,upper=UPPER,
                                         nUnique =  workflow_design$init_reft_size)
                
                set.seed(456L+ii)
                dsimuls <- add_reftable(,Simulate="sample_fn",par.grid=parsp_j,verbose=FALSE)
                # head(dsimuls)
                
                slik_j <- infer_SLik_joint(dsimuls,stat.obs=S_obs,verbose=verboses)
                slik_j <- MSL(slik_j, eval_RMSEs = FALSE, CIs=FALSE)
                for (siz in names(reftable_sizes)) { ## ad hoc as sizes does not contain 1K
                  ## but (later added) explicit control by ntot is more robust. 
                  slik_j <- refine(slik_j, ntot=reftable_sizes[siz]-nrow(slik_j$logLs),
                                   verbose=verboses, eval_RMSEs = FALSE, CIs=FALSE)
                  if (siz %in% names(upsliks)) upsliks[[siz]] <- slik_j
                } 
                save(upsliks,file=paste0("upsliks_",ii,".rda"))
              } else cat(crayon::yellow(paste0(" 'upsliks' loaded from disk")))
            }
            # an upsliks is now available
            more_refines <- as.integer(cmdl_refine)
            save_upsliks <- TRUE
            if (more_refines) {
              cat(crayon::yellow(paste0(" ...augmenting 'upsliks' for 'more_refines'= ",
                                        more_refines,"...\n")))
              curr_size <- nrow((tail(upsliks,1)[[1]])$logLs)
              if (curr_size < workflow_design$final_reft_size) { # older 42K -> 44K
                # 'reftable_sizes', as previously set by the initialization code, should be OK 
              } else {
                more_reftable_sizes <- tail(reftable_sizes,1)+ seq_len(more_refines)*workflow_design$refine_blocksize
                names(more_reftable_sizes) <- paste0(more_reftable_sizes/1000L,"K")
                reftable_sizes <- c(reftable_sizes, more_reftable_sizes)
              }
            }
            
            ## refine pre-existing upsliks: prepare slots in the upsliks list
            slik_siz <- upsliks[[length(upsliks)]]
            ## Fix: Original script produced reference tables of 42K samples,
            ## extended to 44K in a second step. 
            ## Script was then tidied to generate 44K tables in a single step. 
            ## But next line was wrong, as unsaved reftable sizes became new sizes:
            # newsizes <- setdiff(names(reftable_sizes), c("2K","14K", names(upsliks)))
            ## corrected code: 
            max_in <- max(which(names(reftable_sizes) %in%  names(upsliks)))
            newsizes <- reftable_sizes[-seq_len(max_in)]
            ## /Fix
            supp <- vector("list", length(newsizes))
            names(supp) <- newsizes
            upsliks <- c(upsliks, supp) # add empty elements to be filled by next loop
            
            
            it <- 0L
            for (newsize in newsizes) { # loop to get fits beyond the reftable_sizes[1] one
              it <- it+1L
              locverboses <- verboses
              locverboses$proj <- locverboses$movie <- FALSE
              message(paste0("refine for new size ", newsize,":\n"))
              slik_siz <- refine(slik_siz,  
                                 workflow_design=workflow_design,
                                 # control.Simulate=control.Simulate,
                                 verbose=locverboses, 
                                 # map.asp=1, # tempo fix spaMM <= 4.4.0
                                 CIs=FALSE, eval_RMSEs=FALSE, 
                                 cluster_args=cluster_args, precision=0) #
              if (newsize %in% names(upsliks)) upsliks[[newsize]] <- slik_siz
            }
            if (save_upsliks) save(upsliks, file=paste0("saved_fits/upsliks_",ii,".rda"))

          } 
          ##################################################################
          
        }
        nsim <- if (cmdl_fast_reClu) {
          0L
        } else if (length(LOWER)>8L) {
          if (interactive()) {0L} else {199L}
        } else if (final_MAF) {
          0L # 199L
        } else if (interactive()) {19L} else {199L}
        if (exists("S_obs_table_info") && ! is.null(S_obs_table_info$par.grid)) {
          h0s_ii <- S_obs_table_info$par.grid[ii,]
          h0s_ii <- as.list(h0s_ii)
          for (it in seq_along(h0s_ii)) {names(h0s_ii[[it]]) <- names(h0s_ii)[it]}
        } else h0s_ii <- h0s
        
        max_table_size <- tail(names(upsliks),1L)
        resu <- summary_inference_K( # fn def'd in performance_fns.R 
          slik_siz=upsliks[[max_table_size]], 
          # nsim  for bootstrap LRT: 
          nsim=cmdl_boot_nsim, # in interactive session small bootstrap still useful to catch bugs
          # Note two default #cores controls.
          nb_cores_LRboot= if (interactive()) {
            if (nsim>19L) {min(procs_per_job, 50L %/% length(LOWER))} else {1L}
          } else {min(procs_per_job, 50L %/% length(LOWER))}, # 50G 7 params 9 cores -> several failures
          # on win, very roughly for N_7from17, time = C+2.6 B/C for B replicates, C cores ;
          # => time minimized for C= sqrt(2.6 * B)
          # => In this and higher dimensional pbs, in practice this means 
          #   using all available processors except for toy (low B) computations
          h0s=h0s_ii, ret_HPDint=TRUE,
          abcrf_proj_K=NULL, 
          verbose=(interactive() || on_genotoul),
          reftable_abcrf=NULL,
          rstar_R=0L,
          ii =ii
        )
        save(resu, file=paste0("resu_",ii,".rda"))
      })
  }
}

if (FALSE) { # Reformatting and file manipulation functions
  df_extract <- function(df, what, finalFUN=NULL) {
    for (st in what) {
      df <- df[,st]
      if (is.list(df)) df <- do.call(rbind,df)
    }
    df
  }
  
  merge_jobs <- function(upto, from=1L) {
    seqft <- seq(from=from,to=upto)
    merger <- vector("list",length(seqft))
    it <- 0L
    for (jobit in seqft) {
      it <- it+1L
      # jpath <- paste0("job_",jobit,"/")
      summaries_list <- NULL # safety
      try(chk <- load(paste0("resu_",jobit,".rda"))) # loads 'summaries_list'
      if (inherits(chk,"try-error")) {
        merger[it] <- list(NULL)
      } else merger[[it]] <- get(chk[1])
      # print(df_extract(do.call(rbind, merger[[it]]$SLRTs),"basicLRT"))
    }
    names(merger) <- seqft
    invisible(merger)
  }
  
  save_merge_jobs <- function(upto, from=1L) {
    summaries_list <- merge_jobs(upto,from)
    lnames <- names(summaries_list)
    savename <- paste0("summaries_list.",lnames[1],"_",tail(lnames,1), 
                       ".v",packageVersion("Infusion"),".rda")
    save(summaries_list, file=savename)
    invisible(summaries_list)
  }
  
  # resus <- save_merge_jobs(400)
  ## See Rmarkdown file for analysis of results stored in the summaries_list... file

}

