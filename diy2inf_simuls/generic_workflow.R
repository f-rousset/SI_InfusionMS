## This file must be copied in a working directory
## Then, note how 'dir_generic' and 'thisfilepath' are set by the code below.
## you may need to set the urrent R working directory manually.
## Then, the following big block of code should be run (in full) 
## to initialize variables required for a simulation
## in a way controlled by the name of the working directory.

{ #############      General initialization
  {
    library("Infusion")
    if (packageVersion("Infusion")<"2.1.203") stop("Version 2.1.203 or higher of Infusion is required.")
    if (packageVersion("spaMM")<"4.4.5") warning("Version 4.4.5 or higher of spaMM is strongly advised.")
    # Also required 
    #   Rmixmod (only 'suggested' by Infusion but effectively required in all cases),
    # and at least in some cases:  
    #   abcrf for what its name implies,
    #   reticulate and mafR for use of MAFs
    #   likelihoodAsy for some speculative bootstrap computations (use Rstar=0 to prevent them), 
    #   ggplot2 for some plots.
    
    sysinfo  <- Sys.info()
    if (on_oban <- (sysinfo["nodename"] == "oban"))
      reticulate::use_condaenv(condaenv="maf-conda", conda="~/miniconda3/bin/conda")
    
    user <- sysinfo["user"]
    if (user=="frousset") {
      if (on_oban) {
        dir_generic <- "/home/frousset/Infusionplus/diy2inf/package/R/"
        Infusion::Infusion.options(is_devel_session=TRUE)
      } else dir_generic <- "/work/user/frousset/Infusionplus/diy2inf/package/R/"
    } else if (.Platform$OS.type=="windows") {
      # user=="francois.rousset"
      dir_generic <- "D:/home/francois/travail/stats/Infusionplus/diy2inf/package/R/"
    } else if (user=="francois") {
      dir_generic <- "/mnt/d/home/francois/travail/stats/Infusionplus/diy2inf/package/R/"
      Infusion::Infusion.options(is_devel_session=TRUE)
    } else {
      print(Sys.info())
      stop("You need to adapt this generic workflow to your machine.")
    }
    
    on_genotoul <- (user=="frousset" && ! on_oban)
    on_slurm  <- (on_genotoul || on_oban) && ! interactive()
    
    library("Infusion")
    if (packageVersion("Infusion")<"2.1.186.2") stop("Version 2.1.186.2 or higher of Infusion is required.")
    if (packageVersion("spaMM")<"4.4.5") warning("Version 4.4.5 or higher of spaMM is strongly advised.")
    # Also required: Rmixmod (only 'suggested' by Infusion)

    {
      # library(diy2inf)
      source(paste0(dir_generic,"safe_simul_utils.R")) 
      source(paste0(dir_generic,"DIYABC2Infusion_def.R")) # Those suitable for diyabc -o, including DIYABC_reftable_wrapper
      source(paste0(dir_generic,"generic_reparams.R"))
      source(paste0(dir_generic,"abcrf_addl_defs.R"))
      source(paste0(dir_generic,"Infusion_addl_defs.R"))
      source(paste0(dir_generic,"performance_fns.R"))
      source(paste0(dir_generic,"my_rstar.R"))
    }
    
    if (on_slurm) {
      (thisfilepath <- paste0(getwd(),"/"))
      cmdl_args <- commandArgs(trailingOnly = TRUE)
      # bash script has '--args jobNbr=$job threadNbr=$threadNbr 
      #      and perhaps refine=TRUE or =2...
      cmdl_args <- eval(parse(text=paste0("list(", paste(cmdl_args,collapse=", "),")")))
      if (is.null(cmdl_refine <- cmdl_args$refine)) cmdl_refine <- FALSE
      if (is.null(cmdl_resummarize <- cmdl_args$resummarize)) cmdl_resummarize <- FALSE
      if (is.null(cmdl_fast_reClu <- cmdl_args$fast)) cmdl_fast_reClu <- FALSE 
      if (is.null(cmdl_reproject <- cmdl_args$reproject)) cmdl_reproject <- FALSE
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
      procs_per_job <- max(14L, mydetectCores())
      jobsubdir <- "./"
    }
    
    modelpath_nickname <- tail(strsplit(thisfilepath,"/")[[1]],1) 
    final_MAF <- length(grep("_final_MAF", modelpath_nickname)) > 0L
    do_sim_seq_logTh4 <- modelpath_nickname %in% c("C_axy_8pars_logTh4","D_axy_8pars_logTh4")
    do_sim_seq_logN1 <- modelpath_nickname %in% c("O_13from17_logN1")
    do_sim_seq_logNbn34 <- modelpath_nickname %in% c("O_13from17_logNbn34")
    
    # GLOBAL CONTROL PARAMETERS, not generic:
    if ( length(grep("Harmonia", thisfilepath))) {
      templateHeaderName <- "headerRF_Jan24.txt" # <- "Original_headerRF.txt"
    } else { 
      if (modelpath_nickname %in% c("B_11pars","replicate_old_B")) {
        templateHeaderName <- "headerRFtemplate_11pars.txt"
      } else templateHeaderName <- "headerRFtemplate_17pars.txt"
    }
    
    # checks
    if ( ! length(dir(pattern = templateHeaderName))) stop(paste0("Template header file '", templateHeaderName,"' not found"))
    datafile <- readLines(templateHeaderName)[[1]]
    if ( ! length(dir(pattern = datafile))) stop(paste0("Data file '", datafile,"' not found"))
    
    exe_path <-  paste0(thisfilepath,jobsubdir)
    # procs_per_job must have been provided by the previous block
    if (.Platform$OS.type=="windows") {
      Sobs_simulator_call <- 'diyabc-RF-windows-v1.1.36-o.exe -p ./ -R "ALL" -r 1 -g 1'
      reft_simulator_call <- paste0('diyabc-RF-windows-v1.1.36-o.exe -p ./ -R "ALL" -t ',procs_per_job,' -m -o') # suitable -g will be automatically  added 
      executable <- grep("diyabc-RF",  strsplit(Sobs_simulator_call," ")[[1]],value = TRUE)
      if (dir(thisfilepath, executable) != executable) stop("executable not found in working directory")
    } else {
      Sobs_simulator_call <- './diyabc-RF-linux-v1.1.36-o -p ./ -R "ALL" -r 1 -g 1'
      reft_simulator_call <- paste0('./diyabc-RF-linux-v1.1.36-o -p ./ -R "ALL" -t ',procs_per_job,' -m -o')
      executable <- grep("diyabc-RF",  strsplit(Sobs_simulator_call," ")[[1]],value = TRUE)
      executable <- sub("./", "", x=executable,fixed=TRUE)
      if (! length(dir(thisfilepath, executable))) stop("executable not found in working directory")
    }
    
    control.Simulate <- list(                                
      parsTable="essaisimpars.txt",
      exe_path=exe_path,
      inputHeader=templateHeaderName, 
      outputHeader="headerRF.txt",
      simulator_call=reft_simulator_call,
      verbose.=FALSE, # Whether to suppress diyabc stdout and some other info
      tidy=TRUE # Whether to remove temporary paths and files
    )
    
    reparamHeaderName <- sub(".txt", ".reparam.txt", templateHeaderName, fixed=TRUE)
    if (WITH_REPARAM <- length(dir(pattern=reparamHeaderName))) {
      message("reparam header FOUND.")
      control.Simulate$reparamHeaderName <- reparamHeaderName
      trace(parvec2DIYheader, print=TRUE, tracer=quote(print(inputHeader)))
      trace(ranges2DIYheader, print=TRUE, tracer=quote(print(inputHeader)))
      # LOWER and UPPER may contain composite parameters, which should be handled by reparametrize.
      # reparamHeader specifies times as with t2md3 etc. It could probably use some + and -
      # but was originally developed because of doubts about correct implementation 
      # of these operations in diyabc. So it tends to contain eg a singel variable t2md3
      # rather than t2-d3 or t1+t12-d3...
      # templateHeader is a header that may contain + and -, but no DRAW UNTIL. 
      # templateHeader will not be used by diyabc.
      # Instead my R code must contain a reparametrize() function that is able to convert 
      # variables from the variables in templateHeader to the variables in reparamHeader. 
      # reparametrize() may further convert user-level parameters 
      # (with eg composite parameters, or log transformations...) to the variables in reparamHeader.
      # the ranges in templateHeader are important (together with constr_crits) 
      # to control the actual sampling. The ranges in reparamHeader do NOT control the sampling BUT
      # must be wide enough to contain all (meaningful) values produced by Infusion's parameter sampling.
      # (in the diyabc call on a simulpars table, diyabc checks that the vectors in the values in the table
      # are in the ranges given by reparamHeader => it's still useful to give effective 
      # bounds rather than effectively -Inf,Inf, so that the simulpars table 
      # is effectively checked for nonsense).
    } else {
      message("*NO* reparam header found.")
      try(untrace(parvec2DIYheader), silent=TRUE)
      try(untrace(ranges2DIYheader), silent=TRUE)
    }
    
    if (on_slurm) {
      exe_copied <- FALSE
      if ( ! dir.exists(jobsubdir)) {
        dir.create(jobsubdir)
        file.copy(from=templateHeaderName, to =jobsubdir)
        if (WITH_REPARAM) file.copy(from=reparamHeaderName, to =jobsubdir)
        file.copy(from=datafile, to =jobsubdir)
        # I generate a template S_obs in all cases currently, so I need to copy this:
        exe_copied <- file.copy(from=executable, to =jobsubdir)
        file.copy(from="RNG_state_0000.bin", to =jobsubdir)
        # I generate the S_obs_table whenever absent currently (____F I X M E___), so I need to copy this in all cases:
        # if ( ! (cmdl_resummarize)) {
        if (file.exists("S_obs_table.rda")) {
          file.copy(from="S_obs_table.rda", to =jobsubdir)
        } else stop("S_obs_table.rda is missing. Provide it, **remove subdirectories**, and start again.")
        # }
      }
      setwd(jobsubdir)
      if (exe_copied) Sys.chmod(executable, "0777")
      ## Following tests do not seem to perform as expected. Is doc correct?
      # if ( ! file.access(executable,1)) stop("diyabc executable does not have execute permission.")
      # if ( ! file.access(executable,2)) stop("diyabc executable does not have write permission.")
      if ( ! dir.exists("saved_fits")) dir.create("saved_fits") ## in the jobsubdir
      
      if (modelpath_nickname=="D_smallData_30K") { # recycle earlier simulations
        prevfile <- paste0("../../B_smallData/",jobsubdir,"saved_fits/upsliks_",cmdl_jobNbr,".rda")
        chk <- file.copy(from=prevfile, to="saved_fits/")
        cat(crayon::yellow("Copying upsliks file from B_smallData simulation...", chk))
      }
      
    } else {
      if ( ! file.exists("saved_fits")) cat(crayon::underline("'saved_fits/' subdirectory is missing."))
    }
    
    # ! local defs:
    # (savedir <- "/mnt/c/home/francois/travail/stats/Infusionplus/caseStudies/Harmonia/results/") # WSL
    # (savedir <- "C:/home/francois/travail/stats/Infusionplus/caseStudies/Harmonia/results/")
    # ! setwd() to source file location so that next line finds the file:
    (savedir <- thisfilepath) # this makes it the working directory...     with final '/'
    
    
    
    cat(crayon::yellow("Working directory is currently "),getwd(),"\n")
    if ( ! file.exists("saved_fits")) cat(crayon::underline("'saved_fits/' subdirectory is missing."))
  }
  
  
  { # This block is mandatory for finalization onesim_c() and for its final assignments to control.Simulate 
    
    constr_crits <- NULL # may be overwritten below
    if ( length(grep("Harmonia", thisfilepath))) { 
      # The canonize() definition in this block overrides the one from files sourced 
      # in the 'dir_generic' directory.

      (headerPars <- header2pars(templateHeaderName, actual_mut_pars=c("MEANMU","MEANP")))  # from DIYABC2Infusion_def.R : guesses the 'canonical' parameter names
      par_ranges <- get_par_ranges(headerFile=templateHeaderName, format="matrix") #  from DIYABC2Infusion_def.R
      par_ranges <- rbind(par_ranges,
                          get_mut_par_ranges(headerFile=templateHeaderName, 
                                             actual_mut_pars=c("MEANMU","MEANP"), format="matrix"))
      (headerPars <- names(which(par_ranges[,1]!=par_ranges[,2])))
      
      canonize <- function(parvec) {
        logMu <- parvec[["logMu"]]
        parnames <- names(parvec)
        logThposS <- grep("^logTh\\d+$", parnames)
        if (nlogTh <- length(logThposS)) {
          Nlist <- list()
          splits <- strsplit(parnames[logThposS],"logTh")
          pop_inds <- sapply(splits,`[[`,2)
          for (pop_ind in pop_inds) Nlist[[paste0("N",pop_ind)]] <- 10^(parvec[[paste0("logTh",pop_ind)]] - logMu)
        }
        TposS <- grep("^T\\d+$", parnames)
        parvec <- c(unlist(Nlist),
                    parvec[TposS],
                    ar=parvec[["ar"]],
                    MEANMU=10^logMu,
                    MEANP=parvec[["MEANP"]])
        parvec
      }
      
      # determine LOWER, UPPER : user composite parameters and scaling.
      # For the Human admixture scenario this is automated through the reparameterize definition. 
      # But not here.
      scaLOWER <- c(logTh1=-2,logTh2=-2,logTh3=-2,logTh4=-2, T1=3, ar=0.05,logMu=-5, MEANP=0.01)
      scaUPPER <- c(logTh1=1,logTh2=1,logTh3=1,logTh4=1, T1=45, ar=0.95,logMu=-2,  MEANP=0.5)
      
      onesim_c <- get_onesim_c(names(scaLOWER))
    } else { # GENERIC code: header -> caLOWER -> scaLOWER -> onesim_c formal arguments
      
      ###############
      
      (headerPars <- header2pars(templateHeaderName))  # from DIYABC2Infusion_def.R : guesses the 'canonical' parameter names
      par_ranges <- get_par_ranges(headerFile=templateHeaderName, format="matrix") #  from DIYABC2Infusion_def.R
      (headerPars <- names(which(par_ranges[,1]!=par_ranges[,2])))
      
      canonize <- get_canonizefn(headerPars) # from generic_reparams.R
      
      if (WITH_REPARAM) { # always the case in these simulations.
        reparamNames <- header2pars(reparamHeaderName)
        reparam_ranges <- get_par_ranges(headerFile=reparamHeaderName, format="matrix") #  from DIYABC2Infusion_def.R
        # the variable parameters in template reparamHeader must be all those that are affected by the variable 
        # parameters in the basic template header. There may be more reparam names than headerPars.
        # But this also means that  reparametrize() must have info about fixed pars needed to compute the transformations 
        (var_reparamNames <- names(which(reparam_ranges[,1]!=reparam_ranges[,2])))
        control.Simulate$reparamNames <- var_reparamNames
        control.Simulate$fixedPars <- fixedPars <- par_ranges[,1][which(par_ranges[,1]==par_ranges[,2])]
        chkfixed <- par_ranges[,1]
        chkfixed[] <- NA
        chkfixed[names(fixedPars)] <- fixedPars
        chkfixed1 <- na.omit(reparametrize(chkfixed, to=rownames(reparam_ranges)))
        chkfixed2 <- reparam_ranges[,1][which(reparam_ranges[,1]==reparam_ranges[,2])]
        if ( ! identical(names(chkfixed1), names(chkfixed2))) {
          print(chkfixed1)
          print(chkfixed2)
          stop("Fixed pars in reparam-Header not consistent with fixed pars in basic header")
        } else if ( any(chkfixed1[names(chkfixed2)]!=chkfixed2)) {
          print(chkfixed1)
          print(chkfixed2)
          stop("Fixed pars in reparam-Header not consistent with fixed pars in basic header")
        } 
        # # shared parameters in the two parametrizations must be either variable in both or fixed in both:
        if (FALSE) {
          shared_pars <- intersect(rownames(par_ranges),rownames(reparam_ranges))
          fixedinHeaderamongsharedPars <- par_ranges[shared_pars,1][which(par_ranges[shared_pars,1]==par_ranges[shared_pars,2])]
          fixedinreparamHeaderamongsharedPars <- reparam_ranges[shared_pars,1][which(reparam_ranges[shared_pars,1]==reparam_ranges[shared_pars,2])]
          if ( ! setequal(names(fixedinHeaderamongsharedPars),names(fixedinreparamHeaderamongsharedPars))) {
            stop("shared parameters differ by which are fixed and which are variable in the two parametrizations.")
          }
          if ( any(fixedinreparamHeaderamongsharedPars[names(fixedinHeaderamongsharedPars)]!=fixedinHeaderamongsharedPars)) {
            stop("Fixed shared parameters in the two parametrizations differ by their value.")
          }
        }
      } else {
        reparamNames <- reparam_ranges <- NULL
      }
      
      
      # create ranges in canonical parameter space
      caLOWER <- numeric(length(headerPars)) # headerPars from templateHeaderName (: the not reparam one)
      names(caLOWER) <- headerPars
      caUPPER <- caLOWER
      caLOWER[headerPars] <- par_ranges[headerPars,1]
      caUPPER[headerPars] <- par_ranges[headerPars,2]
      
      { # composite parameters
        if (length(grep("O_13from17", modelpath_nickname)) || 
            modelpath_nickname %in% c("P_17from17","O_14from14","O_12from17","O_12from12")) {
          posd3 <- names(caLOWER)=="d3"
          caLOWER[posd3] <- 0
          names(caLOWER)[posd3] <- "d3_Nbn3"
          caUPPER[posd3] <- 1
          names(caUPPER)[posd3] <- "d3_Nbn3"
          
          posd4 <- names(caLOWER)=="d4"
          caLOWER[posd4] <- 0
          names(caLOWER)[posd4] <- "d4_Nbn4"
          caUPPER[posd4] <- 1
          names(caUPPER)[posd4] <- "d4_Nbn4"
          
          posd34 <- names(caLOWER)=="d34"
          caLOWER[posd34] <- 0
          names(caLOWER)[posd34] <- "d34_Nbn34"
          caUPPER[posd34] <- 1
          names(caUPPER)[posd34] <- "d34_Nbn34"
        }
        if (modelpath_nickname %in% c("P_4from17")) {
          posd3 <- names(caLOWER)=="d3"
          caLOWER[posd3] <- 0
          names(caLOWER)[posd3] <- "d3_Nbn3"
          caUPPER[posd3] <- 1
          names(caUPPER)[posd3] <- "d3_Nbn3"
        }
        if (modelpath_nickname %in% c("P_17from17","O_14from14")) {
          posNa <- names(caLOWER)=="Na"
          caLOWER[posNa] <- 0.05
          names(caLOWER)[posNa] <- "Na_N2"
          caUPPER[posNa] <- 1
          names(caUPPER)[posNa] <- "Na_N2"
        }
        
      }  
      
      # create ranges and DGpars in rescaled (eg, log) canonical parameter space
      (tmp <- rescaleLoUp(caLOWER,caUPPER, fac=NA)) # from generic_reparams.R
      
      scaLOWER <- tmp[["scaLOWER"]]
      scaUPPER <- tmp[["scaUPPER"]]
      
      if (modelpath_nickname %in% c("P_17from17","O_14from14")) {
        # constr crits are in the user-level coordinates, in the scaLOWER... ;
        constr_crits <- quote({c(
          10^log.Nbn34.* d34_Nbn34 - 10^log1p.t23.-1,
          10^log.Nbn4.* d4_Nbn4-(10^log1p.t12.-1), # t2-d4 must be > t1 ie d4 < t12
          10^log.Nbn3.* d3_Nbn3-(t1+10^log1p.t12.-1) # t2-d3 must be >0 
        )})
      } else if (modelpath_nickname %in% c("O_13from17","O_13from17_logN1","O_13from17_logNbn34",
                                           "O_12from17")) {
        # constr crits are in the user-level coordinates, in the scaLOWER... ;
        # as.integer() should be applied to values matching the internal times in diyabc 
        # But it's not clear when does diyabc perform integer conversions...
        # There seem to be rather late (int)(0.5+... ) conversions eg particuleC.cpp l.637
        # But I found that diyabc fails on (integers) t2=t2-d3 ie d3=0L fails
        # That means there is an effectively integer t2-d3 rather than simply integers t2, d3
        constr_crits <- quote({
          c(
            as.integer(10^log.Nbn34.* d34_Nbn34- 10^log1p.t23.-1),
            as.integer(10^log.Nbn4.* d4_Nbn4-10^log.t12.), # t2-d4 must be > t1 ie d4 < t12 
            as.integer(10^log.Nbn3.* d3_Nbn3-(t1+10^log.t12.)), # t2-d3 must be >0 
            as.integer((t1+10^log.t12.+10^log1p.t23.-1) - 1006.936), # t3 - [fixed t4]
            - as.integer(10^log.Nbn3.* d3_Nbn3), # I found that diyabc fails on cases where (integers) t2=t2-d3 ie d3=0L fails
            - as.integer(10^log.Nbn4.* d4_Nbn4), # Cautiously applying to other bottleneck durations
            - as.integer(10^log.Nbn34.* d34_Nbn34) # Cautiously applying to other bottleneck durations
          )})
      } else if (modelpath_nickname %in% c("P_4from17")) {
        # t1 t23 d34_Nbn34 log.Nbn4. d4_Nbn4 and log.Nbn34. are now fixed
        constr_crits <- quote({
          c(
            as.integer(19.95262-10^log.t12.), # t2-d4 must be > t1 ie d4 < t12 
            as.integer(10^log.Nbn3.* d3_Nbn3-(6+10^log.t12.)), # t2-d3 must be >0 
            as.integer((6+10^log.t12.+18.95262) - 1006.936), # t3 - [fixed t4]
            - as.integer(10^log.Nbn3.* d3_Nbn3) # I found that diyabc fails on cases where (integers) t2=t2-d3 ie d3=0L fails
          )})
      }
        
      { ## check transfos
        med <- (scaLOWER+scaUPPER)/2
        chk <- reparametrize(canonize(med), to=names(scaLOWER)) -med
        if (max(abs(chk))>1e-6) stop("Problem with parameter transformations.")
        
      }
      
      onesim_c <- get_onesim_c(names(scaLOWER)) # get_onesim_c() from DIYABC2Infusion_def.R 
      
    } 
    
    { # the $caParNames must be the names uses by diyabc in the order assumed by diyabc
      # ideally this info is deduced from the header file, but my generic code does not yet handles all forms of header files
      # so FIX: ise the result of canonize (using itself-fixed local canonize so that it returns params in the order expected by diyabc)
      (caLOWER <- canonize(scaLOWER)) 
      control.Simulate$caParNames <- names(caLOWER) # allow check of colnames of canonical grid. (=> no repetitive message) ( ! proper ordering)
      caUPPER <- canonize(scaUPPER)
    }
    
    
    
    (LOWER <- scaLOWER)  # They have to be the same thing, with different names for different contexts
    UPPER <- scaUPPER
    
  }
  
  {  # This block is mandatory for assignment of control.Simulate$statNames. It also provides DGPars and reftable sizes
    # It contains the fit-real-data subblock
    { 
      (npar <- length(LOWER))
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
    
    
    
    { #I will use names ordered as in simulated Sobs (does it matter?). For that we still need a template Sobs
      rm(list = intersect(ls(), c("template_Sobs")))
      cat(crayon::yellow("\nGenerating template 'Sobs'\n"))
      template_pars <- (scaLOWER+scaUPPER)/2
      if (length(grep("O_13from17", modelpath_nickname)) ||
                 modelpath_nickname %in%  c("P_17from17","O_14from14","O_12from17","O_12from12")) {
        template_pars["log1p.t23."] <- log10(1+50)
        template_pars["d3_Nbn3"] <- 0.1
        template_pars["d4_Nbn4"] <- 0.1
        template_pars["d34_Nbn34"] <- 0.1
      } else if (modelpath_nickname %in%  c("P_4from17")) {
        template_pars["d3_Nbn3"] <- 0.1
      }
      if (! is.null("constr_crits")) {
        chk <- eval(constr_crits,envir = as.list(template_pars))
        if (any(chk>=0)) stop("'template_pars' does not obey 'constr_crits' constraints.")
      }
      set_seeds(123, simulator_call=control.Simulate$simulator_call)
      (template_Sobs <- do.call(onesim_c,
                                c(as.list(template_pars),
                                  list(control.Simulate=control.Simulate,
                                       simulator_call=Sobs_simulator_call) )))
      control.Simulate$statNames <- names(template_Sobs)
      # Turn of the traces once 'inputHeader' name has been checked
      try(untrace(parvec2DIYheader), silent=TRUE)
      try(untrace(ranges2DIYheader), silent=TRUE)
    }
    
    { 
      statobsRFname <-"statobsRF.txt" # diyabc tends to overwrite it?
      chk <- try(actual_obs <- readLines(statobsRFname)) 
      if (inherits(chk, "try-error")) {
        message("no statobsRF.txt file: no check of order of statistics in simulation output.")
      } else {
        actual_obs <- strsplit(actual_obs,"[ ]+")
        statnames <- actual_obs[[1]][-1] # I have used the more complicated expression setdiff(strsplit(actual_obs[1]," ")[[1]],"")
        statnames <- gsub("&",".",statnames) # rename varas when they contain &
        if ( ! all(statnames==names(template_Sobs))) {# => only to see whether there are permutations 
          warning(paste0("Order of statistics in ",statobsRFname," differs form that of template S_obs."))
        } 
        actual_obs <- as.numeric(actual_obs[[3]][-1]) #I have used ... na.omit(as.numeric(strsplit(statobs[3]," ")[[1]]))
        names(actual_obs) <- statnames
        actual_obs <- actual_obs[names(template_Sobs)] # same order as in my workflow
      }
    }  

    # scaDGP and  scaDGparlist:
    fit_real_data <- FALSE # Change this to refit them.
    if (fit_real_data) { # fit real data to obtain the DGpars 'scaDGP'
      set_seeds(123,simulator_call = control.Simulate$simulator_call)
      # debug(calc_sliks) # handy
      Rprof(filename="Rprof.fit_real_data.out", 
            event="elapsed", # non-default under linux, necessary to count time spect in system() calls 
            interval=(length(LOWER)/9)^3) 
      fits <- get_upsliks(ii=NA,  # fn def'd in Infusion_addl_defs.R
                          S_obs = actual_obs,
                          reftable_sizes=reftable_sizes,
                          control.Simulate=control.Simulate,
                          LOWER=LOWER, UPPER=UPPER,
                          cluster_args=list(project=list(num.threads=procs_per_job)),
                          methodArgs=list(importance="none"), # FIXME
                          verbose= list(most=(interactive() || on_slurm),
                                        rparam=(interactive() || on_slurm)),
                          save_upsliks=FALSE,
                          scaDGparlist=NULL,
                          more_refines=FALSE,
                          workflow_design=workflow_design,
                          constr_crits=constr_crits
      )
      Rprof(NULL)
      (max_table_size <- tail(names(fits),1)) #  "14K", say: simple nickname, not numeric value.
      summ_inference_maxK <- summary_inference_K(
        slik_siz=fits[[max_table_size]], S_obs=actual_obs, nsim=0L, # fn def'd in performance_fns.R
        nb_cores_LRboot=1L, ncores_abcrf=procs_per_job,
        h0s=NULL, ret_HPDint=TRUE,
        abcrf_proj_K=NULL, 
        ii="NA",  
        reftable_abcrf=reftable_abcrf_final[seq(na.omit(reftable_sizes[max_table_size])),])
      save(fits, summ_inference_maxK, file=paste0("actual_data_inferences.v",packageVersion("Infusion"),".rda"))
      if (on_slurm) {
        stop("stop() at end of block for analysis of actual data.")
      } else {
        cat(crayon::yellow("browser() at end of block for analysis of actual data."))
        browser()
      }
      
      if ( length(grep("Harmonia", thisfilepath))) { 
        tmpDGP <- fits[[max_table_size]]$MSL$MSLE
        (scaDGP <-  signif(tmpDGP,2))
        scaRANGE <- scaUPPER-scaLOWER
        any(chk_at_bound <- c((scaDGP-scaLOWER)/scaRANGE,
                              (scaUPPER-scaDGP)/scaRANGE)<0.049999999)
        if (modelpath_nickname =="B_axy_8pars") scaDGP["logTh4"] <- 0.85
      } else if ( length(grep("admixtOutOfA", thisfilepath)) ) {
        # Use ML estimates from the likelihood profile for t1=6  
        if (modelpath_nickname %in% c("N_7from17")) {
          load("actual_data_inferences.v2.1.79.2.rda")
          tmpDGP <- fits[[length(fits)]]$MSL$MSLE
          tmpDGP[["t1"]] <- 6 # flat profile
          (scaFixeds <- c(t1=6))
          prof_info <- profile(fits[[length(fits)]], scaFixeds) 
          prof_sol <- attr(prof_info,"solution")
          tmpDGP[names(prof_sol)] <- prof_sol
          (scaDGP <-  signif(tmpDGP,2))
          scaRANGE <- scaUPPER-scaLOWER
          any(chk_at_bound <- c((scaDGP-scaLOWER)/scaRANGE,
                                (scaUPPER-scaDGP)/scaRANGE)<0.049999999)
        } else if (modelpath_nickname %in% c("O_13from17")) {
          load("actual_data_inferences.v2.1.192.rda")
          tmpDGP <- fits[[length(fits)]]$MSL$MSLE
          tmpDGP[["t1"]] <- 6 # flat profile
          (scaFixeds <- c(t1=6))
          prof_info <- profile(fits[[length(fits)]], scaFixeds) 
          prof_sol <- attr(prof_info,"solution")
          tmpDGP[names(prof_sol)] <- prof_sol
          (scaDGP <-  signif(tmpDGP,2))
          scaRANGE <- scaUPPER-scaLOWER
          any(chk_at_bound <- c((scaDGP-scaLOWER)/scaRANGE,
                                (scaUPPER-scaDGP)/scaRANGE)<0.049999999)
        } else if (modelpath_nickname %in% c("B_13from17")) {
          # values for this case appear to be modified ones (for debugging purposes) 
          # from those from another fit which was likely:
          if (FALSE) {
            load("../J_11from17/actual_data_inferences.v2.1.79.2.rda")
            tmpDGP <- fits[[length(fits)]]$MSL$MSLE
            tmpDGP[["t1"]] <- 6 # flat profile
            (scaFixeds <- c(t1=6))
            oldfit <- fits[[length(fits)]]
            oldfit = MSL(oldfit) 
            prof_info <- profile(oldfit, scaFixeds) 
            prof_sol <- attr(prof_info,"solution")
            tmpDGP[names(prof_sol)] <- prof_sol
            tmpDGP <- (signif(tmpDGP,2))
            "log.N1. = 3.5, log.N2. = 3.5, log.N3. = 4.2, log.N4. = 3.5, t1 = 6, ra = 0.17, 
                 log.t12. = 2.2,                  log1p.t23. = 1.8, Nbn34 = 80,            
                     log1p.t34. = 3,   log.Na. = 2.8"
            # vs values for B_13from17 
            "log.N1. = 3.2, log.N2. = 3.5, log.N3. = 3.8, log.N4. = 3.5, t1 = 6, ra = 0.17, 
                 log.t12. = 2.2, d3 = 42, d4 = 9, log1p.t23. = 1.7,              d34 = 24, 
                     log1p.t34. = 2.9, log.Na. = 2.7"
          }
        }
      }
      scaRANGE <- scaUPPER-scaLOWER
      chk_at_bound <- c((scaDGP-scaLOWER)/scaRANGE,
                        (scaUPPER-scaDGP)/scaRANGE)<0.0499999999
      if (any(chk_at_bound)) {
        warning("DG parameters close to bounds or out of range:", immediate. = TRUE)
        print(rep(scaDGP,2)[which(chk_at_bound)])
      }
      paste0(names(scaDGP),"=", scaDGP, collapse=", ") # for easy copy
    }  
    
    # Data-generating values used in performance simulations, often derived by the above fits:
    if ( ! fit_real_data) {
      if ( length(grep("Harmonia", thisfilepath)) ) { # Ladybird invasion scenario
        if ( 
          length(grep("D_axy_8pars", modelpath_nickname)) # grep() to match "D_axy_8pars_logTh4" too
        ) { # values of original 'A_2023_10_8pars' simuls, differents from the new A_2023_10_8pars below
          scaDGP <- c(logTh1=0, logTh2=0, logTh3=0, logTh4=0, T1=13, 
                      ar=0.5, logMu=-3, MEANP=0.25) 
        } else if (modelpath_nickname== "A_2023_10_8pars") { # nickname was reused for different parameters
          # Comments= results from this new A_2023_10_8pars => reason for change in C_axy_8pars
          scaDGP <- c(logTh1=-0.031, # ML ~posterior. No modif                                   
                      logTh2=-0.84, # Bias opposite from prior mean: confirm further from it      
                      logTh3=0.56, # deficiency of high values                                      
                      logTh4=-0.88, # indep of ML and seems superefficient: test on 0 (other side of prior mean)                  
                      T1=19, # No info: push out of central zone                                  
                      ar=0.5, # deficiency of low values                                         
                      logMu=-2.7, # Post ~ML with negative bias. No obvious thing to test        
                      MEANP=0.22) # ML ~posterior. No modif                                       
        } else if ( 
          length(grep("C_axy_8pars", modelpath_nickname)) # grep() to match "C_axy_8pars_logTh4" too
        ) { # as modified from new A_2023_10_8pars
          #                         conclusion   
          scaDGP <- c(logTh1=-0.031, # RAS
                      logTh2=-1.5, # small bias away. Not graphical. weird
                      logTh3=0.75, # strong bias. Effect of prior  
                      logTh4=0, # small bias away from prior. Low variance => low RMSE !
                      T1=9, # in fact this shows there is some info and expected bias of posterior ests 
                      ar=0.2, # confirmed. Commonplace upper bias of posterior ests
                      logMu=-2.7, # Better seen 'commonplace bias'. *Why the defic. of low p-values?*
                      MEANP=0.22 # RAS
          ) 
        }
      }  else if ( length(grep("admixtOutOfA", thisfilepath)) ) { # Human admixture scenario
        scaDGP <- switch(
          modelpath_nickname,
          "N_7from17" = { 
            c(log.N2.=3.5, t1=6,   log.t12.=2.2, 
              log1p.t23.=1.5,         Nbn34=65, log1p.t34.=3, log.Na.=2.7)
          },
          "N_7from17_final_MAF" = { # same values as N_7from17  
            c(log.N2.=3.5, t1=6,   log.t12.=2.2, 
              log1p.t23.=1.5,         Nbn34=65, log1p.t34.=3, log.Na.=2.7)
          },
          "Q_7from17" = { # same values as N_7from17, for simulations with 10000 SNPs
            c(log.N2.=3.5, t1=6,   log.t12.=2.2, 
              log1p.t23.=1.5,         Nbn34=65, log1p.t34.=3, log.Na.=2.7)
          },
          "N_7from17_2_final_MAF" = { # same values as N_7from17  
            c(log.N2.=3.5, t1=6,   log.t12.=2.2, 
              log1p.t23.=1.5,         Nbn34=65, log1p.t34.=3, log.Na.=2.7)
          },
          "Q_7from17" = { # = N_... but with 10000 SNPs
            c(log.N2.=3.5, t1=6,   log.t12.=2.2, 
              log1p.t23.=1.5,         Nbn34=65, log1p.t34.=3, log.Na.=2.7)
          },
          "B_13from17" = { 
            c(log.N1.=3.2, log.N2.=3.5, log.N3.=3.8, log.N4.=3.5, t1=6, ra=0.17, 
              log.t12.=2.2, d3=42,          d4=9,          
              log1p.t23.=1.7, d34=24, log1p.t34.=2.9, log.Na.=2.7) 
          },
          "O_13from17" = { # reparametrization with dx_Nbnx's
            c(log.N1.=4.7, log.N2.=3.4, t1=6, ra=0.21, log.t12.=2.2, d3_Nbn3=0.26, 
              log.Nbn3.=2.1, d4_Nbn4=0.1, log.Nbn4.=2.3, log1p.t23.=1.3, 
              d34_Nbn34=0.44, log.Nbn34.=1.5, log.Na.=2.7)
          },
          "O_13from17_logN1" = { # same as O_13from17 but logN1 will be modified
            c(log.N1.=4.7, log.N2.=3.4, t1=6, ra=0.21, log.t12.=2.2, d3_Nbn3=0.26, 
              log.Nbn3.=2.1, d4_Nbn4=0.1, log.Nbn4.=2.3, log1p.t23.=1.3, 
              d34_Nbn34=0.44, log.Nbn34.=1.5, log.Na.=2.7)
          },
          "O_13from17_logNbn34" = { # same as O_13from17 but logNbn34 will be modified
            c(log.N1.=4.7, log.N2.=3.4, t1=6, ra=0.21, log.t12.=2.2, d3_Nbn3=0.26, 
              log.Nbn3.=2.1, d4_Nbn4=0.1, log.Nbn4.=2.3, log1p.t23.=1.3, 
              d34_Nbn34=0.44, log.Nbn34.=1.5, log.Na.=2.7)
          },
          "P_4from17" = { # fixing parameters from O_13. simulated data are the same (only the saves4Rmd was regenerated)
            c(                                         log.t12.=2.2, d3_Nbn3=0.26, 
                                                       log.Nbn3.=2.1,  
                                                       log.Na.=2.7)
          },
          stop("'scaDGP' not predefined for this model")
        )
        
      }
    }
    
    # For SLRTs, whether 'fit_real_data' is TRUE or not or not:
    scaDGparlist <- as.list(scaDGP)
    h0s <- scaDGparlist; for (it in seq_along(h0s)) {names(h0s[[it]]) <- names(h0s)[it]}
  }
  { # info:
    fixedParsNames <- names(which(par_ranges[,1]==par_ranges[,2]))
    cat("Fixed by the template file:\n")
    print(par_ranges[fixedParsNames,1])
  }

}


#################### Precomputations for abcrf workflow;
# Building the objects takes time, and even simplyloading them may take some time
# Build/load abcrf reftable and provide rf regressions 
# predict.regAbcrf will need both the full reftable and the projections 
#  (rf regressions) as distinct objects, so both must be loaded in memory.  

# reasonable guess that I try to generate the table on my linux machine:
if ( ! on_genotoul &&
     is.null(cmdl_jobNbr) && # If not in series of analyses on cluster 
     ! .Platform$OS.type=="windows" &&
     ! length(grep("_final_MAF", modelpath_nickname)) &&
     length(scaDGP)<14L) { 
  
  ## First check whether the reference table has not already been computed
  (chk <- safe_load(c("reftable_abcrf"), strict=FALSE))
  
  ## Compute it if missing
  if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
    cat(crayon::yellow("Building the abcrf reftable\n"))
    
    set_seeds(456,simulator_call = control.Simulate$simulator_call)
    parsp <- init_reftable(lower=LOWER, upper=UPPER, nUnique = default_indiv_reft_size,
                           constr_crits = constr_crits)
    parsp <- parsp[sample(nrow(parsp)), ] # random ordering useful for debugging at least
    
    reftable_abcrf_final <- add_reftable(    Simulate="DIYABC_reftable_wrapper",
                                             par.grid=parsp,
                                             control.Simulate=control.Simulate
    )
    save(reftable_abcrf_final, file="reftable_abcrf.rda")
  } else cat(crayon::yellow("'reftable_abcrf' file found\n"))
  
  ## Same for the RF regressions: first check whether they have already been computed
  if (exists("abcrf_proj_all_K")) {
    cat(crayon::yellow("'abcrf_proj_all_K' already loaded (BUT IT SHOULD BE REGENERATED if simulation design has been modified)\n"))
  } else {
    cat(crayon::yellow("Loading the abcrf projections may be slow\n"))
    chk <- safe_load(c("abcrf_projs"), strict=FALSE)
    
    if (inherits(chk,"try-error")) {
      cat(crayon::yellow("Building the abcrf projections\n"))
      abcrf_proj_all_K <- list() # empty list: don't assign names to any yet-NULL elements: cf def of 'missing_sizes'
    } else cat(crayon::yellow("'abcrf_projs' file found (may be slow to load)\n"))
  }
  
  ## Identify missing RF regression objects
  if (default_indiv_reft_size>20000L) {
    if ( ! length(abcrf_proj_all_K)) cat(crayon::yellow("Building projections only for largest size, to save a lot of memory:\n"))
    abcrf_reftable_sizes <- tail(reftable_sizes,1L) # to save a lot of memory
  } else abcrf_reftable_sizes <- reftable_sizes
  
  sizes <- names(abcrf_reftable_sizes)
  
  missing_sizes <-  setdiff(sizes, names(abcrf_proj_all_K))
  
  ## Provide the missing ones:
  if (length(missing_sizes)) {
    cat(crayon::yellow(paste("Building the abcrf projections for", paste(missing_sizes, collapse=", "),"\n")))
    library("abcrf")
    ntree <- 1000L  
    min.node.size <- 5L
    for (siz in missing_sizes) {
      abcrf_proj_all_K[[siz]] <- build_abcrf_projections(reftable_abcrf_final[1:reftable_sizes[siz],], # fn def'd in abcrf_addl_defs.R
                                                         parNames=names(LOWER), Sobs=template_Sobs, ntree=ntree, 
                                                         min.node.size=min.node.size, 
                                                         ncores=procs_per_job # was previously set to 10L on francois computer
                                                         # When run on a shared cluster, it is important to control ncores explicitly
                                                         # Otherwise ranger's 'Default is number of CPUs available.
      )
    }
    cat("\n") # ends a line of parameter names
    save(abcrf_proj_all_K, file="abcrf_projs.rda") # FIXME
    
    if (on_slurm) {
      stop("stop() at end of block for computation of 'reftable_abcrf_final' and 'abcrf_proj_all_K'.")
    } else {
      cat(crayon::yellow("browser() at end of block for computation of 'reftable_abcrf_final' and 'abcrf_proj_all_K'."))
      browser()
    }
  }

} # (end) Precomputations for abcrf workflow

if (FALSE) { # O_13from17_logNbn34 and O_13from17_logN1
  source(paste0(dir_generic,"../../generic_workflow/adhoc_run_do_sim_seq.R"))
}

{ # The performance simulations
  
  ## Check directory architecture, and availability of required and optional files.  
  if ( ! file.exists("saved_fits")) stop(paste("No 'saved_fits' directory found in", getwd()))

  if (RUN_ABCRF <- (length(scaDGP)<14L && 
                    ! length(grep("_MAF", modelpath_nickname)) &&
                    ! cmdl_fast_reClu &&
                    ! length(grep("altsampling",thisfilepath)))) {
    library("abcrf")
    
    if (on_slurm && ! is.null(cmdl_jobNbr)) { # (second condition superfluous in routine simulations)
      chk <- try(load("../reftable_abcrf.rda"))
      if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
        stop("reftable_abcrf is missing for slurm job")
      }
      
      chk <- try(load("../abcrf_projs.rda"))
      if (inherits(chk,"try-error")) { 
        stop("abcrf_projs.rda is missing for slurm job") #because parallel execution
      } 
      cat(crayon::yellow(paste0("Precomputed info for 'abcrf' found and loaded.")))
    } else {
      if ( ! exists("abcrf_proj_all_K")) {
        chk <- safe_load(c("abcrf_projs"), strict=FALSE)
        if (inherits(chk,"try-error")) {
          warning("abcrf_proj_all_K not evaluated!", immediate.=TRUE)
          abcrf_proj_all_K <- NULL # so that summary_inference_K() handles it 
        }
      }
      if ( ! exists("reftable_abcrf_final")) {
        (chk <- safe_load(c("reftable_abcrf"), strict=FALSE))
        if (inherits(chk,"try-error")) {
          warning("reftable_abcrf_final not evaluated!", immediate.=TRUE)
          reftable_abcrf_final <- NULL # so that summary_inference_K() handles it 
        }
      }
    } 
  } else  abcrf_proj_all_K <- reftable_abcrf_final <- NULL
  
  { ##### S_obs_table
    if (length(setdiff(names(LOWER), names(scaDGparlist)))) 
      stop("'LOWER' and 'scaDGparlist' elements do not match.")
    rm(list=intersect(ls(), c("S_obs_table")))
    chk <- safe_load(c("S_obs_table"))
    if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
      if (modelpath_nickname == "P_4from17") {
        warning("S_obs_table.rda should be copied from the 'O_13from17' directory")
        browser()
        # below there is some code showing how the original O_13from17 S_obs_table
        # was extended when more smaples were needed for P_4from17.
      }
      cat(crayon::yellow("Building the S_obs_table\n"))
      if (modelpath_nickname %in% c("N_7from17","D_axy_8pars","P_4from17","O_13from17")) {
        NREPL <- 1000L  
      } else NREPL <- 200L 
      set_seeds(567,simulator_call = control.Simulate$simulator_call)
      ###
      # For reasons that belong to the simulator, setting its seed ensures repeatability
      # of the first samples only for constant NREPL. 
      # Thus, some manipulations are needed when extending 
      # a previous S_obs_table by increasing NREPL
      # without modifying the samples from simulations with lower NREPL. 
      # This is relevant for the O_13from17 table (also used by P_4from17) 
      # whose original table of 200 simulations was later extended.
      ###
      par.grid <- data.frame(t(scaDGP)[rep(1,NREPL),])
      S_obs_table_info <- list(WITH_REPARAM=WITH_REPARAM, LOWER=LOWER, UPPER=UPPER,
                               scaDGparlist=scaDGparlist)
      
      if (do_sim_seq_logTh4) { # modify par.grid 
        seq_logTh4 <- seq(LOWER["logTh4"],UPPER["logTh4"],length.out=12)[2:11]
        par.grid$logTh4 <- seq_logTh4[gl(10L,NREPL %/% 10L)]
        S_obs_table_info$par.grid <- par.grid
      } else if (do_sim_seq_logN1) {
        seq_logN1 <- seq(LOWER["log.N1."],UPPER["log.N1."],length.out=12)[2:11]
        par.grid$"log.N1." <- seq_logN1[gl(10L,NREPL %/% 10L)]
        S_obs_table_info$par.grid <- par.grid
      } else if (do_sim_seq_logNbn34) {
        seq_logNbn34 <- seq(LOWER["log.Nbn34."],UPPER["log.Nbn34."],length.out=12)[2:11]
        # Two fixes necessary to respect constraints on parameters:
        par.grid$"log1p.t23." <- 2.5
        seq_logNbn34 <- seq(LOWER["log.Nbn34."],seq_logNbn34[7],length.out=12)[2:11]
        par.grid$"log.Nbn34." <- seq_logNbn34[gl(10L,NREPL %/% 10L)]
        # apply(par.grid,1L, function(v) all(eval(constr_crits, as.list(v))<0))
        S_obs_table_info$par.grid <- par.grid
      }
      
      ## SIMULATION of data sets
      S_obs_table <- add_reftable(Simulate="DIYABC_reftable_wrapper",
                                  par.grid=par.grid,
                                  control.Simulate=control.Simulate)[,-seq(npar)]
      readparms <- as.numeric(strsplit(readLines(control.Simulate$parsTable)[1]," ")[[1]][-1])
      if (WITH_REPARAM) {
        tarfile <- "S_obs_table.header_S_.tar"
        names(readparms) <- var_reparamNames
      } else {
        tarfile <- "S_obs_table.header.tar"
        names(readparms) <- headerPars
      }
      if ( ! (do_sim_seq_logTh4 || do_sim_seq_logN1 || do_sim_seq_logNbn34)) { # check the parameters as they have been written on file read by simulator:
        if (modelpath_nickname == "P_4from17") { # fixed t1 is needed
          chk <- scaDGP - reparametrize(readparms,to=names(scaDGP), fixedPars=fixedPars)
        } else chk <- scaDGP - reparametrize(readparms,to=names(scaDGP))
        if (max(abs(chk))>1e-6) {
          print(chk)
          stop("Parameter values read from file do not seem correct.")
        }
      }
      
      # Save critical files as info about S_obs_table generation.
      tar(tarfile, files=c(control.Simulate$inputHeader,
                           control.Simulate$reparamHeaderName),
          compression="none")
      if (modelpath_nickname=="P_4from17") { # see comments on set_seeds() call above
        new_S_obs_table <- S_obs_table
        new_info <- S_obs_table_info
        load("../O_13from17/S_obs_table.rda") # the first 200 samples originally used for P_4from17
        # This also reads the S_obs_table_info for O_13from17 which replaces the local one.
        # S_obs_table_info$LOWER and $UPPER are confusing, they should not be used
        # The LOWER and UPPER in the .GlobalEnv are always thos deduced from the diyabc headers.
        # => the saves4Rmd.rda elements are not modified and are specific to P_4from17)
        # For safety we remove:
        S_obs_table_info$LOWER <- S_obs_table_info$UPPER <- NULL
        S_obs_table <- rbind(S_obs_table, # first 200 read from O_13from17
                             new_S_obs_table[-(1:200),])
      } 
      save(S_obs_table, S_obs_table_info, file="S_obs_table.rda")
      # Save info used in performance summaries
      if ( ! file.exists("saves4Rmd.rda")) {
        save(LOWER,UPPER,scaDGparlist, file="saves4Rmd.rda")
      } else cat(crayon::yellow("saves4Rmd.rda already exists..."))
    } else {
      cat(crayon::yellow("'S_obs_table' file found; BUT IT SHOULD BE REGENERATED if simulation design has been modified.\n")) # FIXME
      NREPL <- nrow(S_obs_table)
    }
    cat(paste("NREPL=", NREPL, "\n")) # FIXME
    
    if (is.null(control.Simulate$statNames)) stop("'statNames' missing from 'control.simulate': re-run the block which needs 'template_Sobs'.")
    
    cluster_args <- list(project=list(num.threads=procs_per_job)) 
    if (identical(cluster_args$rparam$type,"FORK")) require("progressr") 

    # Default values controlling which simulated data sets will be used. May be manually modified
    if (is.null(cmdl_jobNbr)) { 
      upto <- min(1L, NREPL)
      ii <- prev_latest_one <- 0L 
      # ii <- prev_latest_one <- 175L # continuing simulation: the index of the last simulated 
      cat(crayon::yellow("Starting from the ",prev_latest_one+1L,"th replicated sample.\n"))
    } else {
      upto <- min(cmdl_jobNbr, NREPL)
      ii <- prev_latest_one <- cmdl_jobNbr-1L 
    } 
    
  } #  ! Lengthy inference simulations in the next block ! :
  
  #  ! Lengthy inference simulations in the next block ! :
  
  if (ii < NREPL) { # so that setting prev_latest_one to NREPL prevents the creation of a new summaries_list
    TMP <- summaries_list <- replicate(
      upto-prev_latest_one, {
        ii <<- ii+1L
        cat(crayon::bgYellow$black(ii),"\n")
        S_obs <- unlist(S_obs_table[ii,])
        set_seeds(678L+ii, simulator_call=control.Simulate$simulator_call) #set.seed(678L+ii)
        Infusion.options(is_devel_session= user %in% c("francois.rousset","francois")) # position handling package reinstallations
        if (length(grep("altsampling",thisfilepath))) { # e.g. "[...]/altsampling/N_7from17"
          Infusion.options(rparamFn=Infusion:::.rparam_SLik_j_B_postdens)
          message("***using alternative rparamFn***")
        }
        
        if (cmdl_resummarize || final_MAF || cmdl_fast_reClu) {
          MSL_recomputed <- FALSE
          ## Context for next block:
          # Fit results are saved in the job_XXX/saved_fits/ subdir
          # but are read in the /saved_fits/ subdir
          # so if we forget copying...
          if (exists("upsliks",envir = .GlobalEnv)) {
            cat(crayon::yellow(paste0(" 'upsliks' already exists in .GlobalEnv")))
          } else {
            if (on_genotoul <- (user=="frousset")) {
              (chk <- safe_load(paste0("./../saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
            } else (chk <- safe_load(paste0("./saved_fits/upsliks_",ii), strict=TRUE, objnames="upsliks"))
            
            if (inherits(chk,"try-error")) { # if found neither in global env nor on disk:
              stop("'upsliks' not stored where expected for re-summarizing.")
            } else cat(crayon::yellow(paste0(" 'upsliks' loaded from disk.\n")))
          }
          max_table_size <- tail(names(upsliks),1L)
          if (cmdl_fast_reClu) { # for reclustering without reprojecting
            message("***** cmdl_fast_reClu => nbCluster set to ", cmdl_fast_reClu," ***** ")
            upsliks[[max_table_size]] <- recluster(upsliks[[max_table_size]],
                                                   nbCluster=cmdl_fast_reClu) 
          }
          if (length(grep("_MAF", modelpath_nickname))) {
            config_mafR(torch_device = "cuda")
            if (final_MAF && ! inherits(upsliks[[max_table_size]]$jointdens,"MAF")) {
              message("final_MAF: recluster() called before computation of summaries.")
              # 'recluster' is here a misnomer: MAFs are used instead of multivariate Gaussian clustering.
              if (as_in_PSM19 <- TRUE) { # Use controls from Papamakarios et al 19 by default.
                message("MAF controls as in PSM19.")
                Infusion.options( 
                  MAF_batchsize = function(...) 100L, 
                  MAF_validasize = function(nr, ...) {nr %/% 20L}, 
                  MAF_patience=20,
                  design_hidden_layers= function(projdata, nc=ncol(projdata),
                                                 transforms, # autoregressive layers 'K', 5 by default ['MAF_auto_layers' in Infusion.options()]
                                                 n_hid_layers=2L, # The L ~ 2 of PPM17 ans PSM19
                                                 n_units_per_layer=50L,
                                                 ...
                  ) { 
                    resu <- rep(n_units_per_layer, n_hid_layers) 
                    ## Info: approximate number of parameters according to Table 1 SM PapamakariokPM17:
                    npMAF <- (3L*transforms*nc*n_units_per_layer + 
                                transforms*(n_hid_layers-1)*n_units_per_layer*n_units_per_layer
                    )/2L
                    attr(resu,"info") <- c(npMGM=NA_integer_, npMAF= npMAF)
                    resu
                  },
                  Adam_learning_rate=1e-4 # PSM19
                )
              } else { # Similar to the zuko tutorials
                message("MAF controls similar to the zuko tutorials.")
                Infusion.options( 
                  MAF_batchsize = function(...) 100L, 
                  MAF_validasize = function(nr, ...) {nr %/% 20L}, 
                  MAF_patience=30, # different
                  MAF_auto_layers=3L,
                  design_hidden_layers= function(projdata, nr= nrow(projdata), nc=ncol(projdata),
                                                                       transforms, # autoregressive layers: the K ~ 5 of PPM17; cf MAF_auto_layers for default value
                                                                       n_hid_layers=2L, # The L ~ 2 of PPM17 ans PSM19
                                                                       design_fac, 
                                                                       ...) { 
                    # too high value => very poor fit => not good candidates when sampling pardens
                    # occurred with nr= ~ 6000, nc=4, hidden_units= rep(3000L,3)
                    n1 <- (Infusion.getOption("maxnbCluster"))(nr=nr,nc=nc)
                    n2 <-  nr^0.31
                    nbClu <- max(1L, as.integer(min(n1,n2))) # effectively the # of clusters of MGM fit
                    npMGM <-  ((nc*(nc+3L))%/%2L+1L)*nbClu-1L # exact # of param of MGM fit
                    n_units_per_layer <- (- 3*nc + sqrt(9*nc*nc+8*npMGM/transforms) )/2 # O(npMGM)
                    n_units_per_layer <- n_units_per_layer*design_fac # allows ad-hoc adjustment
                    n_units_per_layer <- 2^as.integer(3+ # so that ~ npMAF is power of two,  O(8 to 16 times npMGM)
                                                        log(n_units_per_layer,base=2))
                    n_units_per_layer <- as.integer(n_units_per_layer +1/(100*npMGM)) # safe rounding
                    resu <- rep(n_units_per_layer,n_hid_layers) # zuko.flows.MAF()'s hidden_features argument:
                    ## => the number of elements gives the number of hidden layers.
                    
                    ## Info: approximate number of parameters according to Table 1 SM PapamakariokPM17:
                    npMAF <- (3L*transforms*nc*n_units_per_layer + 
                                transforms*(n_hid_layers-1)*n_units_per_layer*n_units_per_layer
                    )/2L
                    attr(resu,"info") <- c(npMGM=npMGM, npMAF= npMAF)
                    
                    resu
                  },
                  Adam_learning_rate=1e-3 
                )
                
              }
              upsliks[[max_table_size]] <- recluster(upsliks[[max_table_size]], using="c.mafR") 
              # Thepython objects must be savedseparately from the R ones:
              upsliks[[max_table_size]] <- 
                save_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", 
                        ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))  
              upsliks  <- upsliks[max_table_size]
              if (TRUE) {
                deforest_projectors(upsliks[[max_table_size]]) # otherwise saving may be *slow* and *bulky* 
                save(upsliks, file=paste0("upsliks_",ii,".rda"))
              }
              MSL_recomputed <- TRUE
            } else upsliks[[max_table_size]] <- 
               load_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", 
                         ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))
             
          } 
          if (cmdl_fast_reClu && ! MSL_recomputed) { # for reclustering without reprojecting
            message("***** cmdl_fast_reClu => nbCluster set to ", cmdl_fast_reClu," ***** ")
            upsliks[[max_table_size]] <- recluster(upsliks[[max_table_size]],
                                                   nbCluster=cmdl_fast_reClu) 
            if (length(grep("_MAF",modelpath_nickname))) {
              upsliks[[max_table_size]] <- save_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))
              
            }
            MSL_recomputed <- TRUE
          } 
          
          if (cmdl_resummarize  && ! MSL_recomputed) {
            message("resummarize: MSL() called before computation of summaries.")
            upsliks[[max_table_size]] <- MSL(upsliks[[max_table_size]], CIs=FALSE, eval_RMSEs = FALSE)
          }
        } else { # Full workflow using MAF. Not reported in ms, but doable. 
          # You will need a lot of available GPU if you want to run this on many simulated data sets.
          if (length(grep("_MAF",modelpath_nickname))) {
            using <- "c.mafR"
            config_mafR(torch_device = "cuda")
            message("************** 'c.mafR' is being used. ***************** ")
            # don't forget the export in the bash script ...
          } else using <- "Rmixmod"
          upsliks <- get_upsliks(
            ii=ii,  # fn def'd in Infusion_addl_defs.R
            S_obs = S_obs,
            workflow_design=workflow_design,
            reftable_sizes=reftable_sizes,
            control.Simulate=control.Simulate,
            LOWER=LOWER, UPPER=UPPER,
            methodArgs=list(importance="none"), 
            init_indiv_reft_size=init_indiv_reft_size,
            more_refines=cmdl_refine,
            cluster_args=cluster_args,
            using=using,
            verbose= list(most=interactive() || on_slurm,
                          rparam=interactive() || on_slurm,
                          cloud_parm=if (Infusion:::.is_devel_session()) {
                            intersect(c("log1p.t23.","d3_Nbn3", names(LOWER)),
                                      names(LOWER))[1]
                          } else {NULL}),
            scaDGparlist=scaDGparlist,
            constr_crits=constr_crits
          )
          max_table_size <- tail(names(upsliks),1L)
          if (length(grep("_MAF",modelpath_nickname))) {
             upsliks[[max_table_size]] <- save_MAFs(upsliks[[max_table_size]],prefix="saved_fits/", ext=paste0(".",max_table_size,".",ii,".MAF.pkl"))
          }
          
        }
        cat(paste0("\nmax_table_size=",max_table_size,"\n"))
        
        ## return value of each replicate is 'summaries'
        if (cmdl_reproject) {
          message(paste("Reprojecting with num.trees=",cmdl_reproject))
          upsliks[[max_table_size]] <- reproject(upsliks[[max_table_size]], 
                                                 methodArgs=list(num.trees=cmdl_reproject))
        }
        
        ## Control of number of bootstrap replicates
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
        summaries <- summary_inference_K( # fn def'd in performance_fns.R 
          slik_siz=upsliks[[max_table_size]], 
          # nsim  for bootstrap LRT: 
          nsim=nsim, # in interactive session small bootstrap still useful to catch bugs
          # Note two default #cores controls.
          nb_cores_LRboot= if (interactive()) {
            if (nsim>19L) {min(procs_per_job, 50L %/% length(LOWER))} else {1L}
          } else {min(procs_per_job, 50L %/% length(LOWER))}, # 50G 7 params 9 cores -> several failures
          # on win, very roughly for N_7from17, time = C+2.6 B/C for B replicates, C cores ;
          # => time minimized for C= sqrt(2.6 * B)
          # => In this and higher dimensional pbs, in practice this means 
          #   using all available processors except for toy (low B) computations
          h0s=h0s_ii, ret_HPDint=TRUE,
          abcrf_proj_K=abcrf_proj_all_K[[max_table_size]], 
          verbose=(interactive() || on_slurm),
          reftable_abcrf=reftable_abcrf_final[seq(na.omit(reftable_sizes[max_table_size])),], # ...protect against NA here after refine
          rstar_R=if (length(LOWER)>7L || final_MAF) {0L} else 1000L, # speculative, not used in ms.
          ii =ii
        )
        if (length(grep("_MAF",modelpath_nickname))) {
          summaries$train_times <- c(
            conddens=attr(upsliks[[max_table_size]]$conddens,"train_time"),
            jointdens=attr(upsliks[[max_table_size]]$jointdens,"train_time"),
            pardens=attr(upsliks[[max_table_size]]$pardens,"train_time"),
            postdens=attr(upsliks[[max_table_size]]$postdens,"train_time")
          )
        }
        summaries
      },
      simplify=FALSE
    )
    names(summaries_list) <- seq(prev_latest_one+1L,upto)
    # print(do.call(rbind,do.call(rbind,summaries_list[[1]]$SLRTs)[,"basicLRT"]))
    save(summaries_list, file=paste0(c("summaries_list",cmdl_jobNbr,"rda"), collapse="."))
  }
}

### # merging results from each subdir on genotoul 
if (FALSE) { 
  ## This block of code was run manually after all performance simulations had been run.
  { # Helper functions
    merge_genotoul_jobs <- function(upto, from=1L, interrupted=FALSE) {
      seqft <- seq(from=from,to=upto)
      merger <- vector("list",length(seqft))
      it <- 0L
      for (jobit in seqft) {
        it <- it+1L
        jpath <- paste0("job_",jobit,"/")
        summaries_list <- NULL # SAFETY
        if (interrupted) {
          try(chk <- load(paste0(jpath, "interrupted_summ_",jobit,".rda"))) # loads 'summaries_list'
          if (inherits(chk,"try-error")) {
            merger[it] <- list(NULL)
          } else merger[[it]] <- resu 
        } else {
          try(chk <- load(paste0(jpath, "summaries_list.",jobit,".rda"))) # loads 'summaries_list'
          if (inherits(chk,"try-error")) {
            merger[it] <- list(summaries_list[[1]])
          } else merger[[it]] <- summaries_list[[1]] 
        }
      }
      names(merger) <- seqft
      invisible(merger)
    }

    save_merge_genotoul_jobs <- function(upto, from=1L, interrupted=FALSE) {
      summaries_list <- merge_genotoul_jobs(upto,from, interrupted=interrupted)
      lnames <- names(summaries_list)
      savename <- paste0("summaries_list.",lnames[1],"_",tail(lnames,1), 
                         ".v",packageVersion("Infusion"),".rda")
      save(summaries_list, file=savename)
      invisible(summaries_list)
    }
    
    ## Typically used as 
    # summaries_list <- save_merge_genotoul_jobs( < upto >, < prev_latest_one > +1L)
    
    ##  SEE the Rmarkdown file used to generate tables and Figures for the ms
    ##  for usage of the 'summaries_list' thus produced
  }

}

### Some more ad hoc analyses

if (FALSE) {
  source(paste0(dir_generic,
                "../../generic_workflow/resummarizes_using_synthetic_reftable.R"))
  # On genotoul I used a synthetic_workflow.R, which is saved in the diy2abc/... directory.
  # It was run by running syntParent.sh from the base N_7from17 directory.
  # next one should run (diy2abc::) save_merge_jobs(., prefix="summaries_list.")
  # in the synthetic/ subdirectory.
}


if (FALSE) {
  # Inference from small reftables
  # saved_fits dir and upsliks element may need to be adjusted in the code sourced below
  source(paste0(dir_generic,"../../generic_workflow/resummarize_saved_fits.R"))
}