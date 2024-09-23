safe_load <- function(fileroot,
                      strict=TRUE, # check the name(s) of the object in the file
                      objnames=fileroot
                      ) {  
  chk <- NULL
  for (objname in objnames) {
    if ( ! exists(objname,where = .GlobalEnv)) {
      chk <- suppressWarnings(try(load(paste0(fileroot,".rda"), envir=.GlobalEnv), silent=TRUE))
      if (inherits(chk,"try-error")) {
        message(paste0(fileroot,".rda is not available."))
        break
      } else if (strict && ! objname %in% chk) stop(paste0(objname,".rda file found but does not contain ",objname))
    }
  }
  chk
}

# specifically exclude hyperthreading:
mydetectCores <- function() {
  if (.Platform$OS.type=="windows") {
    max_nb_cores <- parallel::detectCores(logical = FALSE)
  } else if (length(grep("node|genologin",Sys.info()["nodename"]))) {
      max_nb_cores <- 1L # avoid nested parallelisation when running multiple analyses on the cluster
      # ___F I X M E___ this is unclean, as in that I may allow nested parallelization on genotoul 
      warning("mydetectCores() called on genotoul.")
  } else {
    # assuming linux-like OS where this may work:
    # Dirk Eddelbuettel, "Re: [Rd] Detecting physical CPUs in detectCores() on Linux platforms" on R-devel list, 2023/08/08 
    # 4:16 message, and 2:07 message for more heuristic fallback.
    
    # But this apparently returned NA when run thorugh slurm on genotoul
    
    bla <- try(system("hwloc-info", intern=TRUE), silent=TRUE)
    # => This system call presumes that the libhwloc-dev and hwloc (debian) packages have been installed. 
    # Doc: https://www.open-mpi.org/projects/hwloc/doc/
    if (inherits(bla,"try-error")) {
      # then using the more heuristic approach:
      d <- read.dcf("/proc/cpuinfo") 
      max_nb_cores <- as.integer(unique(d[, grep("cpu cores",colnames(d))]))
    } else {
      max_nb_cores <- tail(
        strsplit(
          strsplit(bla[grep("Core", bla)]," Core")[[1]][1],
          " ")[[1]],
        1)  
      max_nb_cores <- as.integer(max_nb_cores)
    }
  }
  max_nb_cores
}

