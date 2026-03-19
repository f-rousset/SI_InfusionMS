safe_load <- function(fileroot,
                      strict=TRUE, # check the name(s) of the object in the file
                      objnames=fileroot
                      ) {  
  chk <- NULL
  for (objname in objnames) {
    if ( ! exists(objname,where = .GlobalEnv)) {
      Rdname <- paste0(fileroot,".rda")
      chk <- suppressWarnings(try(load(Rdname, envir=.GlobalEnv), silent=TRUE))
      if (inherits(chk,"try-error")) {
        if (length(intersect(dir(), Rdname))) { # dir(.,fixed = TRUE) only in recent versions of R
          message(paste0(Rdname," found, but load() failed"))
        } else message(paste0(Rdname," not found"))
        break
      } else if (strict && ! objname %in% chk) stop(paste0(Rdname," file has been read but does not contain object ",objname))
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
      warning("mydetectCores() called on genotoul.")
  } else {
    # assuming linux-like OS where this may work:
    # Dirk Eddelbuettel, "Re: [Rd] Detecting physical CPUs in detectCores() on Linux platforms" on R-devel list, 2023/08/08 
    # 4:16 message, and 2:07 message for more heuristic fallback.
    
    # But this apparently returned NA when run through slurm on genotoul
    bla <- try(system("hwloc-info", intern=TRUE), silent=TRUE)
    # => This system call presumes that the libhwloc-dev and hwloc (debian) packages have been installed. 
    # Doc: https://www.open-mpi.org/projects/hwloc/doc/
    if (inherits(bla,"try-error")) {
      message("mydetectCores() tried to run 'hwloc-info' but this failed.")
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

