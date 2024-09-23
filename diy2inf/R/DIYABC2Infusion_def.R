header2pars <- function(headerFile, actual_mut_pars=NULL) {
  settings <- readLines(headerFile)
  upto <- grep("priors", settings)[1L]-1L
  bla <- settings[1:upto]
  from <- grep("scenario 1", bla)+1L
  bla <- bla[from:upto]
  bla <- paste(bla, collapse=" ")
  bla <- gsub(" sample "," ", bla)
  bla <- gsub(" merge "," ", bla)
  bla <- gsub(" varNe "," ", bla)
  bla <- gsub(" split "," ", bla)
  bla <- gsub("\\+"," ", bla)
  bla <- gsub("\\-"," ", bla)
  while(length(grep(" [0-9]* ", bla))) bla <- gsub(" [0-9]* "," ", bla)
  bla <- strsplit(bla, "\\s+")[[1]]
  bla <- unique(bla)
  
  
  
  ## check of consistency of order of ranges
  line_hist_priors <- grep("historical parameters priors", settings)
  hist_par_nbr <- strsplit(settings[line_hist_priors],"(",fixed=TRUE)[[1]][2]
  hist_par_nbr <- strsplit(hist_par_nbr,",",fixed=TRUE)[[1]][1]
  hist_par_nbr <- as.integer(hist_par_nbr)
  rangenames <- character(hist_par_nbr)
  for (par_it in seq_len(hist_par_nbr)) {
    ligne <- settings[line_hist_priors+par_it]
    parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
    rangenames[par_it] <- parseline[1]
  }
  
  if (hist_par_nbr #+mut_par_nbr    That's the number in 'historical parameters priors (XXX,0)'
      != length(bla)) {
    stop("Number of parameters deduced from scenario description differs\n from that deduced from specification of priors.")
  }
  
  if ( ! all(bla==rangenames)) {
    print(bla)
    print(rangenames)
    stop("order of parameters deduced from scenario description differs\n from that deduced from specification of priors.")
    # HOWEVER the workflow is presumably insensitive to the order of params on the last 
    # line of the *template* bc 
    # * For diyabc -o the parameters in the last line are not copied from the template 
    #   into the header read by diyabc but are those found to be variable in the table of parameters.
    # * Vanilla diyabc is called only for fully fixed parameter values so no parameter is variable
    #   so there is no parameter on the last line of the header read by diyabc 
  } 
  
  
  
  curr_line <- line_group_priors <- grep("group priors", settings) # ___F I X M E___ returns 1 rather than the line...
  mut_pars <- c()
  if (length(curr_line)) {
    grp_nbr <- strsplit(settings[line_group_priors],"(",fixed=TRUE)[[1]][2]
    grp_nbr <- strsplit(grp_nbr,")",fixed=TRUE)[[1]][1]
    grp_nbr <- as.integer(grp_nbr)
    for (grp_it in seq_len(grp_nbr)) {
      curr_line <- curr_line+1L
      parse_grp_line <- strsplit(settings[curr_line]," ",fixed=TRUE)[[1]]
      if (parse_grp_line[3]=="[M]") {
        for (line_it in curr_line+seq_len(6)) {
          ligne <- settings[line_it]
          parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
          mut_pars <- c(mut_pars,parseline[1])
        }
      } else stop("Unhandled marker type")
    }
  }
  mut_pars <- intersect(mut_pars, actual_mut_pars)
  mut_par_nbr <- length(mut_pars)
  rangenames <- c(rangenames, mut_pars)
  
  # line_scenarios <- grep(" scenarios: ", settings)
  # scenar_nbr <- strsplit(settings[line_scenarios]," ",fixed=TRUE)[[1]][1]
  # scenar_nbr <- as.integer(scenar_nbr)
  # => Not sure how multiple scenarios should be handled here, but there is no such case in my simuls.
  last_line <- settings[length(settings)]
  parseline <- strsplit(last_line," ",fixed=TRUE)[[1]]
  parnamesFromLastLine <- parseline[2L:(1L + hist_par_nbr + mut_par_nbr)]
  if ( ! all(rangenames==parnamesFromLastLine)) {
    print(rangenames)
    print(parnamesFromLastLine)
    stop("order of parameters deduced from scenario description differs from that on last line.")
  } 
  
  rangenames
}

# Two uses: called by DIY2Inf_readRefTable(), using the 'header' argument; 
#            also user-level in master script, using the 'format' and 'headerFile' arguments
get_par_ranges <- function(headerFile=NULL, # file name
                           header, # character vector, read from header file as shown below
                           nparamtoth=NULL, format="list") {
  if ( ! is.null(headerFile)) {
    header <- readr::read_file(headerFile)
    header <- gsub("\\r\\n", "\n", header)
  }
  if (is.null(nparamtoth)) {
    reparamtot <- "historical parameters priors \\((\\d+)\\D"
    nparamtoth <- as.integer(stringr::str_match(stringr::str_extract(header, reparamtot)[[1]], 
                                               reparamtot)[2])
  }
  reparamlist <- paste0("\\bhistorical parameters priors.*\\n((?:\\w+\\W[^\\n]*\\n){", 
                       nparamtoth, "})")
  paramlistmatch <- stringr::str_match(stringr::str_extract_all(header, reparamlist), 
                                      reparamlist)[2]
  reparam <- "(\\w+)\\W+\\w\\W+\\w\\w\\[([^,\\]]+),([^,\\]]+)[,\\]][^\\n]*\\n"
  par_ranges <- stringr::str_match_all(stringr::str_extract_all(paramlistmatch, reparam)[[1]], 
                                   reparam)
  if (format=="matrix") {
    par_ranges <- do.call(rbind, par_ranges)
    par_ranges <- matrix(as.numeric(par_ranges[,3:4]), ncol=2, dimnames=list(par_ranges[,2], NULL))
  }
  par_ranges
}

# I deduced this one by trial and error from the previous one...
get_mut_par_ranges <- function(headerFile, actual_mut_pars, format="matrix") {
  lignes <- readLines(headerFile) 
  par_ranges <- list()
  for (mut_parname in actual_mut_pars) {
    reparam <- "(\\w+)\\W+\\w\\w\\[([^,\\]]+),([^,\\]]+)[,\\]][^\\n]*"
    par_ranges <- c(par_ranges,
                    stringr::str_match_all(lignes[grep(mut_parname,lignes)[1]],  reparam) )
  }
  if (format=="matrix") {
    par_ranges <- do.call(rbind, par_ranges)
    par_ranges <- matrix(as.numeric(par_ranges[,3:4]), ncol=2, dimnames=list(par_ranges[,2], NULL))
  }
  par_ranges
}


# called by wrap_read_reftable.bin() hence still used
# From code provided by RL 2023/01/25
DIY2Inf_readRefTable <- function(filename = "reftableRF.bin", headername="headerRF.txt", N=1) {
  header = readr::read_file(headername)
  header <- gsub("\\r\\n", "\n", header)
  rescen = "\\bscenario\\s+(\\d+)\\s+.*\\n((?:(?!(?:scenario|\\n)).*\\n)+)"
  scenprematch = stringr::str_extract_all(header, rescen)[[1]]
  scenprematch
  nscenh = length(scenprematch)
  scenmatch = stringr::str_match_all(scenprematch, rescen)
  scendesc <- vector(mode = "character", length = nscenh)
  for (i in 1:nscenh) {
    scendesc[as.integer(scenmatch[[i]][2])] = scenmatch[[i]][3]
  }
  reparamtot = "historical parameters priors \\((\\d+)\\D"
  nparamtoth = as.integer(stringr::str_match(stringr::str_extract(header, reparamtot)[[1]], 
                                             reparamtot)[2])
  paramsh <- get_par_ranges(header=header, nparamtoth=nparamtoth)
  #
  paramsdesc = list()
  reali <- 1L # count variable parameters (but locally shifted by one for paramsdesc[paramsh[[i]][2]] <- reali assignment)
  for (i in 1:nparamtoth) { 
    mini = as.numeric(paramsh[[i]][3])
    maxi = as.numeric(paramsh[[i]][4])
    if (maxi != 0) {
      crit <- ((maxi - mini)/abs(maxi) > 1e-06)
    } else  crit <- (maxi - mini > 1e-06)
      if (crit) {
        paramsdesc[paramsh[[i]][2]] <- reali
        reali = reali + 1L
      }
  }
  realparamtot <- reali - 1L
  #
  parambyscenh <- vector(mode = "numeric", length = nscenh)
  for (i in 1:nscenh) {
    templist = list()
    listterms = strsplit(scendesc[i], "\\W")[[1]]
    m = 1
    for (j in 1:length(listterms)) {
      if (!is.null(paramsdesc[listterms[j]][[1]])) {
        templist[m] = paramsdesc[listterms[j]][[1]]
        m = m + 1
      }
    }
    parambyscenh[i] = list((unique(unlist(templist))))
  }
  restatsname = "\\n\\nscenario\\s+.*"
  allcolspre = tail(strsplit(stringr::str_extract(header, restatsname), 
                             "\\s+")[[1]], -2)
  to.read <- file(filename, "rb")
  realnrec = readBin(to.read, integer(), endian = "little")
  nrec = if (N > 0) {
    min(N, realnrec)
  }
  else {
    realnrec
  }
  nscen = readBin(to.read, integer(), endian = "little")
  nrecscen = readBin(to.read, integer(), n = nscen, endian = "little")
  nparam = readBin(to.read, integer(), n = nscen, endian = "little")
  nstat = readBin(to.read, integer(), endian = "little")
  paramsname = head(allcolspre, -nstat)
  # OK for short line header
  nmutparams = length(paramsname) - realparamtot
  #nmutparams = length(paramsname) - nparamtoth # was realparamtot
  stats = matrix(nrow = nrec, ncol = nstat)
  colnames(stats) <- tail(allcolspre, nstat)
  params = matrix(nrow = nrec, ncol = realparamtot + nmutparams)
  # OK for short line header
  colnames(params) <- paramsname
  # now (but won't work for all cases of fixed/var params, Be Carefull)
  #colnames(params) <- paramsname[(nparamtoth - realparamtot + 1): (nparamtoth+nmutparams)]
  scenarios = vector(mode = "numeric", length = nrec)
  for (i in 1:nrec) {
    scen = readBin(to.read, integer(), endian = "little")
    scenarios[i] = scen
    lparams = readBin(to.read, numeric(), n = nparam[scen], 
                      size = 4, endian = "little")
    for (j in 1:length(parambyscenh[[scen]])) {
      params[i, parambyscenh[[scen]][j]] = lparams[j]
    }
    if (nmutparams > 0) {
      for (jm in 1:nmutparams) {
        params[i, realparamtot + jm] = lparams[nparam[scen] - 
                                                 nmutparams + jm]
      }
    }
    lstats = readBin(to.read, numeric(), n = nstat, size = 4, 
                     endian = "little")
    for (j in 1:nstat) {
      stats[i, j] = lstats[j]
    }
  }
  close(to.read,"rb")
  
  list(nrec = nrec, nscen = nscen, nrecscen = nrecscen, nparam = nparam, 
       scenarios = as.factor(scenarios), params = params, stats = stats)
}

generate_dir_name <- function(base="tmp", nonbase="[0-9]+") {
  ## pattern <- paste(base,"*",sep="") # bah non, that is regexp which means 
  # tm matches the characters tm literally (case sensitive)
  # p* matches the character p ZERO or more times.
  # so this matches html files....
  # Rather, pattern <- paste(base,".*",sep="") would match tmp files
  pattern <- paste(base,nonbase,sep="")
  allmatches <- dir(pattern=pattern)
  allremainders <- substring(allmatches,nchar(base)+1L)
  # If a remainder is not numeric, NA's are produced by the next line: (FIXME)
  allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  if (length(allremainders) == 0L) {
    num <- 0L
  } else num <- max(allremainders)+1L
  validname <- paste ( base , num,sep="") 
  return(validname)
}

# Takes a standard DIY header file 
# and returns an ad hoc reader with par_vector values interpreted as trivial prior distributions:
# + I added check of consistency of the header file with the par vector
parvec2DIYheader <- function(inputHeader, par_vector, outputHeader, verbose.) {
  settings <- readLines(inputHeader)
  line_scenarios <- grep(" scenarios: ", settings)
  scenar_nbr <- strsplit(settings[line_scenarios]," ",fixed=TRUE)[[1]][1]
  scenar_nbr <- as.integer(scenar_nbr)
  # 
  # Fo consistency check
  line_scenario1 <- grep("scenario 1 ", settings)
  par_nbr_check <- strsplit(settings[line_scenario1],"(",fixed=TRUE)[[1]][2]
  par_nbr_check <- strsplit(par_nbr_check,")",fixed=TRUE)[[1]][1]
  par_nbr_check <- as.integer(par_nbr_check)
  #
  line_hist_priors <- grep("historical parameters priors", settings)
  par_nbr <- strsplit(settings[line_hist_priors],"(",fixed=TRUE)[[1]][2]
  par_nbr <- strsplit(par_nbr,",",fixed=TRUE)[[1]][1]
  par_nbr <- as.integer(par_nbr)
  if (par_nbr_check != par_nbr) {
    warning("Number of parameters specified in 'scenario 1' and 'historical parameters priors' lines should be identical.",
            immediate.=TRUE)
  }
  for (line_it in line_hist_priors+seq_len(par_nbr)) {
    ligne <- settings[line_it]
    parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
    parname <- parseline[1]
    if (! is.na(parval <- par_vector[parname])) parseline[3] <- paste0("UN[",parval,",",parval,",0.0,0.0]")
    settings[line_it] <- paste(parseline, collapse = " ")
  }
  
  curr_line <- line_group_priors <- grep("group priors", settings)
  mut_par_nbr <- 0L
  if (length(curr_line)) { # TRUE for Harmonia microsat example -> muarker parameters
    grp_nbr <- strsplit(settings[line_group_priors],"(",fixed=TRUE)[[1]][2] # markers= marker groups =marker types (?) 
    grp_nbr <- strsplit(grp_nbr,")",fixed=TRUE)[[1]][1]
    grp_nbr <- as.integer(grp_nbr)
    for (grp_it in seq_len(grp_nbr)) {
      curr_line <- curr_line+1L
      parse_grp_line <- strsplit(settings[curr_line]," ",fixed=TRUE)[[1]]
      if (parse_grp_line[3]=="[M]") {
        for (line_it in curr_line+seq_len(6)) {
          ligne <- settings[line_it]
          parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
          parname <- parseline[1]
          parseparse2 <- strsplit(parseline[2],"[",fixed=TRUE)[[1]]
          parseparse22 <- strsplit(parseparse2[2],",",fixed=TRUE)[[1]]
          parvalue <- par_vector[parname]
          if (is.na(parvalue) && parname=="MEANMU") parvalue <- 10^par_vector["logMEANMU"]
          if (is.na(parvalue)) {
            if (verbose.) message(paste("Parameter",parname, "not specified, left as-is."))
          } else {
            parseparse22[1:2] <- parvalue
            parseparse22 <- paste(parseparse22,collapse=",")
            parseparse2[2] <- parseparse22
            parseline[2] <- paste(parseparse2,collapse="[")
            settings[line_it] <- paste(parseline, collapse = " ")
            mut_par_nbr <-   mut_par_nbr + 1L
          }
        }
      } else stop("Unhandled marker type")
      if (par_nbr+mut_par_nbr != length(par_vector)) {
        warning((paste0("parameter number deduced from header file != length of parameter vector (",
                        par_nbr+mut_par_nbr," vs ",length(par_vector),")")
        ), immediate.=TRUE)
      }
      #
    }
  } # but FALSE for SNP bottleneck example
  
  # When a single sample is generated for given parameter values, 
  #   the generated header's last line should contain the statnames but not the parnames
  #   An error here suggests a wrong accounting of the different types of parameters above.
  last_line <- settings[length(settings)]
  parseline <- strsplit(last_line," ",fixed=TRUE)[[1]]
  parseline <- parseline[-c( (scenar_nbr + 1) : (scenar_nbr + par_nbr + mut_par_nbr + 1))]
  settings[length(settings)] <- paste(parseline, collapse = " ")
  writeLines(text = settings,con = outputHeader)
}

# called by DIYABC_wrapper1() hence still used
wrap_read_reftable.bin <- function(outputHeader, exe_path=getwd(), tidy) {
  # the distinct path seems useful here.
  tmpdirname <- generate_dir_name()
  tmpdir <- paste0(exe_path,"/",tmpdirname,"/")
  dir.create(tmpdir)
  file.copy("reftableRF.bin",tmpdir)
  file.copy(outputHeader,tmpdir) 
  setwd(tmpdir)
  reftable <- DIY2Inf_readRefTable(filename = "reftableRF.bin", headername=outputHeader, N=1)
  setwd(exe_path)
  if (tidy) unlink(tmpdir, recursive = TRUE)
  reftable
}

# called by onesim() or onesim_c() hence still used
DIYABC_wrapper1 <- function(par_vector,
                            joint=FALSE, # Whether to join the par_vector to the summary stats in the result
                            control.Simulate,
                            exe_path=control.Simulate$exe_path,
                            inputHeader=control.Simulate$inputHeader,
                            outputHeader=control.Simulate$outputHeader,
                            # NOT the control.Simulate$simulator_call:
                            simulator_call, # ='./diyabc-RF-macos-v1.1.27 -p ./ -R "ALL" -r 1 -g 1 -m -t 1',
                            # as the proper one is without -o here
                            verbose.=FALSE, # Whether to suppress diyabc stdout and some other info
                            tidy=TRUE # Whether to remove temporary paths and files
                            ) {
  par_vector <- unlist(par_vector)
  par_vector <- canonize(par_vector)
  if (is.null(exe_path)) {
    stop("An 'exe_path' to the DIYABC executable must be provided.")
  } else {
    oripath <- setwd(exe_path)
  }
  
  # What for ?: 
  if(length(files <- list.files(pattern = "reftable.*.bin", ignore.case = TRUE)) > 0 ) file.remove(files)
  
  # Create ad hoc reader with par_vector values interpreted as trivial prior distributions:
  if ( ! is.null(reparamHeader <- control.Simulate$reparamHeaderName)) {
    reparam_vector <- reparametrize(par_vector, to=control.Simulate$reparamNames, fixedPars=control.Simulate$fixedPars)
    parvec2DIYheader(inputHeader=reparamHeader, par_vector=reparam_vector, outputHeader=outputHeader, verbose.=verbose.)
  } else parvec2DIYheader(inputHeader=inputHeader, par_vector=par_vector, outputHeader=outputHeader, verbose.=verbose.)
  #### call to diyabc simulator:
  simulator_call <- reformat_simulator_call(simulator_call, nsim=1L)
  chk <- system(simulator_call, ignore.stdout = ! verbose.) # IF NEVER TERMINATES: check that constraints are satisfied on the parameters
  if (chk) {
    warning(paste("Problem with simulator call, system return code=",chk))
    #if (chktime["elapsed"]<3) {
      message("The stdout of the system call is:")
      bla <- system(simulator_call, ignore.stdout = FALSE)
      str(bla)
      #info <- system(simulator_call, intern=TRUE)
      
    #}
  }
  
  #### conversion to reftable in the DIYabc style:
  reftable <- try(wrap_read_reftable.bin(outputHeader=outputHeader, exe_path=exe_path, tidy=tidy))
  if (inherits(reftable,"try-error")) {
    mess <- paste("Run",simulator_call,"to reproduce problem.\n",
                  "If data file appears missing, this refers to file specified on first line of header file.")
    stop(mess)
  }
  

  if (tidy) unlink(dir(pattern = "first_records_of_the_reference_table*"))
  # cat("reftable structure:\n ")
  # str(reftable)
  # if(reftable$nrec < 50 && length(reftable$stats) < 100 ) {
  #   cat("complete reftable:\n ")
  #   print(reftable)
  # }
  setwd(oripath)
  #### conversion to row of reftabel in the Infusion style:
  if (joint) {
    as.data.frame.list(c(par_vector, reftable$stats[1,]))
  } else reftable$stats[1,]
}

if (FALSE) {
  ## Check that this works without of the onesim wrapper defined below:
  DIYABC_wrapper1(exe_path="C:/home/francois/travail/stats/Infusionplus/caseStudies/Harmonia/HA_for_testing_infusion_RF",
                  simulator_call, # ='./diyabc-RF-windows-v1.1.27.exe -p ./ -R "ALL" -r 1 -g 1',
                  par_vector = c(N1=1000, N2=1000, N3=1000, N4=1000, T1=13, ar=0.5, logMEANMU=-4, MEANP=0.25))
}

.DEPARSE <- get(".DEPARSE", asNamespace("spaMM"))

# with canonical parameters (is that used ?)
onesim <- function(N1,N2,N3,N4, ar,logMEANMU, MEANP, ...) {
  # cat(".")
  mc <- match.call()
  mc[[1]] <- NULL
  DIYABC_wrapper1(par_vector = unlist(eval(parse(text=.DEPARSE(mc[1:8])))), ...)
}

get_onesim_c <- function(scaNames) {
  onesim_c <- function(...) {
    mc <- match.call(expand.dots = FALSE)
    mc[[1]] <- NULL
    parNames <- setdiff(names(mc),"...")
    DIYABC_wrapper1(par_vector = unlist(eval(parse(text=.DEPARSE(mc[parNames])))), ...)
  }
  expr <- paste0("alist(",
                 paste0(scaNames, collapse="=, "), 
                 "=, ...=)")
  formals(onesim_c) <-   eval(parse(text=expr))
  onesim_c
}

# Used in up-to-date workflow
ranges2DIYheader <- function(inputHeader, 
                             ranges, # only for variable parameters
                             outputHeader, verbose., update_ranges=FALSE) {
  settings <- readLines(inputHeader)
  line_scenarios <- grep(" scenarios: ", settings)
  scenar_nbr <- strsplit(settings[line_scenarios]," ",fixed=TRUE)[[1]][1]
  scenar_nbr <- as.integer(scenar_nbr)
  #
  line_hist_priors <- grep("historical parameters priors", settings)
  par_nbr <- strsplit(settings[line_hist_priors],"(",fixed=TRUE)[[1]][2]
  par_nbr <- strsplit(par_nbr,",",fixed=TRUE)[[1]][1]
  par_nbr <- as.integer(par_nbr)
  for (line_it in line_hist_priors+seq_len(par_nbr)) {
    if (update_ranges) {
      ligne <- settings[line_it]
      parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
      parname <- parseline[1]
      if (parname %in% colnames(ranges)) { # if variable param
        # in -o case 'ranges' is determined by the parsTable ; 
        # here slightly broadened ranges are written to the output header
        # to make sure that the parameters are numerically ]within[ these ranges.
        # But in the current simulation design the parameters are drawn 
        # ]within[ the template-header ranges so this safety is not required.
        parrange <- ranges[,parname]
        parrange <- parrange + c(-1,1)*1e-2
        parseline[3] <- paste0("UN[",parrange[1L],",",parrange[2L],",0.0,0.0]")
        settings[line_it] <- paste(parseline, collapse = " ")
      }
    } # otherwise only trivial iteration to advance line_it
  }
  
  mut_par_nbr <- 0L
  curr_line <- line_group_priors <- grep("group priors", settings)
  if (length(curr_line)) {
    grp_nbr <- strsplit(settings[line_group_priors],"(",fixed=TRUE)[[1]][2]
    grp_nbr <- strsplit(grp_nbr,")",fixed=TRUE)[[1]][1]
    grp_nbr <- as.integer(grp_nbr)
    parnames <- colnames(ranges)
    for (grp_it in seq_len(grp_nbr)) {
      curr_line <- curr_line+1L
      parse_grp_line <- strsplit(settings[curr_line]," ",fixed=TRUE)[[1]]
      if (parse_grp_line[3]=="[M]") {
        for (line_it in curr_line+seq_len(6)) {
          ligne <- settings[line_it]
          parseline <- strsplit(ligne," ",fixed=TRUE)[[1]]
          parname <- parseline[1]
          if (parname %in% parnames) {
            parseparse2 <- strsplit(parseline[2],"[",fixed=TRUE)[[1]]
            parseparse22 <- strsplit(parseparse2[2],",",fixed=TRUE)[[1]]
            parrange <- ranges[,parname]
            #if (is.na(parvalue) && parname=="MEANMU") parvalue <- 10^par_vector["logMEANMU"]
            parrange <- parrange + c(-1,1)*1e-2
            parseparse22[1:2] <- parrange
            parseparse22 <- paste(parseparse22,collapse=",")
            parseparse2[2] <- parseparse22
            parseline[2] <- paste(parseparse2,collapse="[")
            settings[line_it] <- paste(parseline, collapse = " ")
            mut_par_nbr <-   mut_par_nbr + 1L
          } else {
            if (verbose.) message(paste("Parameter",parname, "not specified, left as-is."))
          }
        }
      } else stop("Unhandled marker type")
      
      #
    }
  }
  last_line <- settings[length(settings)]
  parseline <- strsplit(last_line," ",fixed=TRUE)[[1]]
  # parseline <- parseline[-c( (scenar_nbr + 1) : (scenar_nbr + par_nbr + mut_par_nbr + 1))] 
  #                      => this  forgot to put back (variable) parameters. Rather:
  parseline <- c(parseline[1:scenar_nbr], 
                 colnames(ranges), 
                 parseline[(scenar_nbr + par_nbr + mut_par_nbr + 2L):length(parseline)])
  settings[length(settings)] <- paste(parseline, collapse = " ")
  writeLines(text = settings,con = outputHeader)
}

# Trying to avoid problems with the simulator call:
# -R "ALL" is put at a place where it does not block the executable... 
# -o it set last in the line: input parameter tabel will be passed right after, by the calling code
# When -t is present, (1) a default -g 5 is set when it is missing (2) -m is set.
# But
# -p ./ is always applied (might need to be changed later) 
# -R "ALL" is always applied (might need to be changed later) 
# other options such as -i for output file are ignored (idem)
# For some form of doc: https://github.com/diyabc/diyabc/blob/master/README.md
reformat_simulator_call <- function(opts, nsim) {
  # return(opts)
  if (.Platform$OS.type=="windows") {
    exe_pattern <- "diyabc.*"
  } else exe_pattern <- "./diyabc.*" 
  exe <- na.omit(stringr::str_extract(strsplit(opts, split=" +")[[1]], pattern=exe_pattern))
  
  if (length(grep(" -o", opts))) {
    if (length(grep("-t [0-9]*", opts, value=TRUE))) {
      t_grmatch <- stringr::str_extract(opts, pattern="-t [0-9]*")
      if (length(grep("-g [0-9]*", opts, value=TRUE))) { # obey explicit -g value
        g_grmatch <- stringr::str_extract(opts, pattern="-g [0-9]*")
        tg_opts <- paste(t_grmatch, "-m", g_grmatch)
      } else {
        nThreads <- as.integer(substring(t_grmatch, 3))
        gsize <- nsim %/% nThreads
        if (nsim %% nThreads!=0L) gsize <- gsize+1L
        tg_opts <- paste(t_grmatch, "-m -g", gsize)
      }
    } else tg_opts <- "-t 1 -g 1"
    resu <- paste('-p ./ -R \"ALL\"', tg_opts, "-o")
  } else {
    tg_opts <- "-r 1 -g 1"
    resu <- paste('-p ./ -R \"ALL\"', tg_opts)
  }
  resu <- paste(exe, resu)
  resu
}

# Used in up-to-date workflow
DIYABC_reftable_wrapper <- function(
    # the list of arguments of the call is constructed from control.Simulate, (to which element may be added)
  #so the arguments of this function must include all elements of control.Simulate, say:
  # names(control.Simulate)
  # [1] "parsTable"         "exe_path"          "inputHeader"       "outputHeader"      "simulator_call"    "verbose."         
  # [7] "tidy"              "reparamHeaderName" "reparamNames"      "fixedPars"         "caParNames"        "statNames"
  ## It does not include control.Simulate itself
    parsTable,
    out_parsTableFile="essaisimpars.txt", # name of the file where scenario + par.grid will be written for reading by diyabc
    #    the 'parsTable' argument name is used by Infusion to call a Simulate function with a parameter table as argument, hence required with diyabc -o
    exe_path=Infusion.getOption("DIYABC_path"),
    inputHeader, # must bedistinct from outputHeader
    outputHeader="headerRF.txt",
    simulator_call, # ='./diyabc-RF-macos-v1.1.27 -p ./ -R "ALL" -o',
    verbose.=FALSE, # Whether to suppress diyabc stdout and some other info
    tidy=TRUE, # Whether to remove temporary paths and files
    statNames,
    caParNames, # if different from the colNames of the 'canonized' parsp
    fix_canon.grid=NULL,
    fix_parsTable=NULL,
    reparamHeaderName=NULL,
    reparamNames=NULL,
    fixedPars=NULL
) {
  if (is.null(exe_path)) {
    stop("An 'exe_path' to the DIYABC executable must be provided.")
  } else {
    oripath <- setwd(exe_path)
  }
  if ( ! is.null(fix_parsTable)) parsTable <- fix_parsTable(parsTable)
  # Handling 1D case incl. fact tha simplify=FALSE would drop the param name too...
  canon.grid <- do.call(rbind, apply(parsTable, 1L, canonize, simplify=FALSE)) 
  if ( ! is.null(reparamNames)) {
    canon.grid <- do.call(
      rbind, apply(canon.grid, 1L, reparametrize, to=reparamNames, fixedPars=fixedPars, simplify=FALSE)) 
  } else {
    if ( ! is.null(fix_canon.grid)) canon.grid <- fix_canon.grid(canon.grid)
    if (missing(caParNames)) {
      message(paste("'caParNames' missing -> all colnames of canonical grid used in DIYABC_reftable_wrapper():\n",
                    paste(colnames(canon.grid), collapse=", ")))
    } else canon.grid <- canon.grid[,caParNames, drop=FALSE]
  }
  ranges <- apply(canon.grid, 2L,range)
  canon.grid <- cbind(1,canon.grid) # for 'scenario 1'
  write.table(canon.grid, file=paste0(exe_path,out_parsTableFile), col.names = FALSE, row.names = FALSE)
  
  # What for ?: ... this file was created by the previouscall to diyabc executable. 
  if(length(files <- list.files(pattern = "reftable.*.bin", ignore.case = TRUE)) > 0 ) file.remove(files)
  
  # Create ad hoc header with par_vector values interpreted as trivial prior distributions:
  # ____F I X M E____ currently result is correct but super messy :
  # ranges2DIYheader is called with default argument update_ranges=FALSE so that input 'ranges' values 
  # are ignored (for historical demographic params only !!). The line is thus only copied from inputHeader.
  if ( ! is.null(reparamHeaderName)) { # condition modified 2024/08/07 to handle fixedPars without reparamHeader
    ranges2DIYheader(inputHeader=reparamHeaderName, ranges=ranges, outputHeader=outputHeader, verbose.=verbose.) 
  } else ranges2DIYheader(inputHeader=inputHeader, ranges=ranges, outputHeader=outputHeader, verbose.=verbose.)
  # => the final line of the outputHeader contains the names of the variable parameters only = those that remain in the canon.grid
  
  #### call to diyabc simulator:
  simulator_call <- reformat_simulator_call(simulator_call, nsim=nrow(parsTable))
  simulator_call <- paste(simulator_call, out_parsTableFile)
  file.copy("RNG_state_0000.bin","RNG_state_0000.input.bin",overwrite = TRUE) # save for reproducibility of any problem.
  chktime <- system.time(chk <- system(simulator_call, ignore.stdout = ! verbose.))
  if (chk) {
    warning(paste("Problem with simulator call, system return code=",chk))
    # https://github.com/diyabc/diyabc for doc
    if (chktime["elapsed"]<2) {
      bla <- system(simulator_call, ignore.stdout = FALSE) # simulator output in R console
      str(bla)
    }
  }
  reftable_stats <- read.table(paste0(exe_path,"statfile.txt"), #FIXME hard file name
                               colClasses = c("integer",rep("numeric",length(statNames))))
  if (nrow(reftable_stats) != nrow(parsTable)) {
    if (chktime["elapsed"]<160) { # ___F I X M E___
      warning("Simulator returned fewer simulations than requested (some out of allowed parameter range?)")
      message("Trying to identify suspect parameter values (may take time, simulations are run again):")
      # system(simulator_call, ignore.stdout = FALSE) # simulator output in R console
      info <- system(simulator_call, intern=TRUE)
      message("Suspect parameter values:")
      suspect_parms <- grep("paramvar :", info)
      head(info[suspect_parms])
      stop("Simulator returned fewer simulations than requested.")
      # and to use simulator_call in a console, remove the '\' ...
    }
    
  }
  reftable_stats <- reftable_stats[,-1]
  
  colnames(reftable_stats) <- statNames
  
  setwd(oripath)
  #### conversion to row of reftable in the Infusion style:
  reftable_stats
}

set_seeds <- function(seed, diy_seed=seed, simulator_call) {
  set.seed(seed)
  #
  if (.Platform$OS.type=="windows") {
    exe_pattern <- "diyabc.*"
  } else exe_pattern <- "./diyabc.*" 
  call_splits <- strsplit(simulator_call, split=" +")[[1]]
  exe <- na.omit(stringr::str_extract(call_splits, pattern=exe_pattern))
  diyabc_cores <- as.integer(call_splits[which(call_splits=="-t")+1L])
  simulator_call <- paste0(exe,' -p ./ -n \"t:',diyabc_cores,';c:1;s:',diy_seed,';f:\" ')
  chk <- system(simulator_call, ignore.stdout = FALSE) #=TRUE to see the diyabc output
  if (chk) stop("Problem setting the diyabc seed.")
}
