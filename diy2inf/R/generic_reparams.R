rescaleLoUp <- function(caLOWER,caUPPER, fac) {
  DGP <- caLOWER
  for(it in seq_along(caLOWER)) {
    if (caUPPER[it]>500L) {
      if (caLOWER[it]==0) {
        caLOWER[it] <- log10(1+caLOWER[it])
        caUPPER[it] <- log10(1+caUPPER[it])
        # quotes, whether ` or ', cannot be used because they are not consistently handled (cf as.list(<DGP>) -> doubles the quotes) 
        names(caLOWER)[it] <- names(caUPPER)[it] <- make.names(paste0("log1p(",names(caLOWER)[it],")"))
      } else if (caLOWER[it]>0) {
        caLOWER[it] <- log10(caLOWER[it])
        caUPPER[it] <- log10(caUPPER[it])
        names(caLOWER)[it] <- names(caUPPER)[it] <- make.names(paste0("log(",names(caLOWER)[it],")"))
      }
    }
  }
  DGP <- caLOWER+fac*(caUPPER-caLOWER)
  names(DGP) <- names(caLOWER)
  list(scaLOWER=caLOWER,scaUPPER=caUPPER, DGP=DGP)
} 


# canonize <- function(object, corr=TRUE) reparametrize(object, to=c("Ns","t","dbn","Nbn"), corr=corr)

get_canonizefn <- function(caNames) {
  expr <- paste0("function(object, corr=TRUE) reparametrize(object, to=c('",
                 paste0(caNames, collapse="', '"), 
                 "'), corr=corr)")
  eval(parse(text=expr))
}

############################################

# transformation used on DGparlist 
# and in the primitive workflow building synthetic reftable
#need not be defined on all R^4 but only over the image of reparametrize()
# Maybe largely obsolete now
compose <- function(object, corr=TRUE, to=names(scaLOWER)) {
  
  if ( ! is.null(dim(object))) {
    varnames <- colnames(object)
  } else varnames <- names(object)
  
  if (setequal(varnames,to)) return(object)
  
  for (st in intersect(varnames,to)) assign(st, object[[st]])
  newvars <- setdiff(to, varnames)
  
  logs <- grep("log\\.", newvars)
  log1ps <- grep("log1p\\.", newvars)
  namesplits <- strsplit(newvars, "\\.")
  for (id in logs) {
    namesplit <- namesplits[[id]]
    if (namesplit[2] %in% to) {
      if (id %in% log1ps) {
        assign(namesplit[2], 10^(object[[newvars[id] ]])-1)
      } else  assign(namesplit[2], 10^(object[[newvars[id] ]]))
    }
  }
  
  resu <- mget(to, environment())
  if (inherits(object,"data.frame")) {
    resu <- data.frame(resu)
    #if ( ! is.null(LOWER <- attr(object,"LOWER"))) attr(resu,"LOWER") <- compose(LOWER)
    #if ( ! is.null(UPPER <- attr(object,"UPPER"))) attr(resu,"UPPER") <- compose(UPPER)
    return(resu)
  } else if (is.numeric(object)) {
    unlist(resu)
  } else resu
}

.is_evaluated <- get(".is_evaluated", asNamespace("spaMM"))


# redefine reparametrize and canonize locally
# canonize will be used in the diyabc wrappers
reparametrize <- function(object, to=names(scaLOWER), fixedPars=NULL, corr=TRUE) { # But corr is ignored her
  thisenv <- environment()
  is_evaluated <- function(name) {exists(name, envir=thisenv) && .is_evaluated(name, thisenv)}
  
  # AD HOC PROMISES
  delayedAssign("t12", {
    if (is_evaluated("log.t12.")) {
      10^log.t12.
    } else if (is_evaluated("log1p.t12.")) {
      10^log1p.t12.-1
    } else if (is_evaluated("t1") && is_evaluated("t2")) {
      t2-t1
    } else "t12 missing"
  })
  delayedAssign("t23", {
    if (is_evaluated("log.t23.")) {
      10^log.t23.
    } else if (is_evaluated("log1p.t23.")) {
      10^log1p.t23.-1
    } else if (is_evaluated("t2") && is_evaluated("t3")) {
      t3-t2
    } else "t23 missing"
  })
  delayedAssign("t34", {
    if (is_evaluated("log.t34.")) {
      10^log.t34.
    } else if (is_evaluated("log1p.t34.")) {
      10^log1p.t34.-1
    } else if (is_evaluated("t3") && is_evaluated("t4")) {
      t4-t3
    } else "t34 missing"
  })
  delayedAssign("t4", {
    if (is_evaluated("log.t4.")) {
      10^log.t4.
    } else if (is.character(t34)) { # ie t34 evaluated to "t34 missing"
      t3+dt4
    } else t3+t34
  })
  delayedAssign("t3", {
    if (is_evaluated("log.t3.")) {
      10^log.t3.
    } else if (is.character(t23)) { # as for t4...
      t2+dt3
    } else t2+t23
  })
  delayedAssign("t2", {
    if (is_evaluated("log.t2.")) {
      10^log.t2.
    } else if (is.character(t12)) { # as for t4...
      t1+dt2
    } else t1+t12
  })
  delayedAssign("dt2", {
    if (is_evaluated("log.dt2.")) {
      10^log.dt2.
    } else t2-t1
  })
  delayedAssign("dt3", {
    if (is_evaluated("log.dt3.")) {
      10^log.dt3.
    } else t3-t2
  })
  delayedAssign("dt4", {
    if (is_evaluated("log.dt4.")) {
      10^log.dt4.
    } else t4-t3
  })
  delayedAssign("log.dt2.", log10(dt2))
  delayedAssign("log.dt3.", log10(dt3))
  delayedAssign("log.dt4.", log10(dt4))
  delayedAssign(
    "t2md3", 
    if (is.character(t2)) {
      return(t2)
    } else if (is.character(d3)) {
      return(d3)
    } else t2-d3
  )
  delayedAssign("t2md4", t2-d4)
  delayedAssign("t3md34", t3-d34)

  delayedAssign("d3", 
                if ("d3_Nbn3" %in% varnames) {
                  d3_Nbn3 * Nbn3
                } else t2-t2md3
  )
  delayedAssign("d4", 
                if ("d4_Nbn4" %in% varnames) {
                  d4_Nbn4 * Nbn4
                } else t2-t2md4)
  delayedAssign("d34", 
                if ("d34_Nbn34" %in% varnames) {
                  d34_Nbn34 * Nbn34
                } else t3-t3md34)
  delayedAssign("Na", {
    if (is_evaluated("log.Na.")) {
      10^log.Na.
    } else if (is_evaluated("log1p.Na.")) {
      10^log1p.Na.-1
    } else if (is_evaluated("Na_N2") && is_evaluated("N2")) {
      Na_N2 * N2
    } else "Na missing"
  })
  
  
  delayedAssign("d3_Nbn3", d3/Nbn3)
  delayedAssign("d4_Nbn4", d4/Nbn4)
  delayedAssign("d34_Nbn34", d34/Nbn34)
  delayedAssign("Na_N2", Na/N2)
  
  delayedAssign("logMu", {
    if (is_evaluated("MEANMU")) {
      log10(MEANMU)
    } else "logMu missing"
  })
  
  
  for (st in names(fixedPars)) assign(st, fixedPars[[st]])
  
  #### IDENTIFY newvars
  if ( ! is.null(dim(object))) {
    varnames <- colnames(object)
  } else varnames <- names(object)
  if (setequal(varnames,to)) return(object)
  
  for (st in varnames) assign(st, object[[st]])
  newvars <- setdiff(to, varnames)
  ####
  
  logThs <- newvars[grep("logTh[0-9]", newvars)]
  for (it in seq_along(logThs)) {
    newvar <- logThs[it]
    if ( ! exists(newvar, where=thisenv)) { 
      idx <- strsplit(newvar, "logTh")[[1]][[2]]  
      Nval <- get(paste0("N",idx))
      MEANMU <- get("MEANMU")
      assign(newvar, log10(Nval*MEANMU))
    }
  }
  
  #### WHATEVER IS NOT HANDLED BY THE AD-HOC PROMISES:
  tologs <- grep("log\\.", newvars)
  tolog1ps <- grep("log1p\\.", newvars)
  for (it in seq_along(newvars)) {
    newvar <- newvars[it]
    if ( ! exists(newvar, where=thisenv)) { # testing existence even as unevaluated promise
      namesplit <- strsplit(newvar, "\\.")[[1]] # split at the first '.'
      if (namesplit[1]=="log") { # from X to log10(X)
        assign(newvar, log10(get(namesplit[2])))
      } else if (namesplit[1]=="log1p") { # from X to log10(X)
        assign(newvar, log10(1+get(namesplit[2])))
      } else { # from log10(X) or log10..1p(X) to X 
        chkname <- make.names(paste0("log(",newvar,")"))
        if ( ! exists(chkname, where=thisenv)) {
          chkname <- make.names(paste0("log1p(",newvar,")"))
          if (exists(chkname, where=thisenv)) {
            assign(newvar, 10^get(chkname) -1)
          } else assign(newvar,"is missing")
        } else  assign(newvar, 10^get(chkname))
        # ____F I X M E____ that 'duplicates' the compose() code
      }
    }
  }

  
  
  
  #### 
  resu <- mget(to, environment())
  if (inherits(object,"data.frame")) {
    resu <- data.frame(resu)
    if ( ! is.null(LOWER <- attr(object,"LOWER"))) attr(resu,"LOWER") <- reparametrize(LOWER, to=to)
    if ( ! is.null(UPPER <- attr(object,"UPPER"))) attr(resu,"UPPER") <- reparametrize(UPPER, to=to)
    return(resu)
  } else if (is.numeric(object)) {
    unlist(resu)
  } else resu
}




if (FALSE) {
  
  
  v <- unlist(scaDGparlist)
  
  v <- c(logNs=log10(2),t=200,dbn=1,logNbn=log10(2))
  
  (v <- c(logNs=log10(2),t=200,dbn=1,logNbn=log10(2)))
  
  (v <- c(logNs=1, t=940, dbn=4, logNbn=log10(1.14))) # 1.14 near the lower allwed bound
  
  diff(range(v-unlist(reparametrize(compose(v)))))
  diff(range(compose(v)-unlist(compose(reparametrize(compose(v))))))
  
}


