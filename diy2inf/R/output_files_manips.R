if (FALSE) {
  merge_genotoul_jobs <- function(upto, from=1L) {
    seqft <- seq(from=from,to=upto)
    merger <- vector("list",length(seqft))
    it <- 0L
    for (jobit in seqft) {
      it <- it+1L
      jpath <- paste0("job_",jobit,"/")
      summaries_list <- NULL # safety
      load(paste0(jpath, "summaries_list.",jobit,".rda")) # loads 'summaries_list'
      merger[[it]] <- summaries_list[[1]]
    }
    names(merger) <- seqft
    merger
  }
  
}



df_extract <- function(df, what, finalFUN=NULL) {
  for (st in what) {
    df <- df[,st]
    if (is.list(df)) df <- do.call(rbind,df)
  }
  if ( ! is.null(finalFUN)) {
    # for data frames input at least, lowest level elements loses their attributes: 
    # the output data frame column keep the attr of the last bound dataframe
    # eg, in str(df_extract(ecdfs_df,c("SLRTs",focalpar,"basicLRT")))
    # the last step is to bind the "basicLRT" 1-row data.frames , 
    # whose elements have attribute "solution" which are lost in the process.
    # Hence implementation of finalFUN for use as
    # df_extract(ecdfs_df,c("SLRTs",focalpar), finalFUN=function(v) attr(v$chi2_LR,"solution"))
    # Used for P_4from17 in particular
    df <- lapply(df,finalFUN)
    df <- do.call(rbind,df)
  }
  df
}

merge_jobs <- function(upto, from=1L, prefix="resu_") {
  seqft <- seq(from=from,to=upto)
  merger <- vector("list",length(seqft))
  it <- 0L
  for (jobit in seqft) {
    it <- it+1L
    # jpath <- paste0("job_",jobit,"/")
    summaries_list <- NULL # safety
    try(chk <- load(paste0(prefix,jobit,".rda"))) # loads 'summaries_list'
    if (inherits(chk,"try-error")) {
      merger[it] <- list(NULL)
    } else if (prefix=="summaries_list.") { # synthetic/
      merger[it] <- list(summaries_list[[1]])
    } else merger[[it]] <- get(chk[1]) # old code, probably not the best
    # print(df_extract(do.call(rbind, merger[[it]]$SLRTs),"basicLRT"))
  }
  names(merger) <- seqft
  invisible(merger)
}

# summaries_list <- merge_jobs(50)

save_merge_jobs <- function(upto, from=1L, prefix="resu_") {
  summaries_list <- merge_jobs(upto,from, prefix=prefix)
  lnames <- names(summaries_list)
  savename <- paste0("summaries_list.",lnames[1],"_",tail(lnames,1), 
                     ".v",packageVersion("Infusion"),".rda")
  save(summaries_list, file=savename)
  invisible(summaries_list)
}

# resus <- save_merge_jobs(50)
# resus <- save_merge_jobs(100, 50)
# resus <- save_merge_jobs(150, 100)
# resus <- save_merge_jobs(200, 150)

