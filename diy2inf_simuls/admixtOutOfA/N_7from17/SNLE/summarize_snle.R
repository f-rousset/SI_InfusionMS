if (FALSE) { # this should be run on the simulation server
  jobs <- 1:400 # c(1:80,101:170)
  postm <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_means_10K.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  postm <- t(postm)
  save(postm, file="postm_10K.rda")
  
  q1s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_q1s_10K.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q1s <- t(q1s)
  save(q1s, file="q1s_10K.rda")
  
  q2s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_q2s_10K.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q2s <- t(q2s)
  save(q2s, file="q2s_10K.rda")
  
}

if (FALSE) { # this should be run on the simulation server
  jobs <- 1:400 # c(1:80,101:170)
  postm <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_means.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  postm <- t(postm)
  save(postm, file="postm.rda")
  
  q1s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_q1s.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q1s <- t(q1s)
  save(q1s, file="q1s.rda")
  
  q2s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("job",it,"/posterior_q2s.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q2s <- t(q2s)
  save(q2s, file="q2s.rda")
  
}

#setwd("D:/home/francois/travail/stats/Infusionplus/sbi/paper-infusion-SNLE/admixture")
# setwd("D:/home/francois/travail/stats/Infusionplus/MS/SI_InfusionMS/diy2inf_simuls/admixtOutOfA/N_7from17/snle/")
load("D:/home/francois/travail/stats/Infusionplus/diy2inf_simuls/admixtOutOfA/N_7from17/saves4Rmd.rda")

{
  {
    load("postm.rda")
    postm <- as.data.frame(postm)
    colnames(postm) <- names(LOWER)
    print(nrow(postm))
    
    load("q1s.rda")
    q1s <- as.data.frame(q1s)
    colnames(q1s) <- names(LOWER)
    
    load("q2s.rda")
    q2s <- as.data.frame(q2s)
    colnames(q2s) <- names(LOWER)
  }
  
  {
    load("postm_10K.rda")
    postm <- as.data.frame(postm)
    colnames(postm) <- names(LOWER)
    print(nrow(postm))
    
    load("q1s_10K.rda")
    q1s <- as.data.frame(q1s)
    colnames(q1s) <- names(LOWER)
    
    load("q2s_10K.rda")
    q2s <- as.data.frame(q2s)
    colnames(q2s) <- names(LOWER)
    
  }
  if (the_.csv_contained_unscaled_values <- TRUE) { 
    rescale_py_output <- function(df, # typically a data frame
                                  log10.=c("log.N2.", "log.t12.", "log.Na." ), 
                                  log1p10.=c("log1p.t23.", "log1p.t34.") 
    ) {
      df[,log10.] <- log10(df[,log10.])
      df[,log1p10.] <- log10( 1 + df[,log1p10.])
      df
    }
    postm <- rescale_py_output(postm)
    q1s <- rescale_py_output(q1s)
    q2s <- rescale_py_output(q2s)
  }
  save(postm, q1s, q2s, file="SNLE_summaries.rda")
  
}




for (parm in names(LOWER)) {
  m <- (postm[,parm] - scaDGP[parm])/(UPPER[parm] - LOWER[parm])
  cover <- (q1s[,parm]< scaDGP[parm] & q2s[,parm] > scaDGP[parm] )
  print(c(mean(m),sqrt(mean(m^2)), sum(cover)/length(cover)))
}
