if (FALSE) { # this block should be run on the simulation server
  jobs <- 1:400 
  postm <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("posterior_means.",it-1L,".csv"),header = FALSE,sep=";",
             colClasses = rep("numeric",7)))
  })
  postm <- t(postm)
  save(postm, file="postm.rda")
  
  q1s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("posterior_q1s.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q1s <- t(q1s)
  save(q1s, file="q1s.rda")
  
  q2s <- sapply(jobs, function(it) {
    unlist(read.csv(paste0("posterior_q2s.",it-1L,".csv"),header = FALSE,sep=";",
                    colClasses = rep("numeric",7)))
  })
  q2s <- t(q2s)
  save(q2s, file="q2s.rda")
  
}

{ # conversion to files stored in SI
  setwd("D:/home/francois/travail/stats/Infusionplus/sbi/paper-infusion-SNLE/covariance")
  # setwd("D:/home/francois/travail/stats/Infusionplus/MS/SI_InfusionMS/diy2inf_simuls/admixtOutOfA/N_7from17/snle/")
  
  # original R variables
  load(file="saves4Rmd.rda") # LOWER,UPPER,cholDGPv
  
  { # reordering R objects as in Python code
    Rord_names <- names(LOWER)
    ij <- matrix(paste("V",gl(5,5,25),gl(5,1,25), sep=""),ncol=5,nrow=5)
    Pyord_names <-  ij[lower.tri(ij,diag = TRUE)]
    LOWER <- LOWER[Pyord_names]
    UPPER <- UPPER[Pyord_names]
    cholDGPv <- cholDGPv[Pyord_names]
  }
  
  {
    load("postm.rda")
    postm <- as.data.frame(postm)
    colnames(postm) <- Pyord_names
    postm <- postm[,Rord_names]
    print(nrow(postm))
    
    load("q1s.rda")
    q1s <- as.data.frame(q1s)
    colnames(q1s) <- Pyord_names
    q1s <- q1s[,Rord_names]
    
    load("q2s.rda")
    q2s <- as.data.frame(q2s)
    colnames(q2s) <- Pyord_names
    q2s <- q2s[,Rord_names]
    
    save(postm, q1s, q2s, file="SNLE_summaries.rda")
  }
}

{ # performance summaries
  # load(file="SNLE_summaries.rda"))
  resu <- vector("list", 15L)
  names(resu) <- Rord_names
  for (parm in Rord_names) {
    m <- (postm[,parm] - cholDGPv[parm])/(UPPER[parm] - LOWER[parm])
    cover <- (q1s[,parm]< cholDGPv[parm] & q2s[,parm] > cholDGPv[parm] )
    print(resu[[parm]] <- c(mean(m),sd(m), sum(cover)/length(cover)))
  }
  resu <- do.call(rbind, resu)
  rownames(resu) <- Rord_names
  colMeans(resu)
  
}

