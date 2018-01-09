# The BoutrosLab.statistics.general package is copyright (c) 2017 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

logrank_perm <- function(So, x, no_cores=1) {
	NSo <- nrow(So)
	pSo <- ncol(So)
	N <- length(x)
	k <- sum(x)
  # Checks
  check1 <- NSo==N
  check2 <- pSo==2
  # Parallelization
  library(doParallel,quietly = T,warn.conflicts = F,verbose = F)
  nc <- detectCores()
  check3 <- (nc-1) >= no_cores
  # ----- Error checks ------ #
  if (!check1) {
    print('Error! The number of observations in So does not match x')
    stopifnot(check1)
  }
  if (!check2) {
    print('Error! So should have two columns!')
    stopifnot(check2)
  }
  if (!check3) {
    print('Error! Use one less core that is available!')
    stopifnot(check3)
  }
  
  # ----- (1) Functions ------ # 
  
  # (i) calculate the log-rank stat
  lr.chi <- function(x,So) {
    fit <- survival:::survdiff.fit(y=So,x=x,rho=0)
    otmp <- fit$observed
    etmp <- fit$expected
    stat <- sum((otmp-etmp)^2)/sum(diag(fit$var))
    return(stat)
  }
  
  # (ii) Assign a vector ones to the locations of an index
  idx2ones <- function(idx,N) {
    zz <- rep(0,N)
    zz[idx] <- 1
    return(zz)
  }
  
  # ----- (2) Create the combination list ------ # 
  print('Creating all valid permutations')
  # Time to make the matrix
  time.start <- Sys.time()
  # Put all in a list
  lst.combn <- combn(N,k,simplify = F)
  # Calculate number of combinations
  ncombn <- length(lst.combn)
  # Total seconds
  time.end <- Sys.time()
  # Number of seconds to create combination matrix
  sec.combn <- as.numeric(difftime(time.end,time.start,units='secs'))
  
  # ----- (3) Calculate the statistics for each ------ # 
  
  # Calculate our baseline statistic
  baseline.stat <- lr.chi(x=x,So=So)
  
  # Get the percentiles
  pcs <- round((seq(1,99)/100)*ncombn)
  
  # No parallelization
  if (no_cores==1) {
    print('Running without parallelization')
    # Create a vector to store the stats in
    stat.vec <- rep(NA,ncombn)
    # Start the clock
    time.old <- Sys.time()
    time.start <- Sys.time()
    # Start the loop
    for (k in 1:ncombn) {
      xx <- idx2ones(lst.combn[[k]],N)
      stat.vec[k] <- lr.chi(x=xx,So=So)
      # Check
      if (any(k==pcs)) {
        # Which percent
        wpcs <- which(k==pcs)
        # Update time and trajectory
        time.new <- Sys.time()
        est.min <- as.numeric(difftime(time.new,time.old,units = 'min'))
        est.remaining <- round(est.min*(100-wpcs),1)
        # Update for user
        print(sprintf('%i%% done; (estimated) %0.1f minutes remaining',wpcs,est.remaining))
        # Restart the clock
        time.old <- time.new
      }
    }
    time.end <- Sys.time()
    est.sec <- as.numeric(difftime(time.end,time.start,units = 'sec'))
    print(sprintf('Permutation test took %i seconds',round(est.sec+sec.combn)))
  # Parallelize
  } else {
    print('Running with parallelization')
    # Create the vectorized version
    vec.combn <- lapply(lst.combn,idx2ones,N=N)
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    # use Lapply
    time.old <- Sys.time()
    stat.vec <- parSapply(cl,vec.combn,lr.chi,So=So)
    stopCluster(cl)
    # stat.vec
    time.new <- Sys.time()
    # Get the total time
    est.sec <- as.numeric(difftime(time.new,time.old,units = 'sec'))
    print(sprintf('Permutation test took %i seconds',round(est.sec+sec.combn)))
  }
  
  # ----- (3) p-value and final info ------ # 
  # Number of statistics that are larger (w/ equality)
  ngreater <- sum(stat.vec > baseline.stat)
  pval <- ngreater/ncombn
  
  # Return
  dat <- data.frame(ntests=ncombn,ngeq=ngreater,pval=pval)
  return(dat)
}


