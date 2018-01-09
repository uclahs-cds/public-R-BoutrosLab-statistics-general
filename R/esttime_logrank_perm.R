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

esttime_logrank_perm <- function(N, k, no_cores=1, tform='mins') {
  # Error checks
  check1 <- tform %in% c('secs','mins','hours')
  if (!(check1)) {
    print('Choose a valid tform: secs, mins, or hours')
    stopifnot(check1)
  }
  
  # ----- (2) Functions ------ # 
  
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
  
  # ------- (2) Set-up ------- #
  
  # Load in the microbenchmark package
  # library(microbenchmark,warn.conflicts = F,quietly = T,verbose = F)
  library(rbenchmark,warn.conflicts = F,quietly = T,verbose = F)
  
  # Calculate the number of combinations
  ncombn <- choose(N,k)
  
  # Calculate time needed to create matrices
  time.start <- Sys.time()
  print('Creating all valid permutations')
  # Create a toy combination list
  lst.combn <- combn(N,k,simplify = F)
  # Convert index list to vectors
  vec.combn <- lapply(lst.combn,idx2ones,N=N)
  time.end <- Sys.time()
  # Number of seconds to create combination matrix
  sec.combn <- as.numeric(difftime(time.end,time.start,units='secs'))
  
  # Random survival data
  Tobs <- rexp(N)
  is.event <- rep(1,N)
  So <- matrix(c(Tobs,is.event),ncol=2)
  
  # ----- (3) Use rbenchmark ------ # 
    
  # Use benchmark function for sequence of up to 10K calculations (or lower)
  if (ncombn > 100000) {
    iseq <- seq(100,100000,length.out = 4)
  } else {
    iseq <- seq(1,ncombn,length.out = 5)[-1]
  }
  # Storeage
  sec.store <- rep(NA,length(iseq))
  
  # Loop over and calculate the number of seconds
  for (i in seq_along(iseq)) {
    # Run the benchmark
    temp.bench <- benchmark(
      sapply(vec.combn[1:iseq[i]],lr.chi,So=So),
      replications = 1
    )
    # Number of seconds
    temp.sec <- temp.bench$elapsed
    # Store
    sec.store[i] <- temp.sec
    # Update
    print(sprintf('Simulation %i of %i',i,length(iseq)))
  }
  # Data.frame
  df.sec <- data.frame(ncalc=iseq,seconds=sec.store)
  # Fit a regression model
  mdl <- lm(seconds ~ ncalc,data=df.sec)
  # Make prediction
  pred.secs <- predict(mdl,newdata=data.frame(ncalc=ncombn))
  
  # Use our rough parallelization efficiency:
  nc <- ifelse(no_cores > 6,6,no_cores)
  skedge <- list('1'=1.0000000,'2'=0.5804742,'3'=0.4279533,
       '4'=0.3634352,'5'=0.3190550,'6'=0.3073446)
  para.gain <- skedge[[as.character(nc)]]
  # Update
  pred.secs.parallel <- pred.secs * para.gain
  
  # Combine
  total.secs <- pred.secs.parallel + sec.combn
  # Format
  if (tform=='secs') {
    tt.ret <- total.secs
    print(sprintf('Estimated runtime: %0.1f seconds',tt.ret))
  }
  if (tform=='mins') {
    tt.ret <- total.secs/60
    print(sprintf('Estimated runtime: %0.1f minutes',tt.ret))
  }
  if (tform=='hours') {
    tt.ret <- (total.secs/60)/60
    print(sprintf('Estimated runtime: %0.1f hours',tt.ret))
  }
  # Return
  return(tt.ret)
}

