# Oracle and estimated-filter pipelines share each simulated panel and seed.
INNOVATION_SCENARIOS <- c("null_iid","null_garch95","null_garch99","null_heterogeneous",
                         "null_garch95_t5","null_garch999","alt_sector_shape","alt_one_bank")

innovation_generate <- function(scenario,n,seed) {
  set.seed(seed);J <- 15L;burn <- if(scenario=="null_garch999") 20000L else 3000L;N <- n+burn
  loading <- seq(.15,.65,length.out=J);sector <- rep(1:3,each=5L)
  E <- matrix(rnorm(N*J),N,J)*rep(sqrt(1-loading^2-.4^2),each=N)
  E <- E+outer(rnorm(N),loading)+.4*matrix(rnorm(N*3L),N,3L)[,sector]
  if(scenario=="null_garch95_t5") E <- E*sqrt(3/rchisq(N,5))
  if(scenario=="alt_sector_shape") {
    E[,6:10] <- E[,6:10]*sqrt(3/rchisq(N,5))
    # Mean-zero, unit-variance skew-normal coordinates from correlated normals.
    delta <- .8; extra <- matrix(rnorm(N*5L),N,5L)
    E[,11:15] <- (delta*abs(E[,11:15])+sqrt(1-delta^2)*extra-delta*sqrt(2/pi))/sqrt(1-2*delta^2/pi)
  }
  if(scenario=="alt_one_bank") E[,1L] <- E[,1L]*sqrt(3/rchisq(N,5))
  persistence <- rep(if(scenario=="null_garch99") .99 else if(scenario=="null_garch999") .999 else .95,J)
  alpha <- rep(.1,J);variance <- rep(1,J);mu <- rep(0,J)
  if(scenario=="null_iid") {persistence[] <- 0;alpha[] <- 0}
  if(scenario=="null_heterogeneous") {
    persistence <- seq(.9,.99,length.out=J);alpha <- seq(.06,.14,length.out=J)
    variance <- seq(.6,1.8,length.out=J)^2;mu <- seq(-.1,.1,length.out=J)
  }
  beta <- persistence-alpha;omega <- (1-persistence)*variance
  X <- matrix(0,N,J);h <- previous_x2 <- variance
  for(t in seq_len(N)) {
    h <- omega+alpha*previous_x2+beta*h
    x <- sqrt(h)*E[t,];X[t,] <- mu+x;previous_x2 <- x^2
  }
  index <- seq.int(burn+1L,N)
  list(X=X[index,,drop=FALSE],E=E[index,,drop=FALSE],
       truth=data.frame(stock_index=seq_len(J),true_mean=mu,true_omega=omega,true_alpha=alpha,true_beta=beta,true_persistence=persistence))
}

innovation_array <- function(M,days,d) {
  ids <- which(days);blocks <- length(ids)%/%d
  aperm(array(M[ids[seq_len(blocks*d)],,drop=FALSE],c(d,blocks,ncol(M))),c(2L,3L,1L))
}

innovation_rep <- function(scenario,r,config,input,keep=FALSE) {
  id <- match(scenario,INNOVATION_SCENARIOS);seed <- config$seed0+id*1000000L+r*100L
  generated <- innovation_generate(scenario,nrow(input$valid),seed)
  generated$X[!input$valid] <- NA_real_;generated$E[!input$valid] <- NA_real_
  filters <- lapply(seq_len(15L),function(j) innovation_garch_fit(generated$X[,j],input$fitting_days))
  estimated <- do.call(cbind,lapply(filters,`[[`,"z"))
  stopifnot(identical(dim(estimated),dim(input$valid)),all(is.na(estimated)==!input$valid))
  estimates <- data.frame(scenario=scenario,rep=r,seed=seed,generated$truth,
                          do.call(rbind,lapply(filters,`[[`,"parameters")),row.names=NULL)
  qnames <- c("within_bank","within_pharma","within_machinery","between_sectors")
  results <- audits <- list()
  for(phase in c("oracle","estimated")) {
    M <- if(phase=="oracle") generated$E else estimated
    fit <- innovation_array(M,input$fitting_days,config$d)
    test <- innovation_array(M,input$testing_days,config$d)
    stopifnot(dim(fit)[1L]==length(input$masks$fit),dim(test)[1L]==length(input$masks$test))
    questions <- lapply(seq_len(4L),function(q) {
      members <- if(q<=3L) seq.int((q-1L)*5L+1L,q*5L) else seq_len(15L)
      groups <- if(q<=3L) seq_len(5L) else rep(seq_len(3L),each=5L)
      panel_question(fit,test,input$masks$fit,input$masks$test,members,groups,seed+q)
    })
    rows <- retained <- list()
    for(len in config$lengths) {
      counts <- panel_counts(length(input$masks$test),len,config$B,seed+50L+len)
      calibrated <- lapply(questions,panel_calibrate,counts=counts)
      ps <- vapply(calibrated,`[[`,numeric(1),"p");adjusted <- p.adjust(ps,"holm")
      for(q in seq_len(4L)) rows[[length(rows)+1L]] <- data.frame(
        scenario=scenario,rep=r,seed=seed,d=config$d,block_length=len,B=config$B,question=qnames[q],
        p=ps[q],p_holm=adjusted[q],statistic=calibrated[[q]]$statistic,critical95=calibrated[[q]]$critical95,
        fitting_complete_blocks=sum(input$masks$fit),test_complete_blocks=sum(input$masks$test),
        fitting_rows=questions[[q]]$fit_n,testing_rows=questions[[q]]$test_n,
        learner_warning_count=length(questions[[q]]$warnings),fit_nonconvergence_flags=questions[[q]]$nonconv,
        filter_failed_starts=sum(estimates$failed_starts),filter_near_boundary=sum(estimates$near_boundary))
      if(keep) retained[[as.character(len)]] <- list(counts=counts,calibrated=calibrated)
    }
    results[[phase]] <- do.call(rbind,rows)
    if(keep) audits[[phase]] <- list(questions=questions,retained=retained,scenario=scenario,rep=r,seed=seed)
  }
  list(rows=results,filters=estimates,audit=audits)
}

innovation_job <- function(job,config,input) {
  filename <- sprintf("%s__r%05d-%05d.csv",job$scenario,job$start,job$end)
  targets <- c("oracle","estimated","filters")
  destinations <- setNames(file.path(config$outdir,targets,"chunks",filename),targets)
  if(all(file.exists(destinations))) return(invisible(TRUE))
  rows <- lapply(targets,function(x) list());names(rows) <- targets
  for(r in seq.int(job$start,job$end)) {
    value <- NULL
    for(attempt in 1:3) {
      value <- tryCatch(innovation_rep(job$scenario,r,config,input,keep=r==1L),error=function(e)e)
      if(!inherits(value,"error")) break
    }
    if(inherits(value,"error")) stop(job$scenario," r=",r,": ",conditionMessage(value))
    for(phase in c("oracle","estimated")) {
      rows[[phase]][[length(rows[[phase]])+1L]] <- value$rows[[phase]]
      if(r==1L) {
        path <- file.path(config$outdir,phase,"audits",paste0(job$scenario,"_r00001.rds"))
        if(!file.exists(path)) saveRDS(value$audit[[phase]],path)
      }
    }
    rows$filters[[length(rows$filters)+1L]] <- value$filters
  }
  for(phase in targets) {
    tab <- do.call(rbind,rows[[phase]]);path <- destinations[[phase]]
    if(file.exists(path)) stopifnot(isTRUE(all.equal(read.csv(path),tab,check.attributes=FALSE,tolerance=1e-12))) else {
      write.csv(tab,paste0(path,".part"),row.names=FALSE)
      stopifnot(file.rename(paste0(path,".part"),path))
    }
  }
  cat(sprintf("Completed innovation validation %s reps %d-%d (oracle + re-estimated GARCH)\n",job$scenario,job$start,job$end));flush.console()
  invisible(TRUE)
}
