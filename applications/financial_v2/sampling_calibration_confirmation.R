# One stock per time block, original frozen MCA permutation engine.
SAMPLING_SCENARIOS <- c(INNOVATION_SCENARIOS,"null_equal_mixtures")
sampling_seed <- function(case,r,kind=0L) {
  stopifnot(case<=9L,r<=10000L,kind<100L)
  as.integer((((case-1)*1000000+(r-1)*100+kind)*65537+987654321)%%2147483648)
}

sampling_generate <- function(scenario,n,seed) {
  if(scenario!="null_equal_mixtures") return(innovation_generate(scenario,n,seed))
  set.seed(seed);J <- 15L;burn <- 3000L;N <- n+burn;loading <- seq(.15,.65,length.out=J)
  sector <- rep(1:3,each=5)
  E <- matrix(rnorm(N*J),N,J)*rep(sqrt(1-loading^2-.4^2),each=N)+
       outer(rnorm(N),loading)+.4*matrix(rnorm(N*3),N,3)[,sector]
  component <- c(1,2,3,4,5,3,5,1,2,4,5,4,2,1,3)
  chi5 <- rchisq(N,5);chi8 <- rchisq(N,8)
  for(j in 1:J) {
    type <- component[j]
    if(type==2) E[,j] <- E[,j]*sqrt(3/chi5)
    if(type==3) E[,j] <- E[,j]*sqrt(6/chi8)
    if(type==4) {delta <- .8;E[,j] <- (delta*abs(E[,j])+sqrt(1-delta^2)*rnorm(N)-delta*sqrt(2/pi))/sqrt(1-2*delta^2/pi)}
    if(type==5) E[,j] <- (E[,j]+sample(c(-1,1),N,replace=TRUE))/sqrt(2)
  }
  # Every sector has exactly the same five component laws, in different orders.
  persistence <- seq(.9,.99,length.out=J);alpha <- seq(.06,.14,length.out=J)
  beta <- persistence-alpha;variance <- seq(.6,1.8,length.out=J)^2;omega <- (1-persistence)*variance;mu <- seq(-.1,.1,length.out=J)
  X <- matrix(0,N,J);h <- previous_x2 <- variance
  for(t in seq_len(N)) {h <- omega+alpha*previous_x2+beta*h;x <- sqrt(h)*E[t,];X[t,] <- mu+x;previous_x2 <- x^2}
  ix <- seq.int(burn+1L,N)
  list(X=X[ix,,drop=FALSE],E=E[ix,,drop=FALSE],truth=data.frame(component=component))
}

sampling_records <- function(valid,q,stride,seed) {
  set.seed(seed);k <- if(q<=3L) 5L else 3L
  offset <- if(stride==1L) 0L else sample.int(stride,1L)-1L
  eligible <- which(valid & ((seq_along(valid)-1L)%%stride==offset))
  n <- length(eligible)%/%k
  if(n<4L) stop("Too few observations per population")
  time <- sort(sample(eligible,n*k,replace=FALSE))
  population <- sample(rep(seq_len(k),each=n))
  stock <- if(q<=3L) (q-1L)*5L+population else (population-1L)*5L+sample.int(5L,length(time),replace=TRUE)
  stopifnot(!anyDuplicated(time),all(tabulate(population,k)==n),all(diff(time)>=stride))
  data.frame(time=time,stock=stock,population=population)
}

sampling_matrices <- function(A,records,k) lapply(seq_len(k),function(g) {
  rows <- which(records$population==g)
  do.call(rbind,lapply(rows,function(i) matrix(A[records$time[i],records$stock[i],],nrow=1L)))
})

sampling_fit <- function(fitting,testing,fit_records,test_valid,q,seed) {
  k <- if(q<=3L) 5L else 3L;d <- dim(fitting)[3L]
  members <- if(q<=3L) seq.int((q-1L)*5L+1L,q*5L) else 1:15
  groups <- if(q<=3L) 1:5 else rep(1:3,each=5L)
  fits <- sampling_matrices(fitting,fit_records,k)
  tests <- lapply(seq_len(k),function(g) do.call(rbind,lapply(members[groups==g],function(j)
    matrix(testing[test_valid,j,,drop=FALSE],ncol=d))))
  sp <- lapply(seq_len(k),function(g) list(fit=fits[[g]],test=tests[[g]]))
  set.seed(seed);warnings_seen <- character()
  sc <- withCallingHandlers(fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf")),
    warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  list(S=do.call(cbind,sc$Sovr),targets=rep(seq_len(k),4L),labels=sc$L,k=k,sp=sp,
       time=rep(which(test_valid),length(members)),stock=rep(members,each=sum(test_valid)),
       warning_count=length(unique(warnings_seen)),nonconv=sum(unlist(lapply(sc$info_ovr,function(o)o$nonconv))))
}

sampling_rep <- function(scenario,r,config,input,keep=FALSE) {
  case <- match(scenario,SAMPLING_SCENARIOS)
  generated <- sampling_generate(scenario,nrow(input$valid),sampling_seed(case,r,0L))
  generated$X[!input$valid] <- NA_real_;generated$E[!input$valid] <- NA_real_
  fl <- lapply(1:15,function(j) innovation_garch_fit(generated$X[,j],input$fitting_days))
  Z <- do.call(cbind,lapply(fl,`[[`,"z"))
  qnames <- c("within_bank","within_pharma","within_machinery","between_sectors")
  schedules <- lapply(1:4,function(q) list(fit=sampling_records(input$masks$fit,q,1L,sampling_seed(case,r,10L+q)),
    test=lapply(seq_along(config$strides),function(h) sampling_records(input$masks$test,q,config$strides[h],sampling_seed(case,r,20L+(q-1L)*3L+h)))))
  output <- audit <- list()
  for(phase in "estimated") {
    M <- if(phase=="oracle") generated$E else Z
    fitting <- innovation_array(M,input$fitting_days,config$d);testing <- innovation_array(M,input$testing_days,config$d)
    rows <- qa <- list()
    for(q in 1:4) {
      model_seed <- sampling_seed(case,r,50L+q)
      scored <- sampling_fit(fitting,testing,schedules[[q]]$fit,input$masks$test,q,model_seed)
      for(h in seq_along(config$strides)) {
        rec <- schedules[[q]]$test[[h]]
        ix <- match(paste(rec$time,rec$stock),paste(scored$time,scored$stock))
        stopifnot(!anyNA(ix),all(scored$labels[ix]==rec$population))
        S <- scored$S[ix,,drop=FALSE]
        set.seed(sampling_seed(case,r,60L+(q-1L)*3L+h))
        Lall <- perm_labels(as.integer(rec$population),config$B)
        Tm <- ovr_T(S,scored$targets,Lall);stat <- apply(Tm,2L,max)
        p <- pval_perm(stat)
        # Refit on the actual sparse test matrices once to verify the full-score shortcut.
        shortcut_error <- NA_real_
        if(keep && h==1L) {
          direct_tests <- sampling_matrices(testing,rec,scored$k)
          sp <- lapply(1:scored$k,function(g) list(fit=scored$sp[[g]]$fit,test=direct_tests[[g]]))
          set.seed(model_seed);direct <- fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf"))
          order <- unlist(lapply(1:scored$k,function(g) which(rec$population==g)))
          shortcut_error <- max(abs(do.call(cbind,direct$Sovr)-S[order,,drop=FALSE]))
          stopifnot(shortcut_error<1e-12)
          check <- ovr_T(do.call(cbind,direct$Sovr),scored$targets,Lall[order,,drop=FALSE])
          stopifnot(max(abs(check-Tm))<1e-10)
        }
        rows[[length(rows)+1L]] <- data.frame(scenario=scenario,rep=r,seed=sampling_seed(case,r,0L),d=config$d,block_length=config$strides[h],
          B=config$B,question=qnames[q],p=p,statistic=stat[1L],critical95=as.numeric(quantile(stat[-1L],.95,type=1)),
          fitting_rows=nrow(schedules[[q]]$fit),testing_rows=nrow(rec),population_test_n=nrow(rec)/scored$k,
          fitting_complete_blocks=sum(input$masks$fit),test_complete_blocks=sum(input$masks$test),
          learner_warning_count=scored$warning_count,fit_nonconvergence_flags=scored$nonconv,
          filter_near_boundary=sum(vapply(fl,function(o)o$parameters[["near_boundary"]],numeric(1))),
          sparse_shortcut_max_score_error=shortcut_error)
        if(keep) qa[[paste(q,h,sep="_")]] <- list(S=S,targets=scored$targets,Lall=Lall,Tm=Tm,
          fitting_records=schedules[[q]]$fit,testing_records=rec,shortcut_error=shortcut_error)
      }
    }
    tab <- do.call(rbind,rows);tab$p_holm <- NA_real_
    for(stride in config$strides) {ix <- which(tab$block_length==stride);tab$p_holm[ix] <- p.adjust(tab$p[ix],"holm")}
    output[[phase]] <- tab
    if(keep) audit[[phase]] <- list(scenario=scenario,rep=r,qa=qa)
  }
  list(rows=output,audit=audit)
}

sampling_job <- function(job,config,input) {
  filename <- sprintf("%s__r%05d-%05d.csv",job$scenario,job$start,job$end)
  paths <- setNames(file.path(config$outdir,"estimated","chunks",filename),"estimated")
  if(all(file.exists(paths))) return(invisible(TRUE))
  rows <- list(estimated=list())
  for(r in seq.int(job$start,job$end)) {
    value <- NULL
    for(a in 1:3) {value <- tryCatch(sampling_rep(job$scenario,r,config,input,keep=(r==1L || r==101L)),error=function(e)e);if(!inherits(value,"error"))break}
    if(inherits(value,"error")) stop(job$scenario," r=",r,": ",conditionMessage(value))
    for(phase in names(rows)) {
      rows[[phase]][[length(rows[[phase]])+1L]] <- value$rows[[phase]]
      if(r==1L || r==101L) saveRDS(value$audit[[phase]],file.path(config$outdir,phase,"audits",sprintf("%s_r%05d.rds",job$scenario,r)))
    }
  }
  for(phase in names(rows)) if(!file.exists(paths[[phase]])) {
    write.csv(do.call(rbind,rows[[phase]]),paste0(paths[[phase]],".part"),row.names=FALSE)
    stopifnot(file.rename(paste0(paths[[phase]],".part"),paths[[phase]]))
  }
  cat(sprintf("Completed one-stock-per-time validation %s reps %d-%d (oracle + estimated)\n",job$scenario,job$start,job$end));flush.console()
  invisible(TRUE)
}
