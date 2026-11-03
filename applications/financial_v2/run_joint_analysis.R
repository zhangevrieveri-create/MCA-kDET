args <- commandArgs(TRUE)
arg <- function(key,default){a <- grep(paste0("^--",key,"="),args,value=TRUE);if(length(a))sub(paste0("^--",key,"="),"",a[1])else default}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_joint_runtime <- function(root,src) {
  for(pkg in c("glmnet","ranger")) {
    okay <- FALSE
    for(attempt in 1:3)if(isTRUE(tryCatch(requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))){okay <- TRUE;break}
    if(!okay)stop("Package unavailable in this R process: ",pkg)
  }
  for(f in c("learners.R","learners_lib.R","tests.R"))sys.source(file.path(src,f),envir=.GlobalEnv)
  dyn.load(file.path(root,"native/panel_auc.dll"));sys.source(file.path(root,"panel_calibration.R"),envir=.GlobalEnv)
  install_panel_knn();invisible(TRUE)
}
real_seed <- function(draw,kind)as.integer(((9000000+(draw-1)*100+kind)*65537+987654321)%%2147483648)
joint_job <- function(draw,input,config) {
  path <- file.path(config$out,"jobs",sprintf("draw%02d_q4.rds",draw))
  if(file.exists(path))return(invisible(TRUE))
  recs <- lapply(c(1L,2L,4L),function(stride)input$schedule[input$schedule$draw==draw&input$schedule$stride==stride,])
  bag <- do.call(rbind,recs)
  tests <- lapply(1:3,function(g)input$test[[g]][bag$time[bag$population==g],,drop=FALSE])
  sp <- lapply(1:3,function(g)list(fit=input$fit[[g]],test=tests[[g]]))
  set.seed(real_seed(draw,54L));warnings_seen <- character()
  scored <- withCallingHandlers(fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf")),
    warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  Sbag <- do.call(cbind,scored$Sovr);targets <- rep(1:3,4L)
  bagorder <- unlist(lapply(1:3,function(g)which(bag$population==g)))
  st <- standardize_fit(do.call(rbind,input$fit));rows <- aucs <- audits <- list()
  for(h in 1:3) {
    stride <- c(1L,2L,4L)[h];rec <- recs[[h]]
    ix <- match(paste(rec$time,rec$population),paste(bag$time[bagorder],bag$population[bagorder]))
    stopifnot(!anyNA(ix));S <- Sbag[ix,,drop=FALSE]
    rawX <- do.call(rbind,lapply(seq_len(nrow(rec)),function(i)input$test[[rec$population[i]]][rec$time[i],,drop=FALSE]))
    X <- apply_std(rawX,st)
    set.seed(real_seed(draw,69L+h));Lall <- perm_labels(as.integer(rec$population),config$B)
    Tm <- ovr_T(S,targets,Lall);nb <- ncol(Lall);n <- tabulate(rec$population,3L);G <- group_indicator(Lall,3L)
    D2 <- sqdist(X);D <- sqrt(D2);h2 <- median(D2[upper.tri(D2)]);stopifnot(h2>0)
    K <- exp(-D2/(2*h2));E <- mst_edges(D)
    statistics <- list("MCA-kDET"=apply(Tm,2,max),DISCO=disco_stat(D,G,n,3L,nb),MMD=mmd_stat(K,G,n,3L,nb),
      MST=colSums(Lall[E[,1],,drop=FALSE]==Lall[E[,2],,drop=FALSE]))
    if(draw==1L && stride==4L) {
      sparse <- lapply(1:3,function(g)list(fit=input$fit[[g]],test=input$test[[g]][rec$time[rec$population==g],,drop=FALSE]))
      set.seed(real_seed(draw,54L));direct <- do.call(cbind,fit_scores(sparse,"ovr",lib=c("l1","quad","knn","rf"))$Sovr)
      order <- unlist(lapply(1:3,function(g)which(rec$population==g)))
      stopifnot(max(abs(direct-S[order,,drop=FALSE]))<1e-12)
    }
    for(method in names(statistics)) {
      v <- statistics[[method]]
      rows[[length(rows)+1L]] <- data.frame(draw=draw,question="between_sectors",stride=stride,method=method,d=5L,B=config$B,
        p_exploratory=pval_perm(v),statistic=v[1],critical95=as.numeric(quantile(v[-1],.95,type=1)),populations=3L,
        test_n_per_population=n[1],testing_rows=nrow(rec),fitting_rows=sum(vapply(input$fit,nrow,integer(1))),
        learner_warning_count=length(warnings_seen),fit_nonconvergence_flags=sum(unlist(lapply(scored$info_ovr,`[[`,"nonconv"))),
        primary=(draw==1L&&stride==4L),status="EXPLORATORY_NEW_JOINT_TARGET")
    }
    for(c in 1:ncol(S))aucs[[length(aucs)+1L]] <- data.frame(draw=draw,question="between_sectors",stride=stride,
      population=targets[c],learner=rep(c("l1","quad","knn","rf"),each=3L)[c],
      auc=auc_vec(S[,c],rec$population==targets[c]),T_observed=Tm[c,1])
    audits[[as.character(stride)]] <- list(records=rec,statistics=statistics,h2=h2,
      S=if(draw==1L)S else NULL,Tm=if(draw==1L)Tm else NULL,Lall=if(draw==1L)Lall else NULL,
      X=if(draw==1L)X else NULL,targets=targets,coordinate_scale=st)
  }
  value <- list(rows=do.call(rbind,rows),aucs=do.call(rbind,aucs),audits=audits,draw=draw,q=4L)
  saveRDS(value,paste0(path,".part"));stopifnot(file.rename(paste0(path,".part"),path))
  cat(sprintf("Sorted sector-vector analysis completed draw %d\n",draw));flush.console()
  invisible(TRUE)
}
load_joint_runtime(root,src)
config <- list(B=as.integer(arg("B","4999")),out=arg("out",file.path(root,"joint_analysis")))
draws <- as.integer(arg("draws","10"));cores <- as.integer(arg("cores","10"));stopifnot(draws>=1L,draws<=10L)
dir.create(file.path(config$out,"jobs"),recursive=TRUE,showWarnings=FALSE)
features <- read.csv(file.path(root,"joint_analysis/sorted_sector_vectors.csv"))
schedule <- read.csv(file.path(root,"joint_analysis/fixed_schedules.csv"))
fit <- test <- vector("list",3L)
for(g in 1:3) {
  a <- features[features$split=="fit"&features$sector==g,]
  valid <- tolower(as.character(a$complete_all_stocks))=="true"
  fit[[g]] <- as.matrix(a[valid,paste0("x",1:5)])
  a <- features[features$split=="test"&features$sector==g,];test[[g]] <- as.matrix(a[,paste0("x",1:5)])
}
stopifnot(all(vapply(fit,nrow,integer(1))==116L))
input <- list(fit=fit,test=test,schedule=schedule)
fingerprint <- list(config=config,draws=draws,hashes=tools::md5sum(c(file.path(root,c("run_joint_analysis.R","joint_analysis/protocol.json",
  "joint_analysis/sorted_sector_vectors.csv","joint_analysis/fixed_schedules.csv","panel_calibration.R","native/panel_auc.dll")),file.path(src,c("learners.R","learners_lib.R","tests.R")))))
fp <- file.path(config$out,"fingerprint.rds");if(file.exists(fp))stopifnot(identical(fingerprint,readRDS(fp)))else saveRDS(fingerprint,fp)
writeLines(capture.output(sessionInfo()),file.path(config$out,"sessionInfo.txt"))
if(cores>1L) {
  cl <- parallel::makePSOCKcluster(min(cores,draws),outfile=file.path(config$out,"workers.log"))
  parallel::clusterCall(cl,load_joint_runtime,root,src);parallel::clusterExport(cl,c("joint_job","real_seed"),envir=.GlobalEnv)
  tryCatch(invisible(parallel::parLapplyLB(cl,seq_len(draws),function(draw,input,config)joint_job(draw,input,config),input,config)),finally=parallel::stopCluster(cl))
} else invisible(lapply(seq_len(draws),joint_job,input=input,config=config))
saved <- lapply(seq_len(draws),function(draw)readRDS(file.path(config$out,"jobs",sprintf("draw%02d_q4.rds",draw))))
write.csv(do.call(rbind,lapply(saved,`[[`,"rows")),file.path(config$out,"between_results.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(saved,`[[`,"aucs")),file.path(config$out,"learner_auc.csv"),row.names=FALSE)
cat("Exploratory sorted-vector results saved; no validated calibration claim.\n")
