args <- commandArgs(TRUE)
arg <- function(key,default){a <- grep(paste0("^--",key,"="),args,value=TRUE);if(length(a))sub(paste0("^--",key,"="),"",a[1])else default}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_event_runtime <- function(root,src) {
  for(pkg in c("glmnet","ranger")) {
    okay <- FALSE
    for(attempt in 1:3)if(isTRUE(tryCatch(requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))){okay <- TRUE;break}
    if(!okay)stop("Package unavailable in current process: ",pkg)
  }
  for(f in c("learners.R","learners_lib.R","tests.R"))sys.source(file.path(src,f),envir=.GlobalEnv)
  for(f in c("panel_auc.dll","garch_filter.dll"))dyn.load(file.path(root,"native",f))
  for(f in c("panel_calibration.R","garch_innovations.R"))sys.source(file.path(root,f),envir=.GlobalEnv)
  install_panel_knn();invisible(TRUE)
}
event_seed <- function(job,kind)as.integer(((12000000+(job-1)*100+kind)*65537+987654321)%%2147483648)
event_job <- function(i,input,config) {
  spec <- input$modules[input$modules$job==i,];path <- file.path(config$out,"jobs",sprintf("event_m%02d.rds",i))
  if(file.exists(path))return(invisible(TRUE))
  raw <- read.csv(file.path(input$source,spec$returns_file),check.names=FALSE)
  dates <- as.Date(raw$date);Xraw <- as.matrix(raw[,-1,drop=FALSE])
  baseline <- dates>=as.Date(spec$baseline_start)&dates<=as.Date(spec$baseline_end)
  filters <- lapply(seq_len(ncol(Xraw)),function(j)innovation_garch_fit(Xraw[,j],baseline))
  Z <- do.call(cbind,lapply(filters,`[[`,"z"));H <- do.call(cbind,lapply(filters,`[[`,"h"))
  parameters <- do.call(rbind,lapply(filters,`[[`,"parameters"));rownames(parameters) <- colnames(Xraw)
  phases <- input$phases[input$phases$event==spec$event,]
  sp <- vector("list",3L);fit_dates <- test_dates <- vector("list",3L)
  for(g in 1:3) {
    f <- phases[phases$phase==g&phases$split=="classifier_fit","date"]
    t <- phases[phases$phase==g&phases$split=="test","date"]
    fi <- match(as.Date(f),dates);ti <- match(as.Date(t),dates);stopifnot(!anyNA(fi),!anyNA(ti))
    fi <- fi[apply(is.finite(Z[fi,,drop=FALSE]),1,all)];ti <- ti[apply(is.finite(Z[ti,,drop=FALSE]),1,all)]
    if(min(length(fi),length(ti))<12L)stop("Insufficient prespecified phase data: ",i," phase ",g)
    sp[[g]] <- list(fit=Z[fi,,drop=FALSE],test=Z[ti,,drop=FALSE]);fit_dates[[g]] <- as.character(dates[fi]);test_dates[[g]] <- as.character(dates[ti])
  }
  warnings_seen <- character();set.seed(event_seed(i,1L))
  sc <- withCallingHandlers(fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf")),
    warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  S <- do.call(cbind,sc$Sovr);targets <- rep(1:3,4L);n <- vapply(sp,function(s)nrow(s$test),integer(1))
  st <- standardize_fit(do.call(rbind,lapply(sp,`[[`,"fit")));X <- apply_std(do.call(rbind,lapply(sp,`[[`,"test")),st)
  set.seed(event_seed(i,2L));Lall <- perm_labels(sc$L,config$B);Tm <- ovr_T(S,targets,Lall);nb <- ncol(Lall)
  G <- group_indicator(Lall,3L);D2 <- sqdist(X);D <- sqrt(D2);h2 <- median(D2[upper.tri(D2)]);stopifnot(h2>0)
  K <- exp(-D2/(2*h2));E <- mst_edges(D)
  statistics <- list("MCA-kDET"=apply(Tm,2,max),DISCO=disco_stat(D,G,n,3L,nb),MMD=mmd_stat(K,G,n,3L,nb),
    MST=colSums(Lall[E[,1],,drop=FALSE]==Lall[E[,2],,drop=FALSE]))
  rows <- lapply(names(statistics),function(method) {
    v <- statistics[[method]]
    data.frame(draw=1L,question=paste(spec$event,spec$module,sep="/"),stride=1L,event=spec$event,module=spec$module,
      method=method,d=ncol(Xraw),B=config$B,p_exploratory=pval_perm(v),statistic=v[1],
      critical95=as.numeric(quantile(v[-1],.95,type=1)),test_n=paste(n,collapse="/"),
      fitting_n=paste(vapply(sp,function(s)nrow(s$fit),integer(1)),collapse="/"),
      learner_warning_count=length(warnings_seen),fit_nonconvergence_flags=sum(unlist(lapply(sc$info_ovr,`[[`,"nonconv"))),
      filter_near_boundary=sum(parameters[,"near_boundary"]),status="EXPLORATORY_EVENT_SCREEN_UNVALIDATED")
  })
  record <- data.frame(time=seq_len(nrow(X)),population=sc$L,date=unlist(test_dates))
  audit <- list(records=record,statistics=statistics,h2=h2,S=S,Tm=Tm,Lall=Lall,X=X,targets=targets,coordinate_scale=st)
  value <- list(rows=do.call(rbind,rows),audits=list("1"=audit),draw=1L,q=i,parameters=parameters,
    raw=Xraw,Z=Z,H=H,dates=as.character(dates),baseline=baseline,fitting_dates=fit_dates,testing_dates=test_dates,
    warnings=warnings_seen)
  saveRDS(value,paste0(path,".part"));stopifnot(file.rename(paste0(path,".part"),path))
  cat(sprintf("Event module completed: %s %s\n",spec$event,spec$module));flush.console();invisible(TRUE)
}
load_event_runtime(root,src)
source_folder <- file.path(root,"market_events")
config <- list(B=as.integer(arg("B","4999")),out=arg("out",source_folder))
cores <- as.integer(arg("cores","30"));jobs <- as.integer(strsplit(arg("jobs","1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"),",",fixed=TRUE)[[1]])
dir.create(file.path(config$out,"jobs"),recursive=TRUE,showWarnings=FALSE)
input <- list(source=source_folder,modules=read.csv(file.path(source_folder,"modules_index.csv"),fileEncoding="UTF-8"),
              phases=read.csv(file.path(source_folder,"phase_calendar.csv"),fileEncoding="UTF-8"))
stopifnot(all(jobs%in%input$modules$job))
paths <- c(file.path(root,c("run_event_analysis.R","market_events/protocol.json","market_events/modules_index.csv","market_events/phase_calendar.csv",
  "garch_innovations.R","panel_calibration.R","native/panel_auc.dll","native/garch_filter.dll")),
  file.path(source_folder,input$modules$returns_file),file.path(src,c("learners.R","learners_lib.R","tests.R")))
fingerprint <- list(config=config,jobs=jobs,paths=paths,hashes=tools::md5sum(paths))
fp <- file.path(config$out,"fingerprint.rds");if(file.exists(fp))stopifnot(identical(fingerprint,readRDS(fp)))else saveRDS(fingerprint,fp)
writeLines(capture.output(sessionInfo()),file.path(config$out,"sessionInfo.txt"))
if(cores>1L) {
  cl <- parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(config$out,"workers.log"))
  parallel::clusterCall(cl,load_event_runtime,root,src);parallel::clusterExport(cl,c("event_job","event_seed"),envir=.GlobalEnv)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,function(i,input,config)event_job(i,input,config),input,config)),finally=parallel::stopCluster(cl))
} else invisible(lapply(jobs,event_job,input=input,config=config))
cat(sprintf("%d planned event-module jobs completed. Report/verification use an independent Python parser.\n",length(jobs)))
