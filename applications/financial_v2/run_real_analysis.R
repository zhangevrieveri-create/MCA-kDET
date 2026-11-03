args <- commandArgs(TRUE)
arg <- function(key,default){a <- grep(paste0("^--",key,"="),args,value=TRUE);if(length(a))sub(paste0("^--",key,"="),"",a[1])else default}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_real_runtime <- function(root,src) {
  stopifnot(requireNamespace("glmnet",quietly=TRUE),requireNamespace("ranger",quietly=TRUE))
  for(f in c("learners.R","learners_lib.R","tests.R"))sys.source(file.path(src,f),envir=.GlobalEnv)
  dyn.load(file.path(root,"native/panel_auc.dll"))
  for(f in c("panel_calibration.R","innovation_validation.R","sampling_calibration.R","sampling_calibration_fullfit.R"))sys.source(file.path(root,f),envir=.GlobalEnv)
  install_panel_knn()
  invisible(TRUE)
}
real_seed <- function(draw,kind)as.integer(((9000000+(draw-1)*100+kind)*65537+987654321)%%2147483648)
real_job <- function(job,input,config) {
  q <- job$q;draw <- job$draw;k <- if(q<=3L)5L else 3L
  path <- file.path(config$out,"jobs",sprintf("draw%02d_q%d.rds",draw,q))
  if(file.exists(path))return(readRDS(path)$rows)
  model_seed <- real_seed(draw,50L+q)
  scored <- sampling_fit(input$fit,input$test,NULL,input$masks$test,q,model_seed)
  st <- standardize_fit(do.call(rbind,lapply(scored$sp,`[[`,"fit")))
  rows <- aucs <- list();audits <- list()
  for(h in 1:3) {
    stride <- c(1L,2L,4L)[h]
    rec <- sampling_records(input$masks$test,q,stride,real_seed(draw,20L+(q-1L)*3L+h))
    ix <- match(paste(rec$time,rec$stock),paste(scored$time,scored$stock));stopifnot(!anyNA(ix))
    S <- scored$S[ix,,drop=FALSE]
    X <- do.call(rbind,lapply(seq_len(nrow(rec)),function(i)matrix(input$test[rec$time[i],rec$stock[i],],nrow=1)))
    X <- apply_std(X,st)
    set.seed(real_seed(draw,60L+(q-1L)*3L+h));Lall <- perm_labels(as.integer(rec$population),config$B)
    Tm <- ovr_T(S,scored$targets,Lall)
    n <- tabulate(rec$population,k);nb <- ncol(Lall);G <- group_indicator(Lall,k)
    D2 <- sqdist(X);D <- sqrt(D2);h2 <- median(D2[upper.tri(D2)]);stopifnot(h2>0)
    K <- exp(-D2/(2*h2));E <- mst_edges(D)
    statistics <- list("MCA-kDET"=apply(Tm,2,max),DISCO=disco_stat(D,G,n,k,nb),MMD=mmd_stat(K,G,n,k,nb),
                       MST=colSums(Lall[E[,1],,drop=FALSE]==Lall[E[,2],,drop=FALSE]))
    for(method in names(statistics)) {
      value <- statistics[[method]];p <- pval_perm(value)
      rows[[length(rows)+1L]] <- data.frame(draw=draw,question=input$qnames[q],stride=stride,method=method,d=5L,B=config$B,
        p_exploratory=p,statistic=value[1],critical95=as.numeric(quantile(value[-1],.95,type=1)),
        populations=k,test_n_per_population=n[1],testing_rows=nrow(rec),fitting_rows=if(q<=3L)580L else 1740L,
        learner_warning_count=scored$warning_count,fit_nonconvergence_flags=scored$nonconv,
        primary=(draw==1L && stride==4L),status="EXPLORATORY_UNVALIDATED_CALIBRATION")
    }
    for(c in seq_len(ncol(S))) {
      pop <- scored$targets[c];a <- auc_vec(S[,c],rec$population==pop)
      aucs[[length(aucs)+1L]] <- data.frame(draw=draw,question=input$qnames[q],stride=stride,
        population=pop,learner=rep(c("l1","quad","knn","rf"),each=k)[c],auc=a,T_observed=Tm[c,1])
    }
    audits[[as.character(stride)]] <- list(records=rec,statistics=statistics,h2=h2,
      S=if(draw==1L)S else NULL,Tm=if(draw==1L)Tm else NULL,Lall=if(draw==1L)Lall else NULL,
      X=if(draw==1L)X else NULL,targets=scored$targets,coordinate_scale=st)
  }
  result <- list(rows=do.call(rbind,rows),aucs=do.call(rbind,aucs),audits=audits,draw=draw,q=q)
  saveRDS(result,paste0(path,".part"));stopifnot(file.rename(paste0(path,".part"),path))
  cat(sprintf("Real analysis completed draw %d question %s\n",draw,input$qnames[q]));flush.console()
  result$rows
}
load_real_runtime(root,src)
config <- list(B=as.integer(arg("B","4999")),out=arg("out",file.path(root,"real_analysis")))
draws <- as.integer(arg("draws","10"));cores <- as.integer(arg("cores","30"));qs <- as.integer(strsplit(arg("questions","1,2,3,4"),",",fixed=TRUE)[[1]])
stopifnot(config$B>=19L,draws>=1L,draws<=10L,all(qs%in%1:4))
dir.create(file.path(config$out,"jobs"),recursive=TRUE,showWarnings=FALSE)
real <- readRDS(file.path(root,"innovations/filters.rds"));cal <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"));cal <- cal[cal$dimension==5L,]
valid <- tolower(as.character(cal$complete_all_main_stocks))=="true"
input <- list(fit=innovation_array(real$Z,real$fitting,5L),test=innovation_array(real$Z,real$dates>=as.Date("2018-01-30"),5L),
  masks=list(fit=valid[cal$split=="fit"],test=valid[cal$split=="test"]),qnames=c("within_bank","within_pharma","within_machinery","between_sectors"))
stopifnot(sum(input$masks$fit)==116L,sum(input$masks$test)==311L)
fingerprint <- list(config=config,draws=draws,questions=qs,hashes=tools::md5sum(c(file.path(root,c("run_real_analysis.R","real_analysis/protocol.json",
  "sampling_calibration.R","sampling_calibration_fullfit.R","innovation_validation.R","panel_calibration.R","native/panel_auc.dll","innovations/filters.rds")),file.path(src,c("learners.R","learners_lib.R","tests.R")))))
fp <- file.path(config$out,"fingerprint.rds");if(file.exists(fp))stopifnot(identical(fingerprint,readRDS(fp)))else saveRDS(fingerprint,fp)
writeLines(capture.output(sessionInfo()),file.path(config$out,"sessionInfo.txt"))
jobs <- unlist(lapply(seq_len(draws),function(draw)lapply(qs,function(q)list(draw=draw,q=q))),recursive=FALSE)
parallel_job <- function(job,input,config) real_job(job,input,config)
if(cores>1L) {
  cl <- parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(config$out,"workers.log"))
  parallel::clusterCall(cl,load_real_runtime,root,src)
  parallel::clusterExport(cl,c("real_job","real_seed"),envir=.GlobalEnv)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,parallel_job,input,config)),finally=parallel::stopCluster(cl))
} else invisible(lapply(jobs,real_job,input=input,config=config))
saved <- lapply(jobs,function(job)readRDS(file.path(config$out,"jobs",sprintf("draw%02d_q%d.rds",job$draw,job$q))))
x <- do.call(rbind,lapply(saved,`[[`,"rows"));x$p_holm_exploratory <- NA_real_
for(draw in seq_len(draws))for(stride in c(1L,2L,4L))for(method in unique(x$method)) {
  ix <- which(x$draw==draw & x$stride==stride & x$method==method)
  x$p_holm_exploratory[ix] <- p.adjust(x$p_exploratory[ix],"holm")
}
write.csv(x,file.path(config$out,"results.csv"),row.names=FALSE)
write.csv(x[x$primary,],file.path(config$out,"primary_results.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(saved,`[[`,"aucs")),file.path(config$out,"learner_auc.csv"),row.names=FALSE)
cat(sprintf("Saved %d real-data exploratory method results; no calibrated significance claim.\n",nrow(x)))
