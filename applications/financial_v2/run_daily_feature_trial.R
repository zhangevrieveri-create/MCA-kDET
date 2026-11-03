# Application wrapper only. Frozen learners/statistics/label permutation functions are sourced unchanged.
args <- commandArgs(TRUE)
arg <- function(k,d){x<-grep(paste0('^--',k,'='),args,value=TRUE);if(length(x))sub(paste0('^--',k,'='),'',x[1])else d}
root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),winslash='/'))
src <- normalizePath(file.path(root,'../../R'),winslash='/')
load_daily_runtime <- function(root,src){
  stopifnot(requireNamespace('glmnet',quietly=TRUE),requireNamespace('ranger',quietly=TRUE))
  for(f in c('learners.R','learners_lib.R','tests.R'))sys.source(file.path(src,f),envir=.GlobalEnv)
  dyn.load(file.path(root,'native/panel_auc.dll'))
  sys.source(file.path(root,'panel_calibration.R'),envir=.GlobalEnv);install_panel_knn()
  invisible(TRUE)
}
daily_job <- function(job,input,config){
  sector <- c('银行','医药','机械')[ceiling(job/2)];representation<-if(job%%2==1)'full5'else'price4'
  path<-file.path(config$out,'jobs',sprintf('sector%d_%s.rds',ceiling(job/2),representation))
  if(file.exists(path))return(invisible(TRUE))
  features<-c('overnight_log_return','intraday_log_return','log_high_low_range','close_position','turnover_today')
  if(representation=='price4')features<-features[1:4]
  stocks<-input$panel$order_book_id[input$panel$first_industry_name==sector];stopifnot(length(stocks)==5L)
  records<-list();sp<-lapply(seq_along(stocks),function(g){
    x<-input$data[input$data$stock==stocks[g],];f<-x[x$split=='fit',];t<-x[x$split=='test',]
    if(config$smoke){f<-head(f,60L);t<-head(t,15L)}
    records[[g]]<<-data.frame(stock=stocks[g],population=g,date=t$date)
    list(fit=as.matrix(f[,features,drop=FALSE]),test=as.matrix(t[,features,drop=FALSE]))
  })
  for(s in sp)stopifnot(nrow(s$fit)>=30L,nrow(s$test)>=12L,all(is.finite(s$fit)),all(is.finite(s$test)))
  warnings_seen<-character();set.seed(202610280L+job)
  start<-proc.time()[['elapsed']]
  sc<-withCallingHandlers(fit_scores(sp,'ovr',lib=c('l1','quad','knn','rf')),
       warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
  S<-do.call(cbind,sc$Sovr);targets<-rep(1:5,4L)
  st<-standardize_fit(do.call(rbind,lapply(sp,`[[`,'fit')))
  X<-apply_std(do.call(rbind,lapply(sp,`[[`,'test')),st)
  # Same permutation seed/observed samples for full5 and price4 within a sector.
  set.seed(202610280L+100L+ceiling(job/2));Lall<-perm_labels(sc$L,config$B)
  Tm<-ovr_T(S,targets,Lall);n<-vapply(sp,function(s)nrow(s$test),integer(1));nb<-ncol(Lall)
  G<-group_indicator(Lall,5L);D2<-sqdist(X);D<-sqrt(D2);h2<-median(D2[upper.tri(D2)]);stopifnot(h2>0)
  K<-exp(-D2/(2*h2));E<-mst_edges(D)
  statistics<-list('MCA-kDET'=apply(Tm,2,max),DISCO=disco_stat(D,G,n,5L,nb),
     MMD=mmd_stat(K,G,n,5L,nb),MST=colSums(Lall[E[,1],,drop=FALSE]==Lall[E[,2],,drop=FALSE]))
  elapsed<-proc.time()[['elapsed']]-start
  rows<-do.call(rbind,lapply(names(statistics),function(method){v<-statistics[[method]]
    data.frame(sector=sector,representation=representation,method=method,d=length(features),B=config$B,
      iid_reference_p_unvalidated=pval_perm(v),statistic=v[1],test_n=paste(n,collapse='/'),
      fitting_n=paste(vapply(sp,function(s)nrow(s$fit),integer(1)),collapse='/'),seconds=elapsed,
      status='IID_REFERENCE_ONLY_NO_DEPENDENT_DATA_SIGNIFICANCE_CLAIM')}))
  auc<-do.call(rbind,lapply(seq_len(ncol(S)),function(j)data.frame(stock=stocks[targets[j]],population=targets[j],
    learner=rep(c('l1','quad','knn','rf'),each=5L)[j],auc=auc_vec(S[,j],sc$L==targets[j]),T=Tm[j,1])))
  saveRDS(list(rows=rows,auc=auc,records=do.call(rbind,records),statistics=statistics,S=S,Tm=Tm,Lall=Lall,
    X=X,targets=targets,h2=h2,coordinate_scale=st,features=features,stocks=stocks,warnings=warnings_seen),path)
  cat(sprintf('Completed sector %d %s, %d pooled test observations, %.1f seconds\n',ceiling(job/2),representation,nrow(X),elapsed));flush.console()
  invisible(TRUE)
}
smoke<-identical(arg('smoke','false'),'true');B<-as.integer(arg('B',if(smoke)'39'else'999'))
out<-arg('out',file.path(root,'daily_features_v1',if(smoke)'smoke'else'analysis'))
cores<-as.integer(arg('cores',if(smoke)'1'else'6'));jobs<-as.integer(strsplit(arg('jobs',if(smoke)'1'else'1,2,3,4,5,6'),',',fixed=TRUE)[[1]])
config<-list(out=out,B=B,smoke=smoke);dir.create(file.path(out,'jobs'),recursive=TRUE,showWarnings=FALSE)
input<-list(data=read.csv(gzfile(file.path(root,'daily_features_v1/complete_daily_features.csv.gz')),fileEncoding='UTF-8'),
            panel=read.csv(file.path(root,'design_two_tests/main_panel.csv'),fileEncoding='UTF-8'))
paths<-c(file.path(root,c('run_daily_feature_trial.R','daily_features_v1/protocol.json','daily_features_v1/complete_daily_features.csv.gz','panel_calibration.R','native/panel_auc.dll')),
         file.path(src,c('models.R','calib_params.R','calib_params_v2.R','tests.R','learners.R','learners_lib.R')))
fingerprint<-list(paths=paths,hashes=tools::md5sum(paths),config=config)
fp<-file.path(out,'fingerprint.rds');if(file.exists(fp))stopifnot(identical(fingerprint,readRDS(fp)))else saveRDS(fingerprint,fp)
if(cores>1L){
  cl<-parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(out,'workers.log'))
  parallel::clusterCall(cl,load_daily_runtime,root,src);parallel::clusterExport(cl,'daily_job',envir=.GlobalEnv)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,function(j,input,config)daily_job(j,input,config),input,config)),finally=parallel::stopCluster(cl))
}else{load_daily_runtime(root,src);invisible(lapply(jobs,daily_job,input=input,config=config))}
cat('Daily-feature trial completed. All permutation numbers are IID reference diagnostics, not validated market p-values.\n')
