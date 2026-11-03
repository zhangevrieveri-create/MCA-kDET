args <- commandArgs(TRUE)
arg <- function(key,default){a <- grep(paste0("^--",key,"="),args,value=TRUE);if(length(a))sub(paste0("^--",key,"="),"",a[1])else default}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_runtime <- function(root,src) {
  for(pkg in c("glmnet","ranger")) {
    okay <- FALSE
    for(a in 1:3) if(isTRUE(tryCatch(requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))) {okay <- TRUE;break}
    if(!okay)stop("Package unavailable: ",pkg)
  }
  for(f in c("learners.R","learners_lib.R","tests.R"))sys.source(file.path(src,f),envir=.GlobalEnv)
  for(f in c("panel_auc.dll","garch_filter.dll"))dyn.load(file.path(root,"native",f))
  for(f in c("panel_calibration.R","garch_innovations.R","innovation_validation.R","sampling_calibration.R"))sys.source(file.path(root,f),envir=.GlobalEnv)
  install_panel_knn();invisible(TRUE)
}
cat("Loading one-stock-per-time sampling validation...\n");flush.console();load_runtime(root,src)
config <- list(d=5L,B=as.integer(arg("B","999")),strides=c(1L,2L,4L),
  outdir=arg("out",file.path(root,"sampling_validation_pilot_d5")))
scenarios <- strsplit(arg("scenarios",paste(c(SAMPLING_SCENARIOS[1:6],"null_equal_mixtures"),collapse=",")),",",fixed=TRUE)[[1L]]
stopifnot(all(scenarios %in% SAMPLING_SCENARIOS))
reps <- as.integer(arg("reps","100"));first <- as.integer(arg("first","1"));cores <- as.integer(arg("cores","30"))
real <- readRDS(file.path(root,"innovations/filters.rds"))
cal <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"));cal <- cal[cal$dimension==5L,]
valid <- tolower(as.character(cal$complete_all_main_stocks))=="true"
input <- list(valid=is.finite(real$X),fitting_days=real$fitting,testing_days=real$dates>=as.Date("2018-01-30"),
  masks=list(fit=valid[cal$split=="fit"],test=valid[cal$split=="test"]))
for(phase in c("oracle","estimated"))for(name in c("chunks","audits"))dir.create(file.path(config$outdir,phase,name),recursive=TRUE,showWarnings=FALSE)
fingerprint <- list(config=config,registry=SAMPLING_SCENARIOS,hashes=tools::md5sum(c(file.path(root,
  c("run_sampling_validation.R","sampling_calibration.R","innovation_validation.R","garch_innovations.R","panel_calibration.R",
    "native/panel_auc.dll","native/garch_filter.dll","innovations/filters.rds")),file.path(src,c("learners.R","learners_lib.R","tests.R")))))
path <- file.path(config$outdir,"fingerprint.rds")
if(file.exists(path))stopifnot(identical(fingerprint,readRDS(path)))else saveRDS(fingerprint,path)
writeLines(capture.output(sessionInfo()),file.path(config$outdir,"sessionInfo.txt"))
jobs <- list()
for(case in scenarios)for(start in seq.int(first,first+reps-1L,by=10L)) {
  job <- list(scenario=case,start=start,end=min(start+9L,first+reps-1L))
  filename <- sprintf("%s__r%05d-%05d.csv",case,start,job$end)
  if(!all(file.exists(file.path(config$outdir,c("oracle","estimated"),"chunks",filename))))jobs[[length(jobs)+1L]] <- job
}
cat(sprintf("Sampling pilot: %d cases x %d reps, B=%d, strides=1/2/4 full-calendar blocks, %d pending chunks, cores=%d.\n",length(scenarios),reps,config$B,length(jobs),cores));flush.console()
started <- proc.time()[["elapsed"]]
if(length(jobs)) {
  cl <- parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(config$outdir,"workers.log"))
  parallel::clusterCall(cl,load_runtime,root,src)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,function(job,config,input)sampling_job(job,config,input),config,input)),finally=parallel::stopCluster(cl))
}
cat(sprintf("Sampling pilot completed in %.1f minutes.\n",(proc.time()[["elapsed"]]-started)/60))
