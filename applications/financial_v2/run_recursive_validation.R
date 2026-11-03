args <- commandArgs(TRUE)
arg <- function(key,default) {a <- grep(paste0("^--",key,"="),args,value=TRUE);if(length(a)) sub(paste0("^--",key,"="),"",a[1L]) else default}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_runtime <- function(root,src) {
  for(pkg in c("glmnet","ranger")) {
    available <- FALSE
    for(a in 1:3) if(isTRUE(tryCatch(requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))) {available <- TRUE;break}
    if(!available) stop("Package unavailable: ",pkg)
  }
  for(f in c("learners.R","learners_lib.R","tests.R")) sys.source(file.path(src,f),envir=.GlobalEnv)
  for(f in c("panel_auc.dll","garch_filter.dll")) dyn.load(file.path(root,"native",f))
  for(f in c("panel_calibration.R","garch_innovations.R","innovation_validation.R","recursive_innovation_bootstrap.R")) sys.source(file.path(root,f),envir=.GlobalEnv)
  install_panel_knn();invisible(TRUE)
}
cat("Loading recursive model-bootstrap candidate...\n");flush.console();load_runtime(root,src)
recursive_model_smoke()
config <- list(d=5L,B=as.integer(arg("B","99")),outdir=arg("out",file.path(root,"recursive_innovation_pilot_d5")))
scenarios <- strsplit(arg("scenarios","null_iid,null_garch99,null_heterogeneous"),",",fixed=TRUE)[[1L]]
reps <- as.integer(arg("reps","20"));cores <- as.integer(arg("cores","30"));first <- as.integer(arg("first","1"))
stopifnot(all(scenarios %in% INNOVATION_SCENARIOS),config$B<5000L,reps+first-1L<=1000L)
real <- readRDS(file.path(root,"innovations/filters.rds"))
cal <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"));cal <- cal[cal$dimension==5L,]
valid <- tolower(as.character(cal$complete_all_main_stocks))=="true"
input <- list(valid=is.finite(real$X),fitting_days=real$fitting,testing_days=real$dates>=as.Date("2018-01-30"),
              masks=list(fit=valid[cal$split=="fit"],test=valid[cal$split=="test"]))
dir.create(config$outdir,recursive=TRUE,showWarnings=FALSE)
fingerprint <- list(config=config[setdiff(names(config),"B")],
  hashes=tools::md5sum(c(file.path(root,c("run_recursive_validation.R","recursive_innovation_bootstrap.R","innovation_validation.R",
    "garch_innovations.R","panel_calibration.R","native/panel_auc.dll","native/garch_filter.dll","innovations/filters.rds")),
    file.path(src,c("learners.R","learners_lib.R","tests.R")))))
path <- file.path(config$outdir,"fingerprint.rds")
if(file.exists(path)) stopifnot(identical(fingerprint,readRDS(path))) else saveRDS(fingerprint,path)
writeLines(capture.output(sessionInfo()),file.path(config$outdir,"sessionInfo.txt"))
jobs <- list()
for(case in scenarios) for(r in seq.int(first,first+reps-1L))
  if(!file.exists(file.path(config$outdir,case,sprintf("r%05d/result_B%04d.csv",r,config$B)))) jobs[[length(jobs)+1L]] <- list(scenario=case,rep=r)
write.csv(data.frame(B=config$B,replications_per_case=reps,first_replication=first,cores=cores,
                     scenarios=paste(scenarios,collapse=",")),file.path(config$outdir,sprintf("run_request_B%04d.csv",config$B)),row.names=FALSE)
cat(sprintf("Recursive pilot: %d cases x %d outer reps; B=%d inner null datasets, full filters and classifiers; pending=%d cores=%d.\n",length(scenarios),reps,config$B,length(jobs),cores));flush.console()
started <- proc.time()[["elapsed"]]
if(length(jobs)) {
  cl <- parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(config$outdir,"workers.log"))
  parallel::clusterCall(cl,load_runtime,root,src)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,function(job,config,input) recursive_job(job,config,input),config,input)),
           finally=parallel::stopCluster(cl))
}
cat(sprintf("Recursive pilot finished in %.1f minutes.\n",(proc.time()[["elapsed"]]-started)/60))
