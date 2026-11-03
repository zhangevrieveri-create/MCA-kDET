# End-to-end validation for the explicitly changed innovation-distribution target.
args <- commandArgs(TRUE)
arg <- function(key,default) {
  a <- grep(paste0("^--",key,"="),args,value=TRUE)
  if(length(a)) sub(paste0("^--",key,"="),"",a[1L]) else default
}
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
src <- normalizePath(file.path(root,"../../R"),winslash="/")
load_runtime <- function(root,src) {
  for(pkg in c("glmnet","ranger")) {
    okay <- FALSE
    for(a in 1:3) if(isTRUE(tryCatch(requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))) {okay <- TRUE;break}
    if(!okay) stop("Package unavailable: ",pkg)
  }
  for(f in c("learners.R","learners_lib.R","tests.R")) sys.source(file.path(src,f),envir=.GlobalEnv)
  for(f in c("panel_auc.dll","garch_filter.dll")) dyn.load(file.path(root,"native",f))
  for(f in c("panel_calibration.R","garch_innovations.R","innovation_validation.R")) sys.source(file.path(root,f),envir=.GlobalEnv)
  install_panel_knn()
  invisible(TRUE)
}
cat("Loading frozen learner runtime and application helpers...\n");flush.console()
load_runtime(root,src)
cat("Runtime ready; reading fixed calendar and missingness masks.\n");flush.console()
config <- list(d=5L,B=as.integer(arg("B","999")),lengths=c(1L,4L,8L),seed0=340000000L,
               outdir=arg("out",file.path(root,"innovation_validation_pilot_d5")))
scenarios <- strsplit(arg("scenarios",paste(INNOVATION_SCENARIOS[1:6],collapse=",")),",",fixed=TRUE)[[1L]]
stopifnot(all(scenarios %in% INNOVATION_SCENARIOS))
reps <- as.integer(arg("reps","100"));cores <- as.integer(arg("cores","30"));first_rep <- as.integer(arg("first","1"))
real <- readRDS(file.path(root,"innovations/filters.rds"))
calendar <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"));calendar <- calendar[calendar$dimension==config$d,]
valid <- tolower(as.character(calendar$complete_all_main_stocks))=="true"
input <- list(valid=is.finite(real$X),fitting_days=real$fitting,testing_days=real$dates>=as.Date("2018-01-30"),
              masks=list(fit=valid[calendar$split=="fit"],test=valid[calendar$split=="test"]))
for(phase in c("oracle","estimated","filters")) for(name in c("chunks","audits"))
  dir.create(file.path(config$outdir,phase,name),recursive=TRUE,showWarnings=FALSE)
fingerprint <- list(config=config,registry=INNOVATION_SCENARIOS,
  hashes=tools::md5sum(c(file.path(root,c("run_innovation_validation.R","innovation_validation.R","garch_innovations.R","panel_calibration.R",
                                        "native/garch_filter.cpp","native/garch_filter.dll","native/panel_auc.cpp","native/panel_auc.dll",
                                        "design_two_tests/block_calendar.csv","innovations/filters.rds")),
                         file.path(src,c("learners.R","learners_lib.R","tests.R")))))
path <- file.path(config$outdir,"fingerprint.rds")
if(file.exists(path)) stopifnot(identical(fingerprint,readRDS(path))) else saveRDS(fingerprint,path)
writeLines(capture.output(sessionInfo()),file.path(config$outdir,"sessionInfo.txt"))
jobs <- list()
for(case in scenarios) for(first in seq.int(first_rep,first_rep+reps-1L,by=10L)) {
  job <- list(scenario=case,start=first,end=min(first+9L,first_rep+reps-1L))
  filename <- sprintf("%s__r%05d-%05d.csv",case,job$start,job$end)
  if(!all(file.exists(file.path(config$outdir,c("oracle","estimated","filters"),"chunks",filename)))) jobs[[length(jobs)+1L]] <- job
}
cat(sprintf("Innovation validation: %d scenarios x %d reps, oracle + fitted filters, B=%d, %d pending chunks, %d cores.\n",length(scenarios),reps,config$B,length(jobs),cores));flush.console()
started <- proc.time()[["elapsed"]]
if(length(jobs)) {
  cl <- parallel::makePSOCKcluster(min(cores,length(jobs)),outfile=file.path(config$outdir,"workers.log"))
  parallel::clusterCall(cl,load_runtime,root,src)
  tryCatch(invisible(parallel::parLapplyLB(cl,jobs,function(job,config,input) innovation_job(job,config,input),config,input)),
           finally=parallel::stopCluster(cl))
}
cat(sprintf("Innovation validation finished in %.1f minutes.\n",(proc.time()[["elapsed"]]-started)/60))
