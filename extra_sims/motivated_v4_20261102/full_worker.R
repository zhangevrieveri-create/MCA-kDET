source("motivated_v4_20261102/full_common.R")
action <- commandArgs(TRUE)[1L]
out <- motiv_arg("out",file.path(MOTIV_DIR,"full_B200"))
if(action=="prepare") {
  stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
  for(d in c("rows","worker_logs","errors","report_tables"))dir.create(file.path(out,d))
  grid <- readRDS(file.path(MOTIV_DIR,"study_grid.rds"))
  stopifnot(nrow(grid)>0,!anyDuplicated(grid$cell),!anyDuplicated(grid$seed0),
            all(grid$B==200L),all(grid$reps==ifelse(grid$target<=.5,500L,100L)))
  smoke <- identical(motiv_arg("smoke","false"),"true")
  if(smoke){grid$B <- 19L;grid$reps <- 3L;grid$seed0 <- grid$seed0+2000L}
  saveRDS(grid,file.path(out,"grid.rds"));write.csv(grid,file.path(out,"grid.csv"),row.names=FALSE)
  files <- unique(c(list.files("R",pattern="[.]R$",full.names=TRUE),
    list.files(MOTIV_DIR,pattern="[.](R|ps1|json|rds)$",full.names=TRUE),
    file.path(out,"grid.rds"),file.path(out,"grid.csv")))
  stopifnot(all(file.exists(files)))
  write.csv(data.frame(file=files,md5=as.character(tools::md5sum(files))),file.path(out,"source_md5.csv"),row.names=FALSE)
  writeLines(c("User authorized 2026-11-02: all fixed scenarios and methods; no outcome-based filtering.",
    "Per-rep attempt bound 3; process attempt bound 3 (same seeds, at most 9 calculations after repeated process failure).",
    "Level deviations are retained and reported; no threshold, seed, sample or method adjustment.",
    paste("Smoke only:",smoke)),file.path(out,"run_authorization.txt"))
  cat("PREPARED",nrow(grid),"CELLS;",sum(grid$reps),"REPLICATIONS; B",unique(grid$B),"\n")
} else if(action=="cell") {
  motiv_manifest_check(out)
  grid <- readRDS(file.path(out,"grid.rds"));idx <- as.integer(motiv_arg("cell"))
  first <- as.integer(motiv_arg("start"));last <- as.integer(motiv_arg("end"))
  stopifnot(idx>=1L,idx<=nrow(grid),first>=1L,last<=grid$reps[idx],first<=last)
  cell <- as.list(grid[idx,]);
  for(r in first:last) {
    path <- file.path(out,"rows",sprintf("%s__r%05d.rds",cell$cell,r))
    if(file.exists(path)) {
      old <- readRDS(path);stopifnot(old$seed==cell$seed0+r,old$B==cell$B,
        old$cell==cell$cell,old$rep==r,old$err=="",old$lib==cell$lib,old$libsets==cell$libsets)
      next
    }
    cat("START",format(Sys.time(),"%H:%M:%S"),cell$cell,"rep",r,"\n");flush.console()
    result <- NULL
    for(a in 1:3) {
      result <- tryCatch(motiv_one(cell,r),error=function(e)conditionMessage(e))
      if(is.list(result))break
      msg <- paste(format(Sys.time()),"rep",r,"attempt",a,result)
      cat("ERROR",msg,"\n");flush.console()
      cat(msg,"\n",file=file.path(out,"errors",sprintf("%s__r%05d.log",cell$cell,r)),append=TRUE)
    }
    if(!is.list(result))stop("Same-seed repetition failed three times: ",result)
    result$retries <- a-1L
    if(identical(motiv_arg("repeat_check","false"),"true") && r==first) {
      again <- motiv_one(cell,r)
      stopifnot(identical(result$p,again$p),identical(result$aucs,again$aucs),identical(result$xscr,again$xscr))
      result$same_seed_exact <- TRUE
    }
    motiv_manifest_check(out)
    saveRDS(result,paste0(path,".tmp"))
    stopifnot(!file.exists(path),file.rename(paste0(path,".tmp"),path))
    cat("SAVED",format(Sys.time(),"%H:%M:%S"),r,"seconds",result$seconds,"\n");flush.console()
  }
} else stop("Unknown action")
