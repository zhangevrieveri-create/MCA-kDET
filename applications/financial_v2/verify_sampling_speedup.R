root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
args <- commandArgs(TRUE)
slow <- if(length(args)>=1L) args[1L] else file.path(root,"../../tests/out/financial_sampling_smoke")
fast <- if(length(args)>=2L) args[2L] else file.path(root,"../../tests/out/financial_sampling_fast_smoke")
phases <- if(length(args)>=3L) strsplit(args[3L],",",fixed=TRUE)[[1L]] else c("oracle","estimated")
rows <- list()
for(phase in phases)for(case in c("null_iid","null_garch99","null_equal_mixtures")) {
  path <- file.path(phase,"audits",paste0(case,"_r00001.rds"))
  a <- readRDS(file.path(slow,path));b <- readRDS(file.path(fast,path))
  stopifnot(identical(names(a$qa),names(b$qa)))
  for(key in names(a$qa)) {
    x <- a$qa[[key]];y <- b$qa[[key]]
    stopifnot(identical(x$fitting_records,y$fitting_records),identical(x$testing_records,y$testing_records),
              identical(x$Lall,y$Lall),identical(x$targets,y$targets))
    error <- max(abs(x$S-y$S));stopifnot(error<1e-12)
    ranks_same <- identical(apply(x$S,2L,rank,ties.method="average"),apply(y$S,2L,rank,ties.method="average"))
    stopifnot(ranks_same,identical(x$Tm,y$Tm))
    rows[[length(rows)+1L]] <- data.frame(phase=phase,scenario=case,cell=key,max_score_error=error,
                                        all_score_ranks_identical=ranks_same,all_permutation_statistics_identical=TRUE)
  }
  name <- sprintf("%s__r00001-00001.csv",case)
  x <- read.csv(file.path(slow,phase,"chunks",name));y <- read.csv(file.path(fast,phase,"chunks",name))
  cols <- setdiff(names(x),"sparse_shortcut_max_score_error")
  stopifnot(identical(x[,cols],y[,cols]))
}
tab <- do.call(rbind,rows)
write.csv(tab,file.path(fast,"speedup_verification.csv"),row.names=FALSE)
cat(sprintf("Prediction batching verified: %d question/stride/phase cells, max score error %.3g; all ranks, permutation statistics and reported results identical.\n",nrow(tab),max(tab$max_score_error)))
