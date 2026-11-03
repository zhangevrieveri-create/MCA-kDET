# Independent held-out record and pairwise-comparison audit; no model refitting.
args <- commandArgs(TRUE)
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
folder <- if(length(args)) args[1] else file.path(root,"sampling_fullfit_confirmation_d5/estimated")
fingerprint <- file.path(dirname(folder),"fingerprint.rds")
if(file.exists(fingerprint)) {
  expected <- readRDS(fingerprint)$hashes
  stopifnot(identical(expected,tools::md5sum(names(expected))))
}
paths <- list.files(file.path(folder,"audits"),pattern="_r[0-9]+[.]rds$",full.names=TRUE)
stopifnot(length(paths)>0L)
calendar <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"))
calendar <- calendar[calendar$dimension==5L,]
fit_valid <- tolower(as.character(calendar$complete_all_main_stocks[calendar$split=="fit"]))=="true"
test_valid <- tolower(as.character(calendar$complete_all_main_stocks[calendar$split=="test"]))=="true"
qnames <- c("within_bank","within_pharma","within_machinery","between_sectors")
rows <- list()
for(path in paths) {
  audit <- readRDS(path)
  candidates <- list.files(file.path(folder,"chunks"),pattern=paste0("^",audit$scenario,"__r[0-9]+-[0-9]+[.]csv$"),full.names=TRUE)
  bounds <- lapply(basename(candidates),function(s) as.integer(regmatches(s,gregexpr("[0-9]{5}",s))[[1]]))
  index <- which(vapply(bounds,function(b)length(b)==2L && b[1]<=audit$rep && b[2]>=audit$rep,logical(1)))
  # Audits can precede the atomic write of their containing ten-repetition chunk.
  if(!length(index)) next
  stopifnot(length(index)==1L)
  data <- read.csv(candidates[index])
  for(key in names(audit$qa)) {
    ids <- as.integer(strsplit(key,"_",fixed=TRUE)[[1]]);q <- ids[1];stride <- c(1L,2L,4L)[ids[2]]
    a <- audit$qa[[key]];rec <- a$testing_records;fit <- a$fitting_records
    k <- if(q<=3L)5L else 3L
    stopifnot(!anyDuplicated(rec$time),all(diff(sort(rec$time))>=stride),all(test_valid[rec$time]),
              length(unique((rec$time-1L)%%stride))==1L,
              all(table(rec$population)==nrow(rec)/k),identical(as.integer(a$Lall[,1]),as.integer(rec$population)),
              nrow(fit)==sum(fit_valid)*if(q<=3L)5L else 15L,
              !anyDuplicated(paste(fit$time,fit$stock)),all(fit_valid[fit$time]),
              nrow(a$S)==nrow(rec),ncol(a$S)==k*4L,all(is.finite(a$S)),all(is.finite(a$Tm)))
    if(q<=3L)stopifnot(all(rec$stock==(q-1L)*5L+rec$population)) else
      stopifnot(all((rec$stock-1L)%/%5L+1L==rec$population))
    stopifnot(all(vapply(seq_len(ncol(a$Lall)),function(b)
      all(tabulate(a$Lall[,b],k)==nrow(rec)/k),logical(1))))
    draws <- unique(c(1L,2L,17L,ncol(a$Lall)));error <- 0
    for(b in draws)for(c in seq_len(ncol(a$S))) {
      pos <- a$S[a$Lall[,b]==a$targets[c],c];neg <- a$S[a$Lall[,b]!=a$targets[c],c]
      comparisons <- outer(pos,neg,"-")
      auc <- mean(comparisons>0)+.5*mean(comparisons==0)
      sigma <- sqrt((length(pos)+length(neg)+1)/(12*length(pos)*length(neg)))
      error <- max(error,abs((auc-.5)/sigma-a$Tm[c,b]))
    }
    stopifnot(error<3e-13)
    statistic <- apply(a$Tm,2,max)
    p <- (1+sum(statistic[-1]>=statistic[1]-1e-10))/length(statistic)
    observed <- data[data$rep==audit$rep & data$block_length==stride & data$question==qnames[q],]
    stopifnot(nrow(observed)==1L,abs(observed$p-p)<1e-14,
              abs(observed$statistic-statistic[1])<1e-12,
              observed$testing_rows==nrow(rec),observed$fitting_rows==nrow(fit))
    rows[[length(rows)+1L]] <- data.frame(scenario=audit$scenario,rep=audit$rep,question=qnames[q],stride=stride,
      population_test_n=nrow(rec)/k,pairwise_draws_checked=length(draws),max_statistic_error=error,p=p,
      records_and_labels_valid=TRUE)
  }
}
stopifnot(length(rows)>0L)
result <- do.call(rbind,rows)
write.csv(result,file.path(folder,"independent_sampling_audit.csv"),row.names=FALSE)
cat(sprintf("Independently audited %d saved question/stride cells from %d scenarios; pairwise AUC error %.3g. This is an implementation audit, not a size validation.\n",
  nrow(result),length(unique(result$scenario)),max(result$max_statistic_error)))
