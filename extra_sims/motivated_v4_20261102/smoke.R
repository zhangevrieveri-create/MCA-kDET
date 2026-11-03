# Isolated exploratory scenarios: never edits the frozen generators or tests.
SRC_DIR <- file.path(getwd(),"R"); PROFILE <- "quick"
for(p in c("glmnet","ranger","gTestsMulti"))stopifnot(requireNamespace(p,quietly=TRUE))
for(f in c("models.R","learners.R","learners_lib.R","tests.R","sim.R","garch_ext.R",
           "competitors2.R","rotation.R","learners_lib2.R","models_s.R","experiments.R"))source(file.path(SRC_DIR,f))
folder <- "motivated_v4_20261102"
source(file.path(folder,"models.R"));source("R/oracle_s.R")
calibration <- readRDS(file.path(folder,"calibration.rds"))
out <- file.path(folder,"smoke_B19")
args <- commandArgs(TRUE);action <- args[1L]

xscr_diagnostic <- function(sc,m) {
  d <- ncol(sc$Zfit);pairs <- which(upper.tri(diag(d)),arr.ind=TRUE)
  F <- xscr_features(sc$Zfit,pairs);T <- xscr_features(sc$Zte,pairs)
  label <- c(paste0("x",1:d),paste0("x",1:d,"^2"),paste0("x",pairs[,1],"*x",pairs[,2]))
  ans <- list()
  for(direction in 1:2) {
    A <- if(direction==1L)F else T;C <- if(direction==1L)T else F
    g <- if(direction==1L)sc$gfit else sc$L;h <- if(direction==1L)sc$L else sc$gfit
    mu <- colMeans(A);ss <- colsd(A);ss[!is.finite(ss)|ss<1e-12] <- 1
    Az <- sweep(sweep(A,2L,mu),2L,ss,"/");Cz <- sweep(sweep(C,2L,mu),2L,ss,"/")
    for(i in 1:sc$k) {
      pos <- g==i;n1 <- sum(pos);n0 <- length(g)-n1;s1 <- colSums(Az[pos,,drop=FALSE])
      t <- (s1/n1-(colSums(Az)-s1)/n0)/sqrt(1/n1+1/n0)
      sel <- top_m(abs(t),m);score <- as.numeric(Cz[,sel,drop=FALSE]%*%t[sel])
      auc <- auc_vec(score,h==i)
      ans[[length(ans)+1L]] <- data.frame(m=m,direction=direction,population=i,
        selected=paste(label[sel],collapse=";"),auc=auc,
        T=(auc-.5)/sqrt((length(h)+1)/(12*sum(h==i)*sum(h!=i))))
    }
  }
  ans <- do.call(rbind,ans)
  actual <- xscr_T(sc$Zfit,sc$gfit,sc$Zte,matrix(sc$L,ncol=1L),sc$k,m)[,1L]
  rebuilt <- vapply(1:sc$k,function(i)sum(ans$T[ans$population==i])/sqrt(2),numeric(1))
  stopifnot(max(abs(actual-rebuilt))<1e-10)
  ans
}

one <- function(cell,r) {
  seed <- as.integer(cell$seed0+r);set.seed(seed);t0 <- proc.time()[["elapsed"]]
  data <- mv4_generate(cell$model,cell$target,rep(300L,cell$k),cell$d,calibration)
  sc <- fit_scores(split_data(data,.5),"ovr",lib=parse_lib(cell$lib))
  ct <- classifier_tests(sc,19L,parse_libsets(cell$libsets,parse_lib(cell$lib)))
  comp <- competitor_tests(data,19L,strsplit(cell$groups,",",fixed=TRUE)[[1L]])
  aucs <- do.call(rbind,lapply(names(sc$Sovr),function(a)data.frame(learner=a,population=1:cell$k,
       auc=vapply(1:cell$k,function(i)auc_vec(sc$Sovr[[a]][,i],sc$L==i),numeric(1)))))
  xscr <- rbind(xscr_diagnostic(sc,1L),xscr_diagnostic(sc,3L))
  list(cell=cell$cell,rep=r,seed=seed,B=19L,p=c(ct$p,comp$p),diag=ct$diag,aucs=aucs,xscr=xscr,
       seconds=proc.time()[["elapsed"]]-t0)
}

if(action=="prepare") {
  stopifnot(!dir.exists(out));dir.create(out);dir.create(file.path(out,"rows"));dir.create(file.path(out,"logs"))
  specs <- list()
  for(d in c(20L,100L))for(model in c(MV4_FAMILIES,"N4_mixed"))for(target in if(d==20L)MV4_TARGETS else .70)
    specs[[length(specs)+1L]] <- list(model=model,target=target,d=d,k=if(model=="N4_mixed")4L else 3L)
  for(d in c(20L,100L))for(k in 3:4)specs[[length(specs)+1L]] <- list(model=if(k==3L)"N1_bimodal" else "N4_mixed",target=.5,d=d,k=k)
  grid <- do.call(rbind,lapply(seq_along(specs),function(j) {
    z <- specs[[j]];g <- cell_row("motiv4smoke",z$model,z$target-.5,k=z$k,n=300L,d=z$d,B=19L,reps=3L,
      groups=if(z$k==4L)"ovr,dist,schott,comp2,comp2ovr,songchen" else LP3_GROUPS,
      lib=LIBPILOT3_LIB,libsets=LIBPILOT3_SETS,seed_scheme="v2")
    # Dedicated range strictly between all legacy v1 and v2 seeds.
    g$seed0 <- 150000000L+j*10000L;g$target <- z$target;g
  }))
  stopifnot(nrow(grid)==20L,!anyDuplicated(grid$cell),!anyDuplicated(grid$seed0),
    min(grid$seed0)>PAR$seed+99999*1000,max(grid$seed0)+500<200000000)
  saveRDS(grid,file.path(out,"grid.rds"));write.csv(grid,file.path(out,"grid.csv"),row.names=FALSE)
  files <- c(list.files("R",pattern="[.]R$",full.names=TRUE),list.files(folder,pattern="[.]R$",full.names=TRUE),file.path(folder,"calibration.rds"),file.path(folder,"design.json"))
  write.csv(data.frame(file=files,md5=as.character(tools::md5sum(files))),file.path(out,"source_md5.csv"),row.names=FALSE)
  # Independent calibration validation uses new random numbers, not calibration samples.
  noise <- mv4_noise(100000L,20261104L)
  checks <- do.call(rbind,lapply(MV4_FAMILIES,function(f)do.call(rbind,lapply(MV4_TARGETS,function(a) {
    th <- calibration[[f]][[sprintf("%.2f",a)]];actual <- mv4_oracle(f,th,noise)
    stopifnot(abs(actual-a)<.006)
    data.frame(family=f,target=a,theta=th,independent_auc=actual)
  }))))
  write.csv(checks,file.path(folder,"calibration_check.csv"),row.names=FALSE)
  cat("PREPARED 20 CELLS; INDEPENDENT AUC CHECKS PASS\n")
} else if(action=="cell") {
  grid <- readRDS(file.path(out,"grid.rds"));cell <- as.list(grid[as.integer(args[2L]),])
  manifest <- read.csv(file.path(out,"source_md5_revision1.csv"));stopifnot(identical(as.character(tools::md5sum(manifest$file)),manifest$md5))
  for(r in 1:3) {
    path <- file.path(out,"rows",sprintf("%s__r%05d.rds",cell$cell,r))
    if(file.exists(path))next
    cat("START",cell$cell,r,"\n");flush.console();result <- NULL
    for(attempt in 1:3) {
      result <- tryCatch(one(cell,r),error=function(e)conditionMessage(e))
      if(is.list(result))break
      cat("ERROR",attempt,result,"\n");flush.console()
    }
    if(!is.list(result))stop(result)
    # The old comparator API also returns kmd_eta (an effect statistic) in $p.
    # Preserve it in the raw record, but never validate or report it as a p-value.
    pv <- unlist(result$p[setdiff(names(result$p),"kmd_eta")])
    stopifnot(all(is.finite(pv)),all(pv>=0&pv<=1))
    if(as.integer(args[2L])==1L&&r==1L) {
      again <- one(cell,r)
      stopifnot(identical(result$p,again$p),identical(result$aucs,again$aucs),identical(result$xscr,again$xscr))
      cat("SAME-SEED REPEAT EXACTLY IDENTICAL\n")
    }
    saveRDS(result,paste0(path,".tmp"));stopifnot(!file.exists(path),file.rename(paste0(path,".tmp"),path))
    cat("SAVED",r,"seconds",result$seconds,"\n");flush.console()
  }
} else if(action=="report") {
  grid <- readRDS(file.path(out,"grid.rds"));files <- list.files(file.path(out,"rows"),pattern="[.]rds$",full.names=TRUE)
  rows <- lapply(files,readRDS);stopifnot(length(rows)==60L)
  ids <- vapply(rows,function(x)paste(x$cell,x$rep),character(1));stopifnot(!anyDuplicated(ids))
  auc <- do.call(rbind,lapply(rows,function(x)cbind(cell=x$cell,rep=x$rep,x$aucs)))
  xscr <- do.call(rbind,lapply(rows,function(x)cbind(cell=x$cell,rep=x$rep,x$xscr)))
  p <- do.call(rbind,lapply(rows,function(x) {
    keep <- setdiff(names(x$p),"kmd_eta")
    data.frame(cell=x$cell,rep=x$rep,method=keep,p=unlist(x$p[keep]))
  }))
  write.csv(auc,file.path(out,"aucs.csv"),row.names=FALSE);write.csv(xscr,file.path(out,"xscr_selected.csv"),row.names=FALSE)
  write.csv(p,file.path(out,"pvalues.csv"),row.names=FALSE)
  summary <- aggregate(auc~cell+learner+population,auc,mean)
  write.csv(summary,file.path(out,"auc_summary.csv"),row.names=FALSE)
  manifest <- read.csv(file.path(out,"source_md5_revision1.csv"));stopifnot(identical(as.character(tools::md5sum(manifest$file)),manifest$md5))
  stopifnot(all(table(vapply(rows,`[[`,character(1),"cell"))==3L))
  writeLines(c("60/60 replications, 20/20 cells, B=19; no missing or duplicate replications.",
    "Original statistics, libraries and comparator sources unchanged.",
    "All observed xscr diagnostics match the production statistic; one complete same-seed repetition matched exactly.",
    "Three repeats do not estimate power. Do not select models from these outcomes."),file.path(out,"verification.txt"))
  cat("SMOKE COMPLETE AND VERIFIED\n")
} else stop("Unknown action")
