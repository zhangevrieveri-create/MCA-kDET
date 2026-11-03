# Reproduce existing A-null repetitions and inspect xscr3 permutation ranks.
# No new calibration, no parameter changes, no writes to original results.
args <- commandArgs(TRUE); stopifnot(length(args)==3L)
first <- as.integer(args[1]); last <- as.integer(args[2]); out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for(p in c("glmnet","ranger")) stopifnot(requireNamespace(p,quietly=TRUE))
for(f in c("models.R","learners.R","learners_lib.R","tests.R","sim.R","competitors2.R",
           "rotation.R","learners_lib2.R","models_s.R","experiments.R")) source(file.path(SRC_DIR,f))
original <- "results_v3_pilot_20261101_cpu16"
grid <- readRDS(file.path(original,"formal_grid.rds"))
cell <- as.list(grid[grid$model=="A" & grid$beta==0,]); stopifnot(length(cell$cell)==1L)
for(r in first:last){
  target <- file.path(out,sprintf("r%05d.rds",r))
  if(file.exists(target)) next
  tick <- proc.time()[["elapsed"]]
  set.seed(as.integer(cell$seed0)+r)
  k <- as.integer(cell$k); n <- parse_n(cell$n,k)
  par <- list(cA=cell$cA,cB=cell$cB,sC=cell$sC,cE=cell$cE)
  data <- gen_any(cell$model,cell$beta,n,as.integer(cell$d),par)
  sp <- split_data(data,cell$tau)
  # Preserve the full original learner order and random-number consumption.
  sc <- fit_scores(sp,"ovr",nn_cfg=nn_config(cell$nn),lib=parse_lib(cell$lib))
  labels <- perm_labels(sc$L,as.integer(cell$B))
  tt <- xscr_T(sc$Zfit,sc$gfit,sc$Zte,labels,k,3L)
  v <- apply(tt,2L,max); nb <- length(v)
  p <- pval_perm(v)
  old <- read.csv(file.path(original,"libpilot3",sprintf("%s__r%05d-%05d.csv",cell$cell,r,r)),check.names=FALSE)
  stopifnot(abs(p-old$mca_xscr3)<1e-13, old$seed==as.integer(cell$seed0)+r)
  greater <- sum(v[-1L]>v[1L]+1e-10)
  equal <- sum(abs(v[-1L]-v[1L])<=1e-10)
  stopifnot(abs(p-(1+greater+equal)/nb)<1e-13)
  orbit_p <- vapply(v,function(z) mean(v>=z-1e-10),0)
  # Diagnostic lower bounds only; not proposed p-values for inference.
  no_perm_ties_p <- (1+greater)/nb
  randomized_reject_probability <- max(0,min(1,(0.05*nb-greater)/(equal+1)))
  equivariance_error <- NA_real_
  if(r %in% c(1L,127L,500L)){
    idx <- rev(seq_len(nb))
    permuted <- xscr_T(sc$Zfit,sc$gfit,sc$Zte,labels[,idx],k,3L)
    equivariance_error <- max(abs(permuted-tt[,idx]))
    stopifnot(equivariance_error<1e-10)
  }
  record <- list(rep=r,seed=old$seed,p=p,greater=greater,equal=equal,
    no_perm_ties_p=no_perm_ties_p,orbit_rejections=sum(orbit_p<=0.05),nb=nb,
    randomized_reject_probability=randomized_reject_probability,
    equivariance_error=equivariance_error,v=v)
  saveRDS(record,paste0(target,".tmp"));stopifnot(file.rename(paste0(target,".tmp"),target))
  cat(sprintf("[%s] AUDIT %d/500 | p=%.6f | equal=%d | orbit rejection ranks=%d/%d | %.1fs\n",
    format(Sys.time(),"%H:%M:%S"),r,p,equal,record$orbit_rejections,nb,proc.time()[["elapsed"]]-tick));flush.console()
}
