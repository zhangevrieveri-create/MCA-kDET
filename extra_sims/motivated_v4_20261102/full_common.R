# New isolated experiment wrapper. Existing generators and tests are read only.
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (p in c("glmnet", "ranger", "gTestsMulti")) stopifnot(requireNamespace(p, quietly=TRUE))
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
MOTIV_DIR <- "motivated_v4_20261102"
source(file.path(MOTIV_DIR,"models.R")); source("R/oracle_s.R")
source(file.path(MOTIV_DIR,"models_extra.R"))
MOTIV_CAL <- readRDS(file.path(MOTIV_DIR,"calibration.rds"))
MOTIV_CAL_EXTRA <- readRDS(file.path(MOTIV_DIR,"calibration_extra.rds"))

motiv_arg <- function(key, default=NULL) {
  a <- grep(paste0("^--", key, "="), commandArgs(TRUE), value=TRUE)
  if (!length(a)) return(default)
  stopifnot(length(a)==1L); sub(paste0("^--", key, "="), "", a)
}
motiv_manifest_check <- function(out) {
  z <- read.csv(file.path(out,"source_md5.csv"), stringsAsFactors=FALSE)
  actual <- as.character(tools::md5sum(z$file))
  if (any(is.na(actual)) || !identical(actual,z$md5))
    stop("Source fingerprint mismatch: ", paste(z$file[is.na(actual)|actual!=z$md5],collapse=", "))
  invisible(TRUE)
}
motiv_xscr_diagnostic <- function(sc,m) {
  d <- ncol(sc$Zfit)
  pairs <- if(d<=XSCR_MAX_D_CROSS && d>=2L)which(upper.tri(diag(d)),arr.ind=TRUE) else matrix(integer(),0L,2L)
  F <- xscr_features(sc$Zfit,pairs); T <- xscr_features(sc$Zte,pairs)
  labels <- c(paste0("x",seq_len(d)),paste0("x",seq_len(d),"^2"),
              if(nrow(pairs))paste0("x",pairs[,1],"*x",pairs[,2]))
  ans <- list()
  for(direction in 1:2) {
    A <- if(direction==1L)F else T; C <- if(direction==1L)T else F
    g <- if(direction==1L)sc$gfit else sc$L; h <- if(direction==1L)sc$L else sc$gfit
    mu <- colMeans(A); ss <- colsd(A); ss[!is.finite(ss)|ss<1e-12] <- 1
    Az <- sweep(sweep(A,2L,mu),2L,ss,"/"); Cz <- sweep(sweep(C,2L,mu),2L,ss,"/")
    for(i in seq_len(sc$k)) {
      pos <- g==i; n1 <- sum(pos); n0 <- length(g)-n1; s1 <- colSums(Az[pos,,drop=FALSE])
      tt <- (s1/n1-(colSums(Az)-s1)/n0)/sqrt(1/n1+1/n0)
      sel <- top_m(abs(tt),min(m,length(tt)))
      score <- as.numeric(Cz[,sel,drop=FALSE]%*%tt[sel]); auc <- auc_vec(score,h==i)
      ans[[length(ans)+1L]] <- data.frame(m=m,direction=direction,population=i,
        selected=paste(labels[sel],collapse=";"),auc=auc,
        T=(auc-.5)/sqrt((length(h)+1)/(12*sum(h==i)*sum(h!=i))))
    }
  }
  ans <- do.call(rbind,ans)
  actual <- xscr_T(sc$Zfit,sc$gfit,sc$Zte,matrix(sc$L,ncol=1L),sc$k,m)[,1L]
  rebuilt <- vapply(seq_len(sc$k),function(i)sum(ans$T[ans$population==i])/sqrt(2),numeric(1))
  stopifnot(max(abs(actual-rebuilt))<1e-10)
  ans
}
motiv_one <- function(cell,r) {
  seed <- as.integer(cell$seed0+r); set.seed(seed); t0 <- proc.time()[["elapsed"]]
  k <- as.integer(cell$k); B <- as.integer(cell$B); lib <- parse_lib(cell$lib)
  nvec <- parse_n(cell$n,k)
  data <- study_generate(cell$model,cell$target,nvec,as.integer(cell$d),MOTIV_CAL,MOTIV_CAL_EXTRA)
  stopifnot(length(data)==k,identical(vapply(data,nrow,integer(1)),nvec),
            all(vapply(data,ncol,integer(1))==cell$d))
  sc <- fit_scores(split_data(data,cell$tau),"ovr",lib=lib)
  ct <- classifier_tests(sc,B,parse_libsets(cell$libsets,lib))
  comp <- competitor_tests(data,B,strsplit(cell$groups,",",fixed=TRUE)[[1L]])
  p <- c(ct$p,comp$p[setdiff(names(comp$p),"kmd_eta")])
  effect <- comp$p[intersect(names(comp$p),"kmd_eta")]
  pv <- unlist(p); stopifnot(!anyDuplicated(names(p)),all(is.finite(pv)),all(pv>=0&pv<=1))
  aucs <- do.call(rbind,lapply(names(sc$Sovr),function(a)data.frame(learner=a,population=seq_len(k),
    auc=vapply(seq_len(k),function(i)auc_vec(sc$Sovr[[a]][,i],sc$L==i),numeric(1)))))
  xscr <- rbind(motiv_xscr_diagnostic(sc,1L),motiv_xscr_diagnostic(sc,3L))
  list(cell=cell$cell,model=cell$model,target=cell$target,rep=as.integer(r),seed=seed,
       B=B,k=k,n=cell$n,d=cell$d,tau=cell$tau,lib=cell$lib,libsets=cell$libsets,
       p=p,effect=effect,diag=ct$diag,aucs=aucs,xscr=xscr,err="",
       timing=c(sc$tm,comp$tm),seconds=proc.time()[["elapsed"]]-t0)
}
