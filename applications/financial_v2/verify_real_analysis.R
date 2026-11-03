args <- commandArgs(TRUE)
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
folder <- if(length(args))args[1]else file.path(root,"real_analysis")
fingerprint <- readRDS(file.path(folder,"fingerprint.rds"));stopifnot(identical(fingerprint$hashes,tools::md5sum(names(fingerprint$hashes))))
source(file.path(root,"../../R/learners.R"));source(file.path(root,"../../R/tests.R"))
rows <- list()
kruskal_weight <- function(D) {
  pairs <- which(upper.tri(D),arr.ind=TRUE);order <- order(D[pairs]);parent <- seq_len(nrow(D));total <- 0;count <- 0L
  for(z in order) {
    a <- pairs[z,1];b <- pairs[z,2]
    while(parent[a]!=a)a <- parent[a]
    while(parent[b]!=b)b <- parent[b]
    if(a!=b) {parent[b] <- a;total <- total+D[pairs[z,1],pairs[z,2]];count <- count+1L;if(count==nrow(D)-1L)break}
  }
  stopifnot(count==nrow(D)-1L);total
}
for(q in 1:4) {
  saved <- readRDS(file.path(folder,"jobs",sprintf("draw01_q%d.rds",q)))
  for(key in names(saved$audits)) {
    a <- saved$audits[[key]];k <- max(a$records$population);N <- nrow(a$X)
    D <- as.matrix(dist(a$X));K <- exp(-D^2/(2*a$h2));E <- mst_edges(D)
    stopifnot(abs(sum(D[E])-kruskal_weight(D))<1e-9)
    errors <- c(MCA=0,DISCO=0,MMD=0,MST=0)
    for(b in unique(c(1L,2L,17L,ncol(a$Lall)))) {
      groups <- lapply(1:k,function(i)which(a$Lall[,b]==i))
      directT <- numeric(ncol(a$S))
      for(c in seq_len(ncol(a$S))) {
        pos <- a$S[a$Lall[,b]==a$targets[c],c];neg <- a$S[a$Lall[,b]!=a$targets[c],c]
        delta <- outer(pos,neg,"-");auc <- mean(delta>0)+.5*mean(delta==0)
        directT[c] <- (auc-.5)/sqrt((length(pos)+length(neg)+1)/(12*length(pos)*length(neg)))
      }
      disco <- -sum(vapply(groups,function(g)sum(D[g,g,drop=FALSE])/length(g),numeric(1)))
      pairs <- combn(k,2L);md <- numeric(ncol(pairs))
      for(j in seq_len(ncol(pairs))) {
        x <- groups[[pairs[1,j]]];y <- groups[[pairs[2,j]]]
        md[j] <- (sum(K[x,x])-length(x))/(length(x)*(length(x)-1))+
                 (sum(K[y,y])-length(y))/(length(y)*(length(y)-1))-2*mean(K[x,y])
      }
      mst <- sum(a$Lall[E[,1],b]==a$Lall[E[,2],b])
      errors <- pmax(errors,abs(c(max(directT),disco,mean(md),mst)-vapply(a$statistics,`[`,numeric(1),b)))
    }
    stopifnot(max(errors)<1e-9)
    for(method in names(a$statistics)) {
      value <- a$statistics[[method]];p <- (1+sum(value[-1]>=value[1]-1e-10))/length(value)
      r <- saved$rows[saved$rows$stride==as.integer(key)&saved$rows$method==method,]
      stopifnot(nrow(r)==1L,abs(r$p_exploratory-p)<1e-14)
    }
    rows[[length(rows)+1L]] <- data.frame(question=q,stride=as.integer(key),t(errors),MST_weight_independently_verified=TRUE)
  }
}
write.csv(do.call(rbind,rows),file.path(folder,"independent_method_audit.csv"),row.names=FALSE)
cat("All 12 primary-draw question/stride cells passed independent AUC, distance, kernel and MST-weight checks; p-values recomputed. Calibration limitations remain.\n")
