# Computational batching only: predict requested test records, keeping omega=1/k.
# sampling_seed is affine modulo 2^31, so test schedules can be reconstructed
# from the classifier seed without changing any RNG draw used for fitting.
sampling_fit <- function(fitting,testing,fit_records,test_valid,q,seed) {
  k <- if(q<=3L) 5L else 3L
  records <- lapply(1:3,function(h) {
    kind <- 20L+(q-1L)*3L+h
    schedule_seed <- as.integer((as.double(seed)+(kind-(50L+q))*65537)%%2147483648)
    sampling_records(test_valid,q,c(1L,2L,4L)[h],schedule_seed)
  })
  bag <- do.call(rbind,records)
  stopifnot(length(unique(tabulate(bag$population,k)))==1L)
  fits <- sampling_matrices(fitting,fit_records,k)
  tests <- sampling_matrices(testing,bag,k)
  sp <- lapply(seq_len(k),function(g) list(fit=fits[[g]],test=tests[[g]]))
  set.seed(seed);warnings_seen <- character()
  sc <- withCallingHandlers(fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf")),
    warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  order <- unlist(lapply(seq_len(k),function(g)which(bag$population==g)))
  list(S=do.call(cbind,sc$Sovr),targets=rep(seq_len(k),4L),labels=sc$L,k=k,sp=sp,
       time=bag$time[order],stock=bag$stock[order],
       warning_count=length(unique(warnings_seen)),nonconv=sum(unlist(lapply(sc$info_ovr,function(o)o$nonconv))))
}
