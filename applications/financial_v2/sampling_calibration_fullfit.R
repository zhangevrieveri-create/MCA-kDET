# Full fitting panel is allowed: oracle validity conditions on arbitrary D_fit.
# Only held-out testing records need the one-stock-per-time sampling design.
sampling_fit <- function(fitting,testing,fit_records,test_valid,q,seed) {
  k <- if(q<=3L) 5L else 3L
  members <- if(q<=3L) seq.int((q-1L)*5L+1L,q*5L) else 1:15
  groups <- if(q<=3L) 1:5 else rep(1:3,each=5L)
  fit_valid <- apply(is.finite(fitting),1L,all)
  all_fit <- data.frame(time=rep(which(fit_valid),length(members)),stock=rep(members,each=sum(fit_valid)),
                        population=rep(groups,each=sum(fit_valid)))
  fits <- sampling_matrices(fitting,all_fit,k)
  records <- lapply(1:3,function(h) {
    kind <- 20L+(q-1L)*3L+h
    schedule_seed <- as.integer((as.double(seed)+(kind-(50L+q))*65537)%%2147483648)
    sampling_records(test_valid,q,c(1L,2L,4L)[h],schedule_seed)
  })
  bag <- do.call(rbind,records);tests <- sampling_matrices(testing,bag,k)
  stopifnot(length(unique(tabulate(bag$population,k)))==1L)
  sp <- lapply(seq_len(k),function(g)list(fit=fits[[g]],test=tests[[g]]))
  set.seed(seed);warnings_seen <- character()
  sc <- withCallingHandlers(fit_scores(sp,"ovr",lib=c("l1","quad","knn","rf")),
    warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  order <- unlist(lapply(seq_len(k),function(g)which(bag$population==g)))
  list(S=do.call(cbind,sc$Sovr),targets=rep(seq_len(k),4L),labels=sc$L,k=k,sp=sp,
       time=bag$time[order],stock=bag$stock[order],
       warning_count=length(unique(warnings_seen)),nonconv=sum(unlist(lapply(sc$info_ovr,function(o)o$nonconv))))
}

sampling_rep_sparse_fit <- sampling_rep
sampling_rep <- function(scenario,r,config,input,keep=FALSE) {
  result <- sampling_rep_sparse_fit(scenario,r,config,input,keep)
  for(phase in names(result$rows)) {
    tab <- result$rows[[phase]]
    tab$fitting_rows <- ifelse(tab$question=="between_sectors",15L,5L)*sum(input$masks$fit)
    result$rows[[phase]] <- tab
    if(keep) for(key in names(result$audit[[phase]]$qa)) {
      q <- as.integer(strsplit(key,"_",fixed=TRUE)[[1L]][1L])
      members <- if(q<=3L) seq.int((q-1L)*5L+1L,q*5L) else 1:15
      groups <- if(q<=3L) 1:5 else rep(1:3,each=5L)
      result$audit[[phase]]$qa[[key]]$fitting_records <- data.frame(time=rep(which(input$masks$fit),length(members)),
        stock=rep(members,each=sum(input$masks$fit)),population=rep(groups,each=sum(input$masks$fit)))
    }
  }
  result
}
