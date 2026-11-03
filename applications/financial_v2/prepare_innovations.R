# Per-stock GARCH adjustment authorized after the original-return pilot failed.
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
dyn.load(file.path(root,"native/garch_filter.dll"))
source(file.path(root,"garch_innovations.R"))
innovation_smoke()
out <- file.path(root,"innovations");dir.create(out,showWarnings=FALSE)
panel <- read.csv(file.path(root,"design_two_tests/main_panel.csv"),check.names=FALSE)
input <- read.csv(gzfile(file.path(root,"design_two_tests/daily_returns.csv.gz")),check.names=FALSE)
input <- input[input$date>="2015-01-01",]
dates <- as.Date(input$date)
X <- as.matrix(input[,panel$order_book_id,drop=FALSE])
fitting <- dates<=as.Date("2017-12-31")
stopifnot(all(fitting==(seq_along(fitting)<=sum(fitting))))
fits <- vector("list",ncol(X))
for(j in seq_len(ncol(X))) {
  fits[[j]] <- innovation_garch_fit(X[,j],fitting)
  cat(sprintf("Fitted training-only daily GARCH %d/%d: %s\n",j,ncol(X),colnames(X)[j]));flush.console()
}
Z <- do.call(cbind,lapply(fits,`[[`,"z"));H <- do.call(cbind,lapply(fits,`[[`,"h"))
colnames(Z) <- colnames(H) <- colnames(X)
stopifnot(identical(is.na(Z),is.na(X)))
parameters <- data.frame(stock=colnames(X),do.call(rbind,lapply(fits,`[[`,"parameters")),row.names=NULL)
write.csv(parameters,file.path(out,"garch_parameters.csv"),row.names=FALSE)
write_exact <- function(M,path) {
  tab <- data.frame(date=as.character(dates),check.names=FALSE)
  for(j in seq_len(ncol(M))) tab[[colnames(M)[j]]] <- ifelse(is.finite(M[,j]),sprintf("%.17g",M[,j]),NA_character_)
  con <- gzfile(path,"wt");on.exit(close(con))
  write.csv(tab,con,row.names=FALSE,na="")
}
write_exact(Z,file.path(out,"daily_innovations.csv.gz"))
write_exact(H,file.path(out,"conditional_variance.csv.gz"))
pair_acf <- function(v,lag) {
  a <- v[seq_len(length(v)-lag)];b <- v[seq.int(lag+1L,length(v))]
  valid <- is.finite(a)&is.finite(b)
  if(sum(valid)<10L) return(NA_real_)
  cor(a[valid],b[valid])
}
rows <- list()
testing <- dates>=as.Date("2018-01-30")
for(j in seq_len(ncol(X))) for(representation in c("raw","filtered")) for(power in 1:2) for(lag in c(1L,5L,20L)) {
  v <- if(representation=="raw") X[testing,j] else Z[testing,j]
  rows[[length(rows)+1L]] <- data.frame(stock=colnames(X)[j],representation=representation,power=power,lag=lag,
                                       autocorrelation=pair_acf(v^power,lag))
}
write.csv(do.call(rbind,rows),file.path(out,"autocorrelation.csv"),row.names=FALSE)
complete <- apply(is.finite(X[testing,,drop=FALSE]),1,all)
raw_cor <- cor(X[testing,,drop=FALSE][complete,,drop=FALSE]);filtered_cor <- cor(Z[testing,,drop=FALSE][complete,,drop=FALSE])
write.csv(raw_cor,file.path(out,"daily_correlation_raw.csv"))
write.csv(filtered_cor,file.path(out,"daily_correlation_filtered.csv"))
saveRDS(list(fits=fits,parameters=parameters,dates=dates,X=X,Z=Z,H=H,fitting=fitting,session=sessionInfo()),file.path(out,"filters.rds"))
print(parameters[,c("stock","persistence","convergence","failed_starts","near_boundary")],row.names=FALSE)
cat(sprintf("Prepared GARCH innovations: fitting days=%d, buffer days=%d, testing days=%d. No classifiers or p-values.\n",sum(fitting),sum(!fitting&!testing),sum(testing)))
