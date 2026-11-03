source("R/oracle_s.R")
source("motivated_v4_20261102/models.R")
options(power3.progress=TRUE)
folder <- "motivated_v4_20261102"
cache <- file.path(folder,"calibration.rds")
calibration <- if(file.exists(cache))readRDS(cache) else setNames(vector("list",3L),MV4_FAMILIES)
nz <- mv4_noise()
for(family in MV4_FAMILIES)for(target in MV4_TARGETS) {
  key <- sprintf("%.2f",target)
  if(!is.null(calibration[[family]][[key]]))next
  cat("CALIBRATION",family,"target",target,"\n");flush.console()
  f <- function(theta)mv4_oracle(family,theta,nz)
  endpoint <- f(.995)
  if(endpoint<target)stop("Unreachable fixed target: ",family," ",target," (max=",endpoint,")")
  theta <- oracle_solve(f,target,0,.995,tol=1e-5)
  stopifnot(abs(f(theta)-target)<.001)
  calibration[[family]][[key]] <- theta
  saveRDS(calibration,cache)
}
set.seed(20261103L)
moments <- list()
for(family in MV4_FAMILIES) {
  Z <- matrix(rnorm(200000L*6L),200000L,6L)
  X <- mv4_draw(family,calibration[[family]][["0.70"]],Z)
  err <- max(abs(cov(X)-diag(ncol(X))));mean_err <- max(abs(colMeans(X)))
  stopifnot(err<.04,mean_err<.02)
  moments[[family]] <- c(max_abs_mean=mean_err,max_abs_covariance_error=err)
  cat("MOMENTS OK",family,"mean",mean_err,"covariance",err,"\n")
}
saveRDS(moments,file.path(folder,"moment_checks.rds"))
cat("CALIBRATION AND MOMENT CHECKS COMPLETE\n")
