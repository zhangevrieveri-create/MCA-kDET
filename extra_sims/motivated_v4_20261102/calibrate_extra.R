source("R/oracle_s.R");source("motivated_v4_20261102/models_extra.R")
folder <- "motivated_v4_20261102";options(power3.progress=TRUE)
cache <- file.path(folder,"calibration_extra.rds")
cal <- if(file.exists(cache))readRDS(cache) else setNames(vector("list",2L),EXTRA_FAMILIES)
noise <- extra_noise()
for(f in EXTRA_FAMILIES)for(a in c(.60,.70,.80)) {
  key <- sprintf("%.2f",a)
  if(!is.null(cal[[f]][[key]]))next
  cat("CALIBRATE",f,a,"\n");flush.console()
  fn <- function(theta)extra_oracle(f,theta,noise)
  stopifnot(fn(.995)>=a)
  th <- oracle_solve(fn,a,0,.995,tol=1e-5)
  stopifnot(abs(fn(th)-a)<.001)
  cal[[f]][[key]] <- th;saveRDS(cal,cache)
}
ind <- extra_noise(100000L,20261113L)
check <- do.call(rbind,lapply(EXTRA_FAMILIES,function(f)do.call(rbind,lapply(c(.60,.70,.80),function(a) {
  th <- cal[[f]][[sprintf("%.2f",a)]];actual <- extra_oracle(f,th,ind)
  stopifnot(abs(actual-a)<.006)
  data.frame(family=f,target=a,theta=th,independent_auc=actual)
}))))
write.csv(check,file.path(folder,"calibration_extra_check.csv"),row.names=FALSE)
set.seed(20261114L)
for(f in EXTRA_FAMILIES) {
  Z <- matrix(rnorm(200000L*6L),200000L,6L)
  X <- extra_draw(f,cal[[f]][["0.70"]],Z)
  mean_error <- max(abs(colMeans(X)));cov_error <- max(abs(cov(X)-diag(ncol(X))))
  stopifnot(mean_error<.025,cov_error<.04)
  cat("MOMENTS PASS",f,mean_error,cov_error,"\n")
}
# Independent density calculation by direct mixture / angular quadrature.
x <- seq(-3,3,length.out=31);theta <- .7;sdv <- sqrt(1-theta^2)
X <- cbind(x,0*x,0*x)
direct <- log(.2*dnorm(x,2*theta,sdv)+.8*dnorm(x,-.5*theta,sdv))-dnorm(x,log=TRUE)
zero <- log(.2*dnorm(0,2*theta,sdv)+.8*dnorm(0,-.5*theta,sdv))-dnorm(0,log=TRUE)
stopifnot(max(abs(extra_logratio("N5_skewmix",theta,X)-direct-2*zero))<1e-10)
for(point in list(c(0,0),c(.4,-.9),c(2.1,.3))) {
  integrand <- function(a)dnorm(point[1],sqrt(2)*theta*cos(a),sdv)*dnorm(point[2],sqrt(2)*theta*sin(a),sdv)/(2*pi)
  q <- integrate(integrand,0,2*pi,rel.tol=1e-10)$value
  zz <- matrix(rep(point,3L),nrow=1L)
  actual <- extra_logratio("N6_ring",theta,zz)/3
  expected <- log(q)-sum(dnorm(point,log=TRUE))
  stopifnot(abs(actual-expected)<1e-9)
}
writeLines("All six roots independently validated; moments match; skew mixture density and ring angular-integral density agree with independent formulas.",file.path(folder,"extra_model_verification.txt"))
cat("EXTRA CALIBRATION AND DENSITY CHECKS COMPLETE\n")
