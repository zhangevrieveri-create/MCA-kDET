# Additional exploratory models; original N1--N4 and frozen R/models.R untouched.
EXTRA_FAMILIES <- c("N5_skewmix","N6_ring")

extra_noise <- function(n=200000L,seed=20261112L) {
  had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had)previous <- .Random.seed
  on.exit(if(had)assign(".Random.seed",previous,.GlobalEnv) else rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed)
  list(Zd=matrix(rnorm(n*6L),n,6L),Zr=matrix(rnorm(n*6L),n,6L),U=matrix(runif(n*3L),n,3L))
}

extra_draw <- function(family,theta,Z,U=NULL) {
  stopifnot(family%in%EXTRA_FAMILIES,is.finite(theta),theta>=0,theta<1)
  n <- nrow(Z)
  if(is.null(U))U <- matrix(runif(n*3L),n,3L)
  sdv <- sqrt(1-theta^2)
  if(family=="N5_skewmix")return(theta*ifelse(U<.2,2,-.5)+sdv*Z[,1:3,drop=FALSE])
  X <- Z[,1:6,drop=FALSE]*sdv
  for(j in 1:3) {
    angle <- 2*pi*U[,j]
    X[,2*j-1] <- X[,2*j-1]+sqrt(2)*theta*cos(angle)
    X[,2*j] <- X[,2*j]+sqrt(2)*theta*sin(angle)
  }
  X
}

extra_logratio <- function(family,theta,X) {
  if(theta==0)return(rep(0,nrow(X)))
  v <- 1-theta^2
  if(family=="N5_skewmix") {
    X <- X[,1:3,drop=FALSE]
    return(rowSums(logsumexp2(log(.2)+dnorm(X,2*theta,sqrt(v),log=TRUE),
                             log(.8)+dnorm(X,-.5*theta,sqrt(v),log=TRUE))-dnorm(X,log=TRUE)))
  }
  stopifnot(family=="N6_ring")
  ans <- numeric(nrow(X))
  for(j in 1:3) {
    r2 <- X[,2*j-1]^2+X[,2*j]^2
    z <- sqrt(r2)*sqrt(2)*theta/v
    logI0 <- log(besselI(z,nu=0,expon.scaled=TRUE))+z
    ans <- ans-log(v)-(r2+2*theta^2)/(2*v)+logI0+r2/2
  }
  ans
}

extra_oracle <- function(family,theta,noise) {
  if(theta==0)return(.5)
  X <- extra_draw(family,theta,noise$Zd,noise$U)
  mc_auc(extra_logratio(family,theta,X),extra_logratio(family,theta,noise$Zr))
}

study_generate <- function(model,target,n,d,calibration,calibration_extra) {
  if(!model%in%EXTRA_FAMILIES)return(mv4_generate(model,target,n,d,calibration))
  stopifnot(d>=6L,length(n)>=2L)
  lapply(seq_along(n),function(i) {
    X <- matrix(rnorm(n[i]*d),n[i],d)
    if(target<=.5||i!=length(n))return(X)
    theta <- calibration_extra[[model]][[sprintf("%.2f",target)]]
    stopifnot(length(theta)==1L,is.finite(theta),theta>0,theta<1)
    Z <- matrix(rnorm(n[i]*6L),n[i],6L)
    cols <- seq_len(if(model=="N5_skewmix")3L else 6L)
    X[,cols] <- extra_draw(model,theta,Z)
    X
  })
}
