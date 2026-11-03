# Proposed method-motivated scenarios. All have E[X]=0 and Cov(X)=I.
# Existing data generators and calibration constants are untouched.
MV4_FAMILIES <- c("N1_bimodal","N2_scale_mix","N3_condvar")
MV4_TARGETS <- c(.60,.70,.80)
MV4_S <- 3L

mv4_noise <- function(n=200000L, seed=20261102L) {
  had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had) previous <- .Random.seed
  on.exit(if(had) assign(".Random.seed",previous,.GlobalEnv) else rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed)
  list(Zd=matrix(rnorm(n*6L),n,6L),Zr=matrix(rnorm(n*6L),n,6L),
       signs=matrix(sample(c(-1,1),n*3L,replace=TRUE),n,3L))
}

mv4_draw <- function(family,theta,Z,signs=NULL) {
  if(family %in% c("N1_bimodal","N2_scale_mix")) {
    X <- Z[,1:3,drop=FALSE]
    if(is.null(signs))signs <- matrix(sample(c(-1,1),nrow(X)*3L,replace=TRUE),nrow(X),3L)
    if(family=="N1_bimodal")X <- theta*signs+sqrt(1-theta^2)*X
    else X <- sqrt(1+theta*signs)*X
    return(X)
  }
  stopifnot(family=="N3_condvar")
  X <- Z[,1:6,drop=FALSE]
  for(j in 1:3) {
    a <- 2L*j-1L;b <- 2L*j
    X[,b] <- sqrt(1-theta+theta*X[,a]^2)*X[,b]
  }
  X
}

mv4_logratio <- function(family,theta,X) {
  if(theta==0)return(rep(0,nrow(X)))
  if(family=="N1_bimodal") {
    X <- X[,1:3,drop=FALSE];sd <- sqrt(1-theta^2)
    return(rowSums(logsumexp2(dnorm(X,theta,sd,log=TRUE),dnorm(X,-theta,sd,log=TRUE))-log(2)-dnorm(X,log=TRUE)))
  }
  if(family=="N2_scale_mix") {
    X <- X[,1:3,drop=FALSE]
    return(rowSums(logsumexp2(dnorm(X,sd=sqrt(1+theta),log=TRUE),dnorm(X,sd=sqrt(1-theta),log=TRUE))-log(2)-dnorm(X,log=TRUE)))
  }
  ans <- numeric(nrow(X))
  for(j in 1:3) {
    a <- 2L*j-1L;b <- 2L*j;v <- 1-theta+theta*X[,a]^2
    ans <- ans-.5*log(v)+.5*X[,b]^2*(1-1/v)
  }
  ans
}

mv4_oracle <- function(family,theta,noise) {
  if(theta==0)return(.5)
  X <- mv4_draw(family,theta,noise$Zd,noise$signs)
  mc_auc(mv4_logratio(family,theta,X),mv4_logratio(family,theta,noise$Zr))
}

mv4_generate <- function(model,target,n,d,calibration) {
  k <- length(n);stopifnot(d>=12L)
  mixed <- model=="N4_mixed"
  if(mixed)stopifnot(k==4L)
  lapply(seq_len(k),function(group) {
    X <- matrix(rnorm(n[group]*d),n[group],d)
    if(target<=.5 || (mixed && group==1L) || (!mixed && group!=k))return(X)
    family <- if(mixed)MV4_FAMILIES[group-1L] else model
    cols <- if(!mixed)seq_len(if(family=="N3_condvar")6L else 3L) else switch(as.character(group),"2"=1:3,"3"=4:6,"4"=7:12)
    theta <- calibration[[family]][[sprintf("%.2f",target)]]
    stopifnot(is.finite(theta),theta>0,theta<1)
    Z <- matrix(rnorm(n[group]*6L),n[group],6L)
    X[,cols] <- mv4_draw(family,theta,Z)
    X
  })
}
