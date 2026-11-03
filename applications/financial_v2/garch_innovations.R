# Training-only univariate Gaussian QML. Application layer; frozen E8 unchanged.
innovation_variance <- function(x, p, v0) {
  .C("mca_garch_recursion",as.double(x),as.integer(length(x)),as.double(p[1L]),
     as.double(p[2L]),as.double(p[3L]),as.double(v0),h=double(length(x)),nll=double(1),NAOK=TRUE)
}

innovation_garch_fit <- function(x, fitting, starts=c(.85,.95,.985)) {
  stopifnot(length(x)==length(fitting),any(fitting),!anyNA(fitting))
  mu <- mean(x[fitting],na.rm=TRUE); centered <- x-mu
  xf <- centered[fitting]; v0 <- mean(xf^2,na.rm=TRUE)
  if (!is.finite(v0) || v0<1e-12 || sum(is.finite(xf))<100L) stop("Insufficient nonconstant fitting returns")
  decode <- function(th) {
    persistence <- .9995*plogis(th[2L]); a <- persistence*plogis(th[3L])
    c(omega=exp(th[1L]),alpha=a,beta=persistence-a)
  }
  objective <- function(th) {
    p <- decode(th)
    if (any(!is.finite(p))) return(1e100)
    innovation_variance(xf,p,v0)$nll
  }
  fits <- lapply(starts,function(p) optim(c(log((1-p)*v0),qlogis(p/.9995),qlogis(.1/p)),
    objective,method="Nelder-Mead",control=list(maxit=4000L,reltol=1e-10)))
  converged <- vapply(fits,function(f) f$convergence==0L && is.finite(f$value) && f$value<1e99,logical(1))
  if (!any(converged)) stop("No converged fitting-only GARCH estimate")
  index <- which(converged)[which.min(vapply(fits[converged],`[[`,numeric(1),"value"))]
  best <- fits[[index]]; p <- decode(best$par)
  full <- innovation_variance(centered,p,v0)
  stopifnot(all(is.finite(full$h)),all(full$h>0))
  list(z=centered/sqrt(full$h),h=full$h,
       parameters=c(mean=mu,p,persistence=sum(p[2L:3L]),initial_variance=v0,
                    convergence=best$convergence,negative_log_likelihood=best$value,
                    fitting_observations=sum(is.finite(xf)),failed_starts=sum(!converged),
                    near_boundary=as.integer(p[2L]+p[3L]>.999 || p[2L]<1e-5 || p[3L]<1e-5)),
       optimizer=fits)
}

innovation_smoke <- function() {
  set.seed(37)
  x <- rnorm(800);x[c(25,31,612)] <- NA_real_
  p <- c(.05,.1,.85);v0 <- 1.3
  native <- innovation_variance(x,p,v0)
  h <- numeric(length(x));last_h <- last_x2 <- v0
  for (t in seq_along(x)) {
    h[t] <- p[1L]+p[2L]*last_x2+p[3L]*last_h
    last_h <- h[t]; last_x2 <- if(is.finite(x[t])) x[t]^2 else h[t]
  }
  stopifnot(identical(native$h,h))
  good <- is.finite(x)
  stopifnot(abs(native$nll-.5*sum(log(h[good])+x[good]^2/h[good]))<1e-10)
  fit <- seq_along(x)<=500L
  a <- innovation_garch_fit(x,fit)
  altered <- x;altered[601:800] <- altered[601:800]*4+3
  b <- innovation_garch_fit(altered,fit)
  stopifnot(identical(a$parameters,b$parameters),identical(a$h[1:601],b$h[1:601]),
            identical(is.na(a$z),is.na(x)))
  cat("GARCH smoke: native recursion equals literal R; future perturbations do not change fitting parameters or earlier predictions.\n")
  invisible(TRUE)
}
