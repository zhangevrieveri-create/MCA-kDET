# Candidate model bootstrap: impose a common innovation law and redo both stages.
# Daily innovation vectors are iid in this candidate model; not an exact theorem.
recursive_seed <- function(case,r,b,q=0L) {
  stopifnot(case>=1L,case<=8L,r>=1L,r<=1000L,b>=0L,b<5000L,q>=0L,q<=4L)
  index <- (((case-1)*1000+(r-1))*5000+b)*5+q
  # Odd multiplier permutes integers modulo 2^31. Products remain exact doubles.
  as.integer((index*65537+123456789)%%2147483648)
}

innovation_pipeline <- function(X,input,d,seeds) {
  fit <- lapply(seq_len(ncol(X)),function(j) innovation_garch_fit(X[,j],input$fitting_days))
  Z <- do.call(cbind,lapply(fit,`[[`,"z"))
  parameters <- data.frame(do.call(rbind,lapply(fit,`[[`,"parameters")),row.names=NULL)
  fitting <- innovation_array(Z,input$fitting_days,d);testing <- innovation_array(Z,input$testing_days,d)
  questions <- lapply(1:4,function(q) {
    members <- if(q<=3L) seq.int((q-1L)*5L+1L,q*5L) else 1:15
    groups <- if(q<=3L) 1:5 else rep(1:3,each=5L)
    panel_question(fitting,testing,input$masks$fit,input$masks$test,members,groups,seeds[q])
  })
  statistic <- vapply(questions,function(question) {
    counts <- matrix(1L,length(input$masks$test),1L)
    auc <- panel_auc(question$S,question$groups,question$targets,question$time_id,counts)[,1L]
    component <- (auc-.5)/question$scale
    original <- ovr_T(question$S,question$targets,matrix(question$groups,ncol=1L))[,1L]
    stopifnot(max(abs(component-original))<1e-10)
    max(component)
  },numeric(1))
  list(Z=Z,parameters=parameters,questions=questions,statistic=statistic)
}

innovation_null_model <- function(observed,input) {
  Z <- observed$Z;complete <- apply(is.finite(Z),1L,all)
  fit <- Z[complete&input$fitting_days,,drop=FALSE]
  test <- Z[complete&input$testing_days,,drop=FALSE]
  stopifnot(min(nrow(fit),nrow(test))>100L)
  values <- as.vector(rbind(fit,test))
  values <- sort((values-mean(values))/sqrt(mean((values-mean(values))^2)))
  copula <- function(M) list(lower=apply(M,2L,rank,ties.method="min")-1L,
                            span=apply(M,2L,rank,ties.method="max")-apply(M,2L,rank,ties.method="min")+1L)
  list(common_law=values,fit=copula(fit),test=copula(test),parameters=observed$parameters,
       hypothesis="Common unknown 1D innovation law with iid daily vectors; joint empirical copulas by period",
       initialization="Conditional recursion from original fitting variance, with no stationary burn-in")
}

innovation_null_draw <- function(model,input,seed) {
  set.seed(seed);N <- nrow(input$valid);J <- ncol(input$valid)
  E <- matrix(NA_real_,N,J)
  for(part in c("fit","test")) {
    days <- if(part=="fit") which(input$fitting_days) else which(!input$fitting_days)
    source <- model[[part]];n <- nrow(source$lower)
    indices <- sample.int(n,length(days),replace=TRUE)
    U <- (source$lower[indices,,drop=FALSE]+matrix(runif(length(days)*J),length(days),J)*source$span[indices,,drop=FALSE])/n
    quantile <- pmax(1L,pmin(length(model$common_law),as.integer(ceiling(U*length(model$common_law)))))
    E[days,] <- matrix(model$common_law[quantile],length(days),J)
  }
  p <- model$parameters;X <- matrix(0,N,J);h <- previous_x2 <- p$initial_variance
  for(t in seq_len(N)) {
    h <- p$omega+p$alpha*previous_x2+p$beta*h
    centered <- sqrt(h)*E[t,];X[t,] <- p$mean+centered;previous_x2 <- centered^2
  }
  X[!input$valid] <- NA_real_
  list(X=X,E=E)
}

recursive_model_smoke <- function() {
  set.seed(73);M <- matrix(sample(-2:2,15*500,replace=TRUE),500,15)
  copula <- list(lower=apply(M,2,rank,ties.method="min")-1,
                 span=apply(M,2,rank,ties.method="max")-apply(M,2,rank,ties.method="min")+1)
  # Repeated tied intervals have weight equal to their width, giving uniform U.
  for(j in 1:15) {
    low <- copula$lower[,j];span <- copula$span[,j]
    for(a in unique(low)) stopifnot(sum(low==a)==unique(span[low==a]))
  }
  seeds <- unlist(lapply(1:6,function(case) lapply(1:20,function(r) sapply(0:99,function(b) sapply(0:4,function(q) recursive_seed(case,r,b,q))))))
  stopifnot(!anyDuplicated(seeds),all(seeds>=0))
  cat("Recursive model smoke: tie-rank copula intervals and distinct scalar seed mapping verified.\n")
}

recursive_job <- function(job,config,input) {
  case <- match(job$scenario,INNOVATION_SCENARIOS);r <- job$rep
  folder <- file.path(config$outdir,job$scenario,sprintf("r%05d",r))
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  final <- file.path(folder,sprintf("result_B%04d.csv",config$B))
  if(file.exists(final)) return(invisible(TRUE))
  observed_path <- file.path(folder,"observed.rds")
  if(file.exists(observed_path)) observed <- readRDS(observed_path) else {
    generated <- innovation_generate(job$scenario,nrow(input$valid),recursive_seed(case,r,0L,0L))
    generated$X[!input$valid] <- NA_real_
    observed <- innovation_pipeline(generated$X,input,config$d,sapply(1:4,function(q) recursive_seed(case,r,0L,q)))
    saveRDS(observed,observed_path)
  }
  model_path <- file.path(folder,"null_model.rds")
  if(file.exists(model_path)) model <- readRDS(model_path) else {model <- innovation_null_model(observed,input);saveRDS(model,model_path)}
  inner <- list()
  for(first in seq.int(1L,config$B,by=5L)) {
    last <- min(first+4L,4999L)
    path <- file.path(folder,sprintf("inner_b%04d-%04d.csv",first,last))
    if(file.exists(path)) tab <- read.csv(path) else {
      tab <- do.call(rbind,lapply(seq.int(first,last),function(b) {
        value <- NULL
        for(attempt in 1:3) {
          value <- tryCatch({
            draw <- innovation_null_draw(model,input,recursive_seed(case,r,b,0L))
            pipeline <- innovation_pipeline(draw$X,input,config$d,sapply(1:4,function(q) recursive_seed(case,r,b,q)))
            data.frame(b=b,within_bank=pipeline$statistic[1L],within_pharma=pipeline$statistic[2L],
                       within_machinery=pipeline$statistic[3L],between_sectors=pipeline$statistic[4L],
                       filter_near_boundary=sum(pipeline$parameters$near_boundary),
                       learner_nonconvergence=sum(vapply(pipeline$questions,`[[`,numeric(1),"nonconv")))
          },error=function(e)e)
          if(!inherits(value,"error")) break
        }
        if(inherits(value,"error")) stop(job$scenario," r=",r," b=",b,": ",conditionMessage(value))
        value
      }))
      write.csv(tab,paste0(path,".part"),row.names=FALSE);stopifnot(file.rename(paste0(path,".part"),path))
    }
    inner[[length(inner)+1L]] <- tab
    cat(sprintf("Recursive %s r=%d: %d/%d null datasets, filters + classifiers re-estimated\n",job$scenario,r,last,config$B));flush.console()
  }
  inner <- do.call(rbind,inner);inner <- inner[inner$b<=config$B,]
  stopifnot(nrow(inner)==config$B,!anyDuplicated(inner$b))
  statistic <- as.matrix(inner[,c("within_bank","within_pharma","within_machinery","between_sectors")])
  p <- (1+colSums(sweep(statistic,2L,observed$statistic,"-")>=-1e-10))/(config$B+1)
  result <- data.frame(scenario=job$scenario,rep=r,d=config$d,B=config$B,question=colnames(statistic),
                       statistic=observed$statistic,p=p,p_holm=p.adjust(p,"holm"),
                       filter_near_boundary=sum(observed$parameters$near_boundary))
  write.csv(result,paste0(final,".part"),row.names=FALSE);stopifnot(file.rename(paste0(final,".part"),final))
  cat(sprintf("Completed recursive innovation pilot %s r=%d\n",job$scenario,r));flush.console()
  invisible(TRUE)
}
