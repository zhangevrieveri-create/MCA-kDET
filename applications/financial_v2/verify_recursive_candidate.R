# Reproduce one full inner pipeline and check the actual null-generating law.
args <- commandArgs(TRUE)
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),winslash="/"))
folder <- if(length(args)) args[1L] else file.path(root,"../../tests/out/financial_recursive_smoke/null_garch99/r00001")
for(pkg in c("glmnet","ranger")) stopifnot(requireNamespace(pkg,quietly=TRUE))
for(f in c("learners.R","learners_lib.R","tests.R")) source(file.path(root,"../../R",f))
for(f in c("panel_auc.dll","garch_filter.dll")) dyn.load(file.path(root,"native",f))
for(f in c("panel_calibration.R","garch_innovations.R","innovation_validation.R","recursive_innovation_bootstrap.R")) source(file.path(root,f))
install_panel_knn()
real <- readRDS(file.path(root,"innovations/filters.rds"))
cal <- read.csv(file.path(root,"design_two_tests/block_calendar.csv"));cal <- cal[cal$dimension==5L,]
valid <- tolower(as.character(cal$complete_all_main_stocks))=="true"
input <- list(valid=is.finite(real$X),fitting_days=real$fitting,testing_days=real$dates>=as.Date("2018-01-30"),
              masks=list(fit=valid[cal$split=="fit"],test=valid[cal$split=="test"]))
model <- readRDS(file.path(folder,"null_model.rds"))
stopifnot(abs(mean(model$common_law))<1e-12,abs(mean(model$common_law^2)-1)<1e-12)
for(part in c("fit","test")) for(j in 1:15) {
  low <- model[[part]]$lower[,j];span <- model[[part]]$span[,j]
  starts <- sort(unique(low));ends <- numeric(length(starts))
  for(a in seq_along(starts)) {
    selected <- low==starts[a]
    stopifnot(length(unique(span[selected]))==1L,sum(selected)==unique(span[selected]))
    ends[a] <- starts[a]+unique(span[selected])
  }
  stopifnot(starts[1L]==0,tail(ends,1L)==length(low),all(head(ends,-1L)==tail(starts,-1L)))
}
case <- match("null_garch99",INNOVATION_SCENARIOS)
draw <- innovation_null_draw(model,input,recursive_seed(case,1L,1L,0L))
repeated <- innovation_null_draw(model,input,recursive_seed(case,1L,1L,0L))
stopifnot(identical(draw,repeated),all(is.na(draw$X)==!input$valid))
pipeline <- innovation_pipeline(draw$X,input,5L,sapply(1:4,function(q) recursive_seed(case,1L,1L,q)))
saved <- read.csv(file.path(folder,"inner_b0001-0003.csv"))
expected <- as.numeric(saved[saved$b==1L,c("within_bank","within_pharma","within_machinery","between_sectors")])
error <- max(abs(pipeline$statistic-expected));stopifnot(error<1e-12)
write.csv(data.frame(common_law_mean=mean(model$common_law),common_law_second_moment=mean(model$common_law^2),
                     fit_test_copula_rank_partitions_verified=TRUE,missingness_preserved=TRUE,
                     repeated_draw_bit_identical=TRUE,full_pipeline_statistic_reproduction_error=error),
          file.path(folder,"independent_verification.csv"),row.names=FALSE)
cat(sprintf("Recursive candidate verified: common null marginals, missingness, deterministic draws, full-pipeline reproduction error %.3g.\n",error))
