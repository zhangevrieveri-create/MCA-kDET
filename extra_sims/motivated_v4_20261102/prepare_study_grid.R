SRC_DIR <- file.path(getwd(),"R");PROFILE <- "quick";source("R/experiments.R")
source("motivated_v4_20261102/models.R");source("motivated_v4_20261102/models_extra.R")
folder <- "motivated_v4_20261102";stopifnot(file.exists(file.path(folder,"calibration_extra.rds")))
families <- c(MV4_FAMILIES,"N4_mixed",EXTRA_FAMILIES);specs <- list()
for(d in c(20L,100L))for(model in families)for(target in if(d==20L)c(.6,.7,.8) else .7)
  specs[[length(specs)+1L]] <- list(model=model,target=target,d=d,k=if(model=="N4_mixed")4L else 3L)
for(d in c(20L,100L))for(k in 3:4)specs[[length(specs)+1L]] <- list(model=if(k==3L)"N1_bimodal" else "N4_mixed",target=.5,d=d,k=k)
grid <- do.call(rbind,lapply(seq_along(specs),function(i) {
  s <- specs[[i]]
  g <- cell_row("motiv4full",s$model,s$target-.5,k=s$k,n=300L,d=s$d,B=200L,
      reps=if(s$target==.5)500L else 100L,
      groups=if(s$k==4L)"ovr,dist,schott,comp2,comp2ovr,songchen" else LP3_GROUPS,
      lib=LIBPILOT3_LIB,libsets=LIBPILOT3_SETS,seed_scheme="v2")
  g$target <- s$target;g$seed0 <- 160000000L+i*10000L;g
}))
stopifnot(nrow(grid)==28L,sum(grid$reps)==4400L,!anyDuplicated(grid$cell),!anyDuplicated(grid$seed0),
  min(grid$seed0)>PAR$seed+99999*1000,max(grid$seed0)+2500<200000000L)
path <- file.path(folder,"study_grid.rds")
if(file.exists(path))stopifnot(identical(grid,readRDS(path))) else saveRDS(grid,path)
csv <- file.path(folder,"study_grid.csv")
if(!file.exists(csv))write.csv(grid,csv,row.names=FALSE)
cat("FIXED GRID: 28 cells; 4400 formal repetitions; B=200\n")
