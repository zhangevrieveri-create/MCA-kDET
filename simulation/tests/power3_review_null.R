# Explicitly authorized continuation: report unchanged null-size thresholds.
# Does not change any statistic, calibration, source manifest or saved raw row.
args <- commandArgs(TRUE)
argument <- function(name) sub(paste0("^--",name,"="),"",grep(paste0("^--",name,"="),args,value=TRUE)[1L])
out <- argument("out"); id <- argument("cell")
stopifnot(file.exists(file.path(out,"continue_authorization.json")))
source("R/summarize.R")
s <- summarize_exp("libpilot3",out)
z <- s[s$cell==id,]
stopifnot(nrow(z)==1L,z$errors==0L,z$reps==500L)
cols <- c("mca_xscr1","mca_xscr3","mca_g5f")
rates <- unlist(z[cols],use.names=TRUE)
stopifnot(all(is.finite(rates)))
flag <- rates<.031 | rates>.069
review <- file.path(out,"null_level_reviews")
dir.create(review,showWarnings=FALSE)
write.csv(data.frame(cell=id,method=names(rates),rate=unname(rates),lower=.031,upper=.069,
                     outside_range=unname(flag),authorized_continue=TRUE),
          file.path(review,paste0(id,".csv")),row.names=FALSE)
cat("NULL REVIEW",id,paste(sprintf("%s=%.4f",names(rates),rates),collapse=", "),"\n")
if(any(flag)) {
  cat("FLAG RETAINED: outside the original range; continue under explicit user authorization.\n")
} else cat("All three monitored components are within the original range.\n")
