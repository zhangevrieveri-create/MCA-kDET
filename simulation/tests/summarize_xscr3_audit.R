root <- "tests/out/xscr3_audit_20261101"
fs <- list.files(file.path(root,"rows"),pattern="\\.rds$",full.names=TRUE)
z <- lapply(fs,readRDS)
stopifnot(length(z)>0)
v <- function(name) vapply(z,`[[`,0,name)
n <- length(z); p <- v("p"); g <- v("greater"); e <- v("equal")
out <- list(n=n,reproduced_p_values=n,standard_rejections=sum(p<=0.05),
  repetitions_with_observed_ties=sum(e>0),max_observed_ties=max(e),
  ignore_permutation_ties_rejections=sum(v("no_perm_ties_p")<=0.05),
  expected_randomized_rejections=sum(v("randomized_reject_probability")),
  possible_randomized_extra_cases=sum(p>0.05 & g/201<=0.05),
  mean_eligible_rejection_ranks=mean(v("orbit_rejections")),
  orbit_rejection_fraction=mean(v("orbit_rejections")/v("nb")),
  full_ten_rejection_ranks=sum(v("orbit_rejections")==10),
  tested_column_equivariance_errors=v("equivariance_error")[is.finite(v("equivariance_error"))])
if(n==500){
  out$unadjusted_binomial_lower_tail <- pbinom(out$standard_rejections,500,10/201)
  out$monte_carlo_interval <- c(qbeta(.025,out$standard_rejections,501-out$standard_rejections),
                               qbeta(.975,out$standard_rejections+1,500-out$standard_rejections))
}
print(out)
saveRDS(out,file.path(root,"summary.rds"))
writeLines(capture.output(print(out)),file.path(root,"summary.txt"))
