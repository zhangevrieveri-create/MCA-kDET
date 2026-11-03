# Run from sim/. Execute only deterministic sections 1-4 of test_power3.R.
# This deliberately cannot report POWER3 TESTS OK or replace the full test.
ee <- parse("tests/test_power3.R")
for (j in seq.int(3L, length(ee))) {
  e <- ee[[j]]
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name("lib5"))) break
  eval(e, envir = .GlobalEnv)
}
stopifnot(all_ok)
cat("SMOKE ONLY OK; full validation not run\n")
