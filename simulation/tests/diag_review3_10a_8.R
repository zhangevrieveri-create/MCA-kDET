# Diagnostic for test_review3.R checks 8 and 10a (writes nothing)
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", fa[1L])))
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R",
            "competitors2.R", "rotation.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
source(file.path(SRC_DIR, "summarize.R"))
suppressPackageStartupMessages(requireNamespace("glmnet"))
suppressPackageStartupMessages(requireNamespace("ranger", quietly = TRUE))

par0 <- list(cA = 0.6, cB = 4, sC = 2L, cE = 0.2)
cat("R", as.character(getRversion()), " glmnet", as.character(packageVersion("glmnet")), "\n")

# ---- 8: Song-Chen raw output
set.seed(11)
dd <- gen_data("A", 1, c(40L, 40L, 40L), 10L, modifyList(par0, list(cA = 3)))
D2 <- sqdist(do.call(rbind, dd))
E <- kmst_edges(sqrt(D2), SC_KMST)
r <- gTestsMulti::gtestsmulti(E, dd, perm = 99L)
str(r$pval)
set.seed(12)
dd0 <- gen_data("A", 0, c(40L, 40L, 40L), 10L, modifyList(par0, list(cA = 3)))
E0 <- kmst_edges(sqrt(sqdist(do.call(rbind, dd0))), SC_KMST)
str(gTestsMulti::gtestsmulti(E0, dd0, perm = 99L)$pval)
