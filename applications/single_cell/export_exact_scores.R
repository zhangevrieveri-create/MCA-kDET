# Preserve all binary64 score distinctions for independent rank verification.
# The original analysis.rds and original CSV are kept unchanged.
for (folder in commandArgs(TRUE)) {
  result <- readRDS(file.path(folder, "analysis.rds"))
  tab <- result$test_metadata[, c("sampleID", "verified_animal_block", "group", "group_number", "geo_tumor", "geo_Sex")]
  for (a in result$sc$lib) for (i in 1:3)
    tab[[paste0(a, "_target", i)]] <- sprintf("%.17g", result$sc$Sovr[[a]][, i])
  write.csv(tab, file.path(folder, "test_scores_exact.csv"), row.names = FALSE)
}
