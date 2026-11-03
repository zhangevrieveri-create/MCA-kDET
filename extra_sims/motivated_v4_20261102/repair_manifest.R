out <- "motivated_v4_20261102/smoke_B19"
original <- read.csv(file.path(out,"source_md5.csv"))
actual <- as.character(tools::md5sum(original$file))
changed <- original$file[actual!=original$md5]
stopifnot(identical(changed,"motivated_v4_20261102/smoke.R"))
stopifnot(as.character(tools::md5sum("motivated_v4_20261102/backup_20261102_0956/smoke.R"))==original$md5[original$file==changed])
original$md5 <- actual
stopifnot(!file.exists(file.path(out,"source_md5_revision1.csv")))
write.csv(original,file.path(out,"source_md5_revision1.csv"),row.names=FALSE)
writeLines(c("2026-11-02: smoke runner p-value assertion mistakenly included kmd_eta, a signed effect statistic returned inside the legacy competitor $p list.",
"Fix only the assertion and p-value export: exclude the explicitly named kmd_eta; preserve every raw record and all initial error logs.",
"All original simulation, learner and comparator sources verified unchanged. Backup reproduces the original smoke runner hash.",
"The initial manifest is preserved; revision1 is an explicit audited source update. No data, model, seed, statistic, library or calibration change.",
"Failed cells resume identical seeds, up to three process attempts after this diagnosed wrapper repair."),file.path(out,"repair_record.txt"))
