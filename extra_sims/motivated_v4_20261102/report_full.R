source("motivated_v4_20261102/full_common.R")
out <- motiv_arg("out",file.path(MOTIV_DIR,"full_B200"));motiv_manifest_check(out)
grid <- readRDS(file.path(out,"grid.rds"));smoke <- all(grid$B==19L)
tabdir <- file.path(out,"report_tables");dir.create(tabdir,showWarnings=FALSE)
files <- list.files(file.path(out,"rows"),pattern="[.]rds$",full.names=TRUE)
stopifnot(length(files)==sum(grid$reps));rows <- lapply(files,readRDS)
keys <- vapply(rows,function(z)paste(z$cell,z$rep),character(1))
expected <- unlist(lapply(seq_len(nrow(grid)),function(j)paste(grid$cell[j],seq_len(grid$reps[j]))))
stopifnot(!anyDuplicated(keys),setequal(keys,expected))
integrity <- list()
for(j in seq_len(nrow(grid))) {
  g <- as.list(grid[j,]);rr <- rows[vapply(rows,function(z)z$cell==g$cell,logical(1))]
  method_names <- sort(names(rr[[1L]]$p))
  for(z in rr) {
    stopifnot(z$seed==g$seed0+z$rep,z$B==g$B,z$err=="",z$model==g$model,z$target==g$target,
      z$k==g$k,z$d==g$d,z$tau==g$tau,z$lib==g$lib,z$libsets==g$libsets,
      identical(sort(names(z$p)),method_names),!anyDuplicated(names(z$p)),!("kmd_eta"%in%names(z$p)),
      all(is.finite(unlist(z$p))),all(unlist(z$p)>=0&unlist(z$p)<=1))
  }
  integrity[[j]] <- data.frame(cell=g$cell,expected=g$reps,saved=length(rr),B=g$B,
                              unique_keys=TRUE,seeds=TRUE,libraries=TRUE,record_errors=0L)
}
write.csv(do.call(rbind,integrity),file.path(tabdir,"integrity.csv"),row.names=FALSE)
pvalues <- do.call(rbind,lapply(rows,function(z)data.frame(cell=z$cell,rep=z$rep,seed=z$seed,
                      method=names(z$p),p=unlist(z$p),check.names=FALSE)))
aucs <- do.call(rbind,lapply(rows,function(z)cbind(cell=z$cell,rep=z$rep,z$aucs)))
xscr <- do.call(rbind,lapply(rows,function(z)cbind(cell=z$cell,rep=z$rep,z$xscr)))
timing <- do.call(rbind,lapply(rows,function(z)data.frame(cell=z$cell,rep=z$rep,seconds=z$seconds,retries=z$retries)))
write.csv(pvalues,file.path(tabdir,"pvalues.csv"),row.names=FALSE)
write.csv(aucs,file.path(tabdir,"learner_auc_per_population.csv"),row.names=FALSE)
write.csv(xscr,file.path(tabdir,"xscr_selected_and_auc.csv"),row.names=FALSE)
write.csv(timing,file.path(tabdir,"timing_per_rep.csv"),row.names=FALSE)
ci <- function(y,n) c(if(y==0)0 else qbeta(.025,y,n-y+1),if(y==n)1 else qbeta(.975,y+1,n-y))
splitp <- split(pvalues,paste(pvalues$cell,pvalues$method,sep="|"))
rates <- do.call(rbind,lapply(splitp,function(z) {
  count <- sum(z$p<=.05);n <- nrow(z);interval <- ci(count,n)
  data.frame(cell=z$cell[1L],method=z$method[1L],rejections=count,reps=n,rate=count/n,
             lower95=interval[1L],upper95=interval[2L],stringsAsFactors=FALSE)
}))
meta <- grid[,c("cell","model","target","k","n","d","B")]
rates <- merge(meta,rates,by="cell",sort=FALSE)
rates <- rates[order(match(rates$cell,grid$cell),rates$method),]
write.csv(rates,file.path(tabdir,"all_method_rates_and_intervals.csv"),row.names=FALSE)
null <- rates[rates$target<=.5,,drop=FALSE]
null$ci_excludes_005 <- null$upper95<.05|null$lower95>.05
null$outside_003_007 <- null$rate<.03|null$rate>.07
write.csv(null,file.path(tabdir,"null_rates_and_flags.csv"),row.names=FALSE)
auc_summary <- aggregate(auc~cell+learner+population,aucs,function(x)c(mean=mean(x),sd=sd(x),n=length(x)))
auc_summary <- cbind(auc_summary[,1:3],as.data.frame(auc_summary$auc))
write.csv(auc_summary,file.path(tabdir,"learner_auc_summary.csv"),row.names=FALSE)
timing_summary <- aggregate(seconds~cell,timing,function(x)c(mean=mean(x),median=median(x),total=sum(x)))
timing_summary <- cbind(timing_summary[1L],as.data.frame(timing_summary$seconds))
write.csv(timing_summary,file.path(tabdir,"timing_summary.csv"),row.names=FALSE)

sets <- strsplit(LIBPILOT3_SETS,";",fixed=TRUE)[[1L]];libcols <- paste0("mca_set_",sets)
components <- paste0("mca_",parse_lib(LIBPILOT3_LIB))
alternatives <- grid$cell[grid$target>.5]
rate_at <- function(cell,method) {
  z <- rates$rate[rates$cell==cell&rates$method==method]
  if(length(z)!=1L)stop("Missing or repeated method ",method," ",cell)
  z
}
library_rule <- do.call(rbind,lapply(seq_along(libcols),function(j) {
  perf <- vapply(alternatives,rate_at,numeric(1),method=libcols[j])
  best <- vapply(alternatives,function(c)max(vapply(components,function(m)rate_at(c,m),numeric(1))),numeric(1))
  gap <- best-perf
  data.frame(library=paste0("L",j-1L),components=sets[j],size=length(strsplit(sets[j],"+",fixed=TRUE)[[1L]]),
    mean_power=mean(perf),worst_single_component_deficit=max(gap),failing_cells=sum(gap>.05+1e-12),
    eligible=all(gap<=.05+1e-12))
}))
write.csv(library_rule,file.path(tabdir,"fixed_library_rule.csv"),row.names=FALSE)
comparators <- c("maxmean","clx_mean","clx_cov","disco","mmd","mst","schott","schott_perm",
                 "mmd_max","kmd","songchen","songchen_A","maxmean_ovr","clx_cov_ovr")
paired <- list()
for(cell in alternatives)for(lib in c(libcols[1L],libcols[7L]))for(comp in comparators) {
  a <- pvalues[pvalues$cell==cell&pvalues$method==lib,c("rep","p")]
  b <- pvalues[pvalues$cell==cell&pvalues$method==comp,c("rep","p")]
  if(!nrow(b))next
  z <- merge(a,b,by="rep");stopifnot(nrow(z)==nrow(a))
  diff <- as.numeric(z$p.x<=.05)-as.numeric(z$p.y<=.05)
  estimate <- mean(diff);se <- sd(diff)/sqrt(length(diff))
  paired[[length(paired)+1L]] <- data.frame(cell=cell,library=if(lib==libcols[1L])"L0" else "L6",comparator=comp,
    reps=length(diff),power_difference=estimate,paired_se=se,normal_lower95=max(-1,estimate-1.96*se),
    normal_upper95=min(1,estimate+1.96*se),library_only=sum(z$p.x<=.05&z$p.y>.05),comparator_only=sum(z$p.x>.05&z$p.y<=.05))
}
write.csv(do.call(rbind,paired),file.path(tabdir,"paired_power_differences.csv"),row.names=FALSE)

display_methods <- c(libcols[1L],libcols[7L],"mca_rf","maxmean","clx_cov","mmd_max","kmd","songchen","songchen_A")
display_names <- c("L0","L6","RF分量","均值max","协方差max","MMD-max","KMD","SC-S","SC-SA")
comparison <- meta
for(j in seq_along(display_methods))comparison[[display_names[j]]] <- vapply(comparison$cell,rate_at,numeric(1),method=display_methods[j])
write.csv(comparison,file.path(tabdir,"L0_L6_and_competitors.csv"),row.names=FALSE)
motiv_manifest_check(out)
nullflags <- null[null$ci_excludes_005| null$outside_003_007,,drop=FALSE]
error_logs <- list.files(file.path(out,"errors"),pattern="[.]log$",full.names=TRUE)
smoke_exact <- sum(vapply(rows,function(z)isTRUE(z$same_seed_exact),logical(1)))
verification <- list(status="verified",cells=nrow(grid),replications=length(rows),expected=sum(grid$reps),
  unique_cell_rep=TRUE,all_seeds_correct=TRUE,all_B_correct=TRUE,all_libraries_correct=TRUE,
  saved_record_errors=0L,retained_calculation_error_logs=length(error_logs),
  null_rate_flags=nrow(nullflags),source_fingerprints_unchanged=TRUE,
  same_seed_exact_checks=smoke_exact,smoke_only=smoke)
saveRDS(verification,file.path(out,"verification.rds"))
writeLines(c(paste("CELLS",nrow(grid)),paste("REPLICATIONS",length(rows)),paste("EXPECTED",sum(grid$reps)),
  "Unique (cell,rep), every seed=seed0+rep, B, library, parameters and p-value checks passed.",
  "All source fingerprints unchanged; all observed xscr diagnostics agree with the production statistic.",
  paste("Saved record errors:",0,"; retained calculation error logs:",length(error_logs)),
  paste("Null method/cell flags:",nrow(nullflags)),paste("Smoke same-seed exact checks:",smoke_exact)),
  file.path(out,"verification.txt"))
report <- c("# 新场景全方法比较报告", "",paste("生成时间：",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),"",
  "## Material Passport", "",
  "任务：用户授权补充适合现有分类器的机制场景，并比较全部现有方法。状态：原始结果、完整性与汇总核验完成。来源：本地固定网格模拟；未使用外部数据。",
  "",
  if(smoke)"**这是机械小试：每格仅3次、B=19。下面的拒绝比例不能用于可靠功效排序，不据此删格或调参。**" else
    "这是探索性功效试点：备择每格100次，原假设每格500次，B=200，拒绝阈值0.05。备择100次的误差仍较大，接近的比例不能直接解释为可靠优胜。",
  "",
  "## 网格与核验", "",
  sprintf("%d格、%d次重复全部保存；(cell,rep)唯一，种子、B、原模型参数、库和源码指纹核验通过，已保存结果中的错误字段均为空。",nrow(grid),length(rows)),
  sprintf("重试诊断日志保留%d份；详见 errors/ 与 worker_logs/，不因重试而更换种子。",length(error_logs)),
  "",
  "原四类 N1–N4 全部保留，并加入 N5_skewmix、N6_ring。主维数20、强度理想AUC 0.60/0.70/0.80；维数100保留AUC0.70。各组n=300，拟合/测试各半。四个原假设格覆盖k3/k4与d20/d100。",
  "",
  "校准使用真实密度比对基准正态的理想AUC，与任何检验拒绝率无关。混合模型各偏离分别对同一个基准校准，标签不是一对其余混合分布的理想AUC。公式和校准详见同目录设计文件与 calibration_check.csv、calibration_extra_check.csv（如该独立复核文件采用其他命名，以校准脚本输出为准）。",
  "",
  "## 全格对比", "",
  "L0为旧库 l1+quad+knn+rf；L6为 xscr1+xscr3+knn+rf+g5f。RF分量也对全部总体取最大值。这里始终展示固定的L0/L6，没有逐格挑库。SC-S/SC-SA为Song–Chen的两个统计量。",
  "",
  paste0("| 模型 | 理想AUC | d | k | ",paste(display_names,collapse=" | ")," |"),
  paste0("|---|---:|---:|---:|",paste(rep("---:",length(display_names)),collapse="|"),"|"))
for(i in seq_len(nrow(comparison))) {
  z <- comparison[i,]
  report <- c(report,paste0("| ",z$model," | ",sprintf("%.2f",z$target)," | ",z$d," | ",z$k," | ",
    paste(sprintf("%.3f",as.numeric(unlist(z[,display_names],use.names=FALSE))),collapse=" | ")," |"))
}
report <- c(report,"", "完整方法（包括全部12个分量、8个固定库、均值与协方差最大值及其适用的一对其余版本、DISCO、MMD、MST、Schott、置换Schott、最大两两MMD、KMD和Song–Chen）见 report_tables/all_method_rates_and_intervals.csv，每项附Clopper–Pearson 95%区间。",
  "",
  "L0/L6与每个对照的同次数据拒绝差见 paired_power_differences.csv。该表正态区间仅作试点诊断；没有对大量模型和方法比较作确认性多重校正。",
  "",
  "## 固定选库规则", "",
  "按此前固定规则：每个备择格都不能比最好的单分量低超过0.05，然后优先分量少、再看平均功效。本批也输出同一规则的结果，不因新模型而放宽。",
  "",
  if(any(library_rule$eligible))paste("本批网格满足该规则的候选：",paste(library_rule$library[library_rule$eligible],collapse=", "),"。这不自动推翻上一轮全网格的失败，也不自动冻结最终库。") else "本批网格无合格候选；不临时放宽阈值。",
  "",
  "全部L0–L7结果和失败格数见 fixed_library_rule.csv。LIB_MAIN_V3仍由作者决定，程序没有冻结最终库，也未启动E1–E13最终库模拟。",
  "",
  "## 原假设水平", "",
  sprintf("按95%%区间不含0.05或比例超出[0.03,0.07]的描述性规则，共标记%d个方法×原假设格。全部记录保留；这不是大量诊断后已校正的显著性结论。",nrow(nullflags)),
  "",
  "全表见 null_rates_and_flags.csv。不能通过改阈值、增加样本或换种子消除告警，也不能把不拒绝原假设解释为总体已证实相同。",
  "",
  "## 学习与筛选诊断", "",
  "learner_auc_per_population.csv及其summary逐总体记录留出AUC，不只记最后一个总体；xscr_selected_and_auc.csv保存两方向所选变量及AUC，每次均与原xscr观察统计量重建核对。这样可以区分学习失败与最大值校准代价，不能仅凭拒绝率猜测筛选机制。",
  "",
  "## 探索边界与后续", "",
  "这些场景是应用户要求、在已看过原四类3次小试后补充的探索设计。新网格与强度在新增正式结果前固定，全部输赢保留；它们不能被包装成独立于探索过程的预注册普遍优越性证据。",
  "",
  "总体均值和协方差相同，所以超过均值/协方差专门方法本身不能证明优于一般分布检验；应同时查看核、距离与图方法。N3同时改变边际形状和依赖，不能作为纯copula差异的证明。不同模型的学习难度也不能仅由理想AUC相同视为完全相同。",
  "",
  "原A–E、S1–S3及不利结果未删除；上一轮完整结果见 ../results_v3_pilot_20261101_cpu16/，中文报告为 sim/试点报告_第三轮_功效与场景.md。上一轮A格xscr3=0.028、D格g5f=0.024水平偏离继续保留，不因本轮模型改变而抹去。",
  "",
  "尚未修改论文、替作者决定正文/补充场景或冻结最终库。下一步先共同审阅各格差异、水平与置信区间，再确定需要独立确认的科研设定。")
writeLines(report,file.path(out,"模拟报告_新场景全方法.md"),useBytes=TRUE)
if(!smoke)writeLines(report,file.path(MOTIV_DIR,"全方法比较报告.md"),useBytes=TRUE)
cat("REPORT COMPLETE;",nrow(grid),"CELLS;",length(rows),"REPLICATIONS VERIFIED\n")
