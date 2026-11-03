# Produce the task's Chinese pilot report from all fixed-design raw repetitions.
args <- commandArgs(TRUE);stopifnot(length(args)==1L);out <- args[1L]
source("R/summarize.R")
g <- readRDS(file.path(out,"formal_grid.rds"));d <- read_exp("libpilot3",out)
d$err[is.na(d$err)] <- ""
stopifnot(nrow(d)==5100L,all(d$err==""),all(d$B==200L),!anyDuplicated(paste(d$cell,d$rep)))
for(id in g$cell){
  z <- d[d$cell==id,];q <- g[g$cell==id,]
  stopifnot(nrow(z)==q$reps,identical(sort(z$rep),seq_len(q$reps)),
            all(z$seed==q$seed0+z$rep),all(z$lib==q$lib),all(z$libsets==q$libsets))
}
s <- summarize_exp("libpilot3",out)
folder <- file.path(out,"report_tables");dir.create(folder,showWarnings=FALSE)
singles <- paste0("mca_",strsplit(g$lib[1],",",fixed=TRUE)[[1]])
sets <- strsplit(g$libsets[1],";",fixed=TRUE)[[1]]
libs <- paste0("mca_set_",sets);ids <- paste0("L",seq_along(sets)-1L)
sizes <- lengths(strsplit(sets,"+",fixed=TRUE))
opponents <- intersect(c("maxmean","clx_mean","clx_cov","maxmean_ovr","clx_cov_ovr",
  "disco","mmd","mst","schott","schott_perm","mmd_max","kmd","songchen","songchen_A"),names(s))
alt <- s[s$beta>0,];nul <- s[s$beta==0,]
label <- function(x) paste0(x$model," b=",sprintf("%.2f",x$beta)," k=",x$k," d=",x$d,
                           ifelse(x$rot>0,paste0(" rot=",x$rot),""))
power_table <- function(cols,names_out=cols){z <- data.frame(cell=label(alt),alt[,cols,drop=FALSE],check.names=FALSE);names(z)[-1] <- names_out;z}
single_tab <- power_table(singles,sub("^mca_","",singles))
library_tab <- power_table(libs,ids)
opponent_tab <- power_table(opponents)
write.csv(single_tab,file.path(folder,"power_components.csv"),row.names=FALSE)
write.csv(library_tab,file.path(folder,"power_libraries.csv"),row.names=FALSE)
write.csv(opponent_tab,file.path(folder,"power_competitors.csv"),row.names=FALSE)
best_single <- apply(as.matrix(alt[,singles]),1,max,na.rm=TRUE)
selection <- do.call(rbind,lapply(seq_along(libs),function(j){
  gap <- best_single-alt[[libs[j]]]
  data.frame(library=ids[j],components=sets[j],size=sizes[j],max_gap=max(gap),
             cells_gap_over_005=sum(gap>.05+1e-12),mean_power=mean(alt[[libs[j]]]),eligible=all(gap<=.05+1e-12))
}))
write.csv(selection,file.path(folder,"library_rule.csv"),row.names=FALSE)
eligible <- selection[selection$eligible,]
chosen <- if(nrow(eligible)) eligible$library[order(eligible$size,-eligible$mean_power)][1L] else NA_character_
methods <- unique(c("mca",singles,libs,opponents))
level_long <- do.call(rbind,lapply(nul$cell,function(id){
  z <- d[d$cell==id,]
  do.call(rbind,lapply(methods,function(m){
    p <- z[[m]];n <- sum(is.finite(p));x <- sum(p<=.05,na.rm=TRUE)
    ci <- if(n) c(if(x==0)0 else qbeta(.025,x,n-x+1),if(x==n)1 else qbeta(.975,x+1,n-x)) else c(NA,NA)
    data.frame(cell=id,method=m,rejections=x,n=n,rate=if(n)x/n else NA,
      lower95=ci[1],upper95=ci[2],outside_original_band=if(n)x/n<.031||x/n>.069 else NA)
  }))
}))
write.csv(level_long,file.path(folder,"null_rates_intervals.csv"),row.names=FALSE)
watch <- level_long[level_long$method%in%c("mca_xscr1","mca_xscr3","mca_g5f"),]
write.csv(watch,file.path(folder,"monitored_null_flags.csv"),row.names=FALSE)
diag_cols <- grep("^(auc_.+_dev|l1|l1th|l1se|quad|scrm|scr2|knn|rf)_?",names(s),value=TRUE)
diag_cols <- setdiff(diag_cols,c(singles,libs))
diagnostics <- s[,unique(c("cell","model","beta","k","d",diag_cols)),drop=FALSE]
write.csv(diagnostics,file.path(folder,"learner_diagnostics.csv"),row.names=FALSE)
write.csv(s[,c("cell","model","beta","k","d","reps","sec_per_rep")],file.path(folder,"timing.csv"),row.names=FALSE)
best_opponent <- apply(as.matrix(alt[,opponents]),1,max,na.rm=TRUE)
losses <- do.call(rbind,lapply(seq_along(libs),function(j){
  data.frame(library=ids[j],cell=label(alt),power=alt[[libs[j]]],best_competitor=best_opponent,
             gap=best_opponent-alt[[libs[j]]])
}))
losses <- losses[losses$gap>0,];losses <- losses[order(-losses$gap),]
write.csv(losses,file.path(folder,"cells_below_competitors.csv"),row.names=FALSE)
md <- function(z,digits=3L){
  z <- as.data.frame(z,check.names=FALSE)
  for(j in seq_along(z)){
    if(is.numeric(z[[j]]))z[[j]] <- ifelse(is.finite(z[[j]]),formatC(z[[j]],digits=digits,format="f"),"NA")
    z[[j]] <- gsub("|","/",as.character(z[[j]]),fixed=TRUE)
  }
  c(paste0("| ",paste(names(z),collapse=" | ")," |"),paste0("| ",paste(rep("---",ncol(z)),collapse=" | ")," |"),
    apply(z,1,function(r)paste0("| ",paste(r,collapse=" | ")," |")))
}
lines <- c("# 第三轮功效与场景：完整试点报告", "",paste("生成时间：",format(Sys.time())),"",
  "## 范围与完整性","",
  "固定 libpilot3 设计：27 格，21 个备择格各 100 次、6 个原假设格各 500 次，共 5100 次；B=200，alpha=0.05。",
  "全部逐次文件的格号、重复编号、种子、库和库集合均已核对；无重复记录、无记录错误。",
  "原假设阶段曾因 A 格 xscr3=0.028 触发预设停止；用户随后明确授权保留异常并继续完整试点。原停止记录保留，阈值未调整。",
  "统计定义、冻结模型、校准常数及置换引擎未因本轮水平结果而修改。运行采用 R 4.3.3、单进程单线程、关闭 JIT，调节并行度使 CPU 约 70%。",
  "","## 水平结果与偏离记录","",
  "下表保留原先 [0.031,0.069] 的监测范围。置信区间为逐项 Clopper–Pearson 95% 区间，不是所有格和方法的同时置信区间；越界不能直接定位原因。所有单分量、八个库及对照方法的完整水平表见 report_tables/null_rates_intervals.csv。",
  "",md(watch),"","## 备择功效：单分量","",md(single_tab),
  "","## 备择功效：预先固定的八个库","",md(library_tab),"",md(selection),"",
  "上述资格仅按预设功效差距规则计算；水平异常仍需单列复核。用户授权继续计算不等于水平验证已经通过。",
  if(is.na(chosen))"**没有候选库满足所有备择格距最佳单分量不超过 0.05 的预设条件，因此不能按该规则给出合格推荐。未放宽规则，未冻结新库。**"
  else paste0("按预设规则得到的建议为 **",chosen,"**；这是试点建议，最终库仍由用户决定，尚未冻结。"),
  "","## 对照方法及不利结果","",md(opponent_tab),"",
  "所有候选库低于该格最佳对照方法的格子均保存在 report_tables/cells_below_competitors.csv。以下列出经验差距最大的最多 20 项；每格仅 100 次，不把经验大小关系直接表述为显著优劣。",
  "",if(nrow(losses))md(head(losses,20)) else "这批试点中未出现所列库低于最佳对照的格子。",
  "","## AUC、调参诊断与耗时","",
  "完整数据见 report_tables/learner_diagnostics.csv、timing.csv 和 libpilot3_timing.csv。",
  "已有 AUC 列仅对应总体 k 的测试 AUC；S2 对应总体 4，不能代表总体 2、3。xscr 分量原输出未保存所选特征编号及双向 AUC，不能凭这些结果断言变量筛选命中或失败。",
  "l1th 在部分原假设格可能输出空模型；应与其非零系数记录共同解读，不把零拒绝率隐藏。",
  "筛选失误、库临界值代价和不同偏离结构均可能影响功效；现有诊断不足以给每个输掉的格子确认因果机制。",
  "",md(s[,c("model","beta","k","d","reps","sec_per_rep")]),
  "","## 后续决定","",
  "到达任务书停止点 2：由用户选择库及正文/补充材料场景，再冻结 LIB_MAIN_V3。E1–E13 的快速版与精确版尚未启动。",
  "本轮未实现可选 ksrf/deprf，不能把它们列为已经比较的方法。",
  "","## 附录：第三轮既有 E9/E10 表",""
)
appendix <- file.path("results_v2_quick","review3_tables.txt")
if(file.exists(appendix)) {
  lines <- c(lines,"```text",readLines(appendix,warn=FALSE),"```")
} else lines <- c(lines,"review3_tables.txt 未生成；需查看运行日志，不能声称该附录已完成。")
report <- "试点报告_第三轮_功效与场景.md"
if(file.exists(report))file.copy(report,paste0(report,".before_",format(Sys.time(),"%Y%m%d_%H%M%S")),overwrite=FALSE)
writeLines(lines,report,useBytes=TRUE)
writeLines(lines,file.path(out,"试点报告_第三轮_功效与场景.md"),useBytes=TRUE)
cat("COMPLETE PILOT REPORT:",normalizePath(report),"\n")
