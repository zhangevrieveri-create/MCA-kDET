// Research workbook: verified financial results, without credentials or API calls.
import fs from "node:fs/promises";
import path from "node:path";
import os from "node:os";
import { createRequire } from "node:module";
import { fileURLToPath, pathToFileURL } from "node:url";
const root = path.dirname(fileURLToPath(import.meta.url));
const runtime = path.join(os.tmpdir(), "mca-singlecell-artifact-runtime");
await fs.mkdir(runtime, { recursive: true });
try { await fs.symlink("C:/Users/zhang/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules", path.join(runtime, "node_modules"), "junction"); }
catch (error) { if (error.code !== "EEXIST") throw error; }
const require = createRequire(path.join(runtime, "package.json"));
const { Workbook, SpreadsheetFile } = await import(pathToFileURL(require.resolve("@oai/artifact-tool")).href);
const input = JSON.parse(await fs.readFile(path.join(root, "results/workbook_data.json"), "utf8"));
const plan = JSON.parse(await fs.readFile(path.join(root, "study_plan.json"), "utf8"));
const objects = data => data.data.map(row => Object.fromEntries(data.columns.map((key, i) => [key, row[i]])));
const wb = Workbook.create();
const names = { MCA: "MCA-kDET", MCA_l1: "L1", MCA_quad: "Quad", MCA_knn: "kNN", MCA_rf: "RF", Two_sided: "Two-sided", Sum: "Sum", Comb: "Comb", DISCO: "DISCO", MMD: "MMD", MST: "MST" };
const mainMethods = ["MCA", "DISCO", "MMD", "MST"];
const configNames = { primary: "主分析10股", dimension_20: "20股", dimension_30: "30股", block_13: "13周分块", block_52: "52周分块", buffer_8: "8周缓冲", raw_blocks: "未过滤分块", filtered_random: "过滤后随机拆分", vix_15_25: "VIX阈值15/25" };
for (let i = 0; i < 4; i++) configNames[`thin4_offset${i}`] = `每4周测试偏移${i}`;
for (let i = 1; i <= 5; i++) configNames[`split_${String(i).padStart(2, "0")}`] = `额外拆分${i}`;
function sheet(name, title) {
  const ws = wb.worksheets.add(name); ws.showGridLines = false;
  ws.getRange("A1:K40").format.font = { name: "Arial", size: 10, color: "#202A35" };
  ws.getRange("A2").values = [[title]]; ws.getRange("A2").format.font = { name: "Arial", size: 14, bold: true };
  ws.getRange("A3:J3").format.borders = { bottom: { style: "thin", color: "#667788" } };
  return ws;
}
function table(ws, headers, rows) {
  const range = ws.getRangeByIndexes(5, 0, rows.length + 1, headers.length);
  range.values = [headers, ...rows];
  range.format = { columnWidth: 16, rowHeight: 21, verticalAlignment: "center", font: { name: "Arial", size: 10, color: "#202A35" } };
  ws.getRangeByIndexes(5, 0, 1, headers.length).format = { fill: "#243B53", rowHeight: 25, horizontalAlignment: "center", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" } };
  for (let i = 1; i < rows.length; i += 2) ws.getRangeByIndexes(i + 6, 0, 1, headers.length).format.fill = "#F2F5F7";
}
const main = sheet("主结果", "A股周收益：外部VIX状态的分布检验");
main.tabColor = "#243B53";
main.getRange("A4").values = [["10股；GARCH过滤；26周分块；4周缓冲；4999次标签置换，近似校准"]];
const pvalues = objects(input.main);
const order = [...mainMethods, ...pvalues.map(r => r.method).filter(m => !mainMethods.includes(m))];
table(main, ["方法", "p值", "蒙特卡洛标准误", "结果用途"], order.map(m => {
  const r = pvalues.find(row => row.method === m);
  return [names[m], r.p_value, r.MC_SE, mainMethods.includes(m) ? "四方法主比较" : "补充诊断"];
}));
main.getRange("B7:C17").setNumberFormat("0.0000");
main.getRange("C6:C17").format.columnWidth = 24;
main.getRange("D6:D17").format.columnWidth = 20;
main.getRange("F6:G12").values = [
  ["训练 / 测试", "212 / 211周"], ["有效日期", "2015-01-09至2025-01-31"],
  ["校准", "依赖数据下的近似周标签置换"], ["主结论", "四方法均未检出显著差异"],
  ["敏感性", "20/30股时部分对比方法显著，MCA未显著"],
  ["数据缺失", "573个目标日历周中121周无完整股价向量"],
  ["股票名单", "只用2014年底权重及此前数据质量固定"],
];
main.getRange("F6:F12").format.columnWidth = 18;
main.getRange("G6:G12").format.columnWidth = 63;
main.getRange("A20").values = [["未拒绝不等于证明同分布；敏感性结果未作多重校正，不能替换预定主结果。"]];
const sensitivity = sheet("全部分析", "全部18项预定分析");
const comparisons = objects(input.comparison), checks = objects(input.checks);
table(sensitivity, ["设定", "股票数", "训练周", "测试周", "MCA-kDET", "DISCO", "MMD", "MST"], plan.specifications.map(s => {
  const r = comparisons.find(row => row.configuration === s.name), c = checks.find(row => row.configuration === s.name);
  return [configNames[s.name], s.dimension, c.fitting_weeks, c.testing_weeks, ...mainMethods.map(m => r[m])];
}));
sensitivity.getRange("A6:A24").format.columnWidth = 30;
sensitivity.getRange("E7:H24").setNumberFormat("0.0000");
sensitivity.getRange("E7:H24").conditionalFormats.add("cellIs", { operator: "lessThanOrEqual", formula: .05, format: { fill: "#E1EEE9", font: { color: "#155E49" } } });
sensitivity.getRange("A4").values = [["加股票会改变可用周交集；反复拆分共享日期，不能解释为独立功效重复。"]];
sensitivity.freezePanes.freezeRows(6);
const stocks = sheet("股票及过滤", "历史股票名单与预样本GARCH参数");
table(stocks, ["排名", "股票代码", "名称", "2014年底权重", "样本均值", "omega", "alpha", "beta", "持续性", "收敛码", "基准期有效周"], objects(input.stocks).map(r =>
  [r.selected_rank, r.order_book_id, r.symbol, r.weight_20141231, r.mean, r.omega, r.alpha, r.beta, r.persistence, r.convergence, r.baseline_observations]));
stocks.getRange("B6:B36").format.columnWidth = 22;
stocks.getRange("D7:D36").setNumberFormat("0.0000%");
stocks.getRange("E7:I36").setNumberFormat("0.000000");
stocks.getRange("A4").values = [["参数仅用2010–2014周收益拟合；均值和omega的单位基于100×复权对数收益。"]];
stocks.freezePanes.freezeRows(6);
const sample = sheet("样本与说明", "每状态的样本删失与分析解释");
table(sample, ["VIX组", "日历周", "缺失周", "训练周", "测试周", "缓冲周"], objects(input.attrition).map(r =>
  [["低", "中", "高"][r.group_number - 1], r.calendar_weeks, r.missing, r.fitting, r.testing, r.buffer]));
const notes = [
  ["状态", "上一周美国VIX；基准期阈值14.876667、18.43，未用当周A股收益分组。"],
  ["日期", "计划2015–2025，历史股票后续停牌/退市使有效完整向量止于2025-01-31。"],
  ["收益", "后复权周对数收益×100；每周中国最后交易日，两个端点必须可观测且未全天停牌。"],
  ["过滤", "逐股GARCH(1,1)参数仅用2015年前估计，向前递推仅用过去收益。"],
  ["缺失", "保留缺失日历周，不前填价格，不把多周复牌收益当作一周。"],
  ["检验", "4999随机标签加1个观测分配；依赖数据下近似校准，未保证有限样本精确水平。"],
  ["限制", "金融业股票集中、VIX是美国风险、缺失随状态不同、残差可能仍有依赖。"],
  ["核验", "30条GARCH递推及18项分块、学习器、置换和对比统计量经独立复算。"],
  ["米筐文档", "https://www.ricequant.com/doc/rqdata/python/generic-api"],
  ["VIX来源", "https://www.cboe.com/tradable-products/vix/vix-historical-data/"],
];
sample.getRange("H6:I15").values = notes;
sample.getRange("H6:H15").format.columnWidth = 18;
sample.getRange("I6:I15").format.columnWidth = 115;
wb.recalculate();
console.log((await wb.inspect({ kind: "table", range: "主结果!A6:D17", tableMaxRows: 12, tableMaxCols: 4, maxChars: 2000 })).ndjson);
const preview = await wb.render({ sheetName: "主结果", range: "A2:G17", scale: 1.5, format: "png" });
const previewDir = process.argv[2] || runtime;
await fs.mkdir(previewDir, { recursive: true });
await fs.writeFile(path.join(previewDir, "financial_workbook.png"), new Uint8Array(await preview.arrayBuffer()));
await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(root, "results/股票实证结果.xlsx"));
console.log("Exported financial research workbook.");
