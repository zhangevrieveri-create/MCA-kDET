// Build the research workbook from verified result tables, using bundled Artifact Tool.
import fs from "node:fs/promises";
import path from "node:path";
import os from "node:os";
import { createRequire } from "node:module";
import { fileURLToPath, pathToFileURL } from "node:url";

const root = path.dirname(fileURLToPath(import.meta.url));
const modules = "C:/Users/zhang/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules";
const runtime = path.join(os.tmpdir(), "mca-singlecell-artifact-runtime");
await fs.mkdir(runtime, { recursive: true });
try { await fs.symlink(modules, path.join(runtime, "node_modules"), "junction"); }
catch (error) { if (error.code !== "EEXIST") throw error; }
const require = createRequire(path.join(runtime, "package.json"));
const { Workbook, SpreadsheetFile } = await import(pathToFileURL(require.resolve("@oai/artifact-tool")).href);
const data = JSON.parse(await fs.readFile(path.join(root, "results/workbook_data.json"), "utf8"));
const plan = JSON.parse(await fs.readFile(path.join(root, "study_plan.json"), "utf8"));
const workbook = Workbook.create();
const font = "Arial";
const tableObjects = (table) => table.data.map((row) => Object.fromEntries(table.columns.map((column, i) => [column, row[i]])));
const short = { MCA: "MCA-kDET", MCA_l1: "L1", MCA_quad: "Quad", MCA_knn: "kNN", MCA_rf: "RF", Two_sided: "Two-sided", Sum: "Sum", Comb: "Comb", Pair: "Pair", DISCO: "DISCO", MMD: "MMD", MST: "MST" };
const primaryNames = ["MCA", "DISCO", "MMD", "MST"];
const label = { primary: "主分析", genes_50: "50基因", genes_200: "200基因", cells_20: "每鼠20细胞", cells_40: "每鼠40细胞", middle_20w_only: "中间组仅20周", one_tumor: "每只晚期鼠一个肿瘤" };

function sheet(name, title) {
  const ws = workbook.worksheets.add(name);
  ws.showGridLines = false;
  ws.getRange("A1:L40").format.font = { name: font, size: 10, color: "#202A35" };
  ws.getRange("A2").values = [[title]];
  ws.getRange("A2").format.font = { name: font, size: 14, bold: true };
  ws.getRange("A3:J3").format.borders = { bottom: { style: "thin", color: "#667788" } };
  return ws;
}
function table(ws, headerRow, headers, rows) {
  ws.getRangeByIndexes(headerRow - 1, 0, rows.length + 1, headers.length).values = [headers, ...rows];
  const full = ws.getRangeByIndexes(headerRow - 1, 0, rows.length + 1, headers.length);
  full.format.font = { name: font, size: 10, color: "#202A35" };
  full.format.rowHeight = 20;
  full.format.verticalAlignment = "center";
  full.format.columnWidth = 15;
  const head = ws.getRangeByIndexes(headerRow - 1, 0, 1, headers.length);
  head.format.fill = "#243B53";
  head.format.font = { name: font, size: 10, color: "#FFFFFF", bold: true };
  head.format.horizontalAlignment = "center";
  head.format.rowHeight = 25;
  for (let i = 0; i < rows.length; i++) {
    if (i % 2 === 1) ws.getRangeByIndexes(headerRow + i, 0, 1, headers.length).format.fill = "#F2F5F7";
  }
  ws.getRangeByIndexes(headerRow + rows.length - 1, 0, 1, headers.length).format.borders = { bottom: { style: "thin", color: "#BAC5CF" } };
}
const main = sheet("主结果", "KP肺癌阶段：单细胞分布检验");
main.tabColor = "#243B53";
main.getRange("A4").values = [["12 / 18–20 / 30周；100基因；每鼠30细胞；完整整鼠标签枚举"]];
const pvals = tableObjects(data.main);
const order = [...primaryNames, ...pvals.map(r => r.method).filter(m => !primaryNames.includes(m))];
const mainRows = order.map(method => {
  const r = pvals.find(row => row.method === method);
  return [short[method], Math.round(r.p_value * r.assignments), r.assignments, null,
    primaryNames.includes(method) ? "四方法主比较" : "补充诊断"];
});
table(main, 6, ["方法", "同样或更极端的分配", "全部分配数", "p值", "结果用途"], mainRows);
main.getRange("B6:B18").format.columnWidth = 24;
main.getRange("A6:A18").format.columnWidth = 18;
main.getRange("E6:E18").format.columnWidth = 20;
main.getRange("B7:D18").format.horizontalAlignment = "right";
main.getRange("B7:C18").setNumberFormat("0");
main.getRange("D7:D18").formulas = mainRows.map((_, i) => [`=B${i + 7}/C${i + 7}`]);
main.getRange("D7:D18").setNumberFormat("0.000000");
main.getRange("D7:D18").conditionalFormats.add("cellIs", { operator: "lessThanOrEqual", formula: 0.05,
  format: { fill: "#E1EEE9", font: { color: "#155E49" } } });
const notes = [
  ["样本", "14只小鼠；2583个可用细胞"],
  ["训练 / 测试", "6鼠180细胞 / 8鼠240细胞"],
  ["原假设", "抽样动物细胞块联合分布可交换"],
  ["水平", "0.05；颜色仅标示p≤0.05"],
  ["补充检验", "诊断用途；不按最小p选择主方法"],
  ["解释", "横断面阶段关联；存在取材和来源混杂"],
  ["稳定性", "10次额外拆分中MCA有3次p>0.05"],
  ["核验", "17项分析均经独立复算"],
];
main.getRange("G6:H13").values = notes;
main.getRange("G6:G13").format = { columnWidth: 17, font: { name: font, size: 10, bold: true } };
main.getRange("H6:H13").format.columnWidth = 53;
main.getRange("A21").values = [["Comb的分配计数来自联合排名的下尾；其他方法为统计量上尾。"]];

const sensitivity = sheet("稳健性", "完整敏感性与额外拆分结果");
sensitivity.getRange("A4").values = [["每个设定均保留。重复拆分共享动物，次数占比不能解释为功效。"]];
const comparisons = tableObjects(data.comparison);
const allResults = tableObjects(data.all_results);
const sensitivityRows = plan.specifications.map(spec => {
  const c = comparisons.find(r => r.configuration === spec.name);
  const raw = allResults.find(r => r.configuration === spec.name);
  return [label[spec.name] || spec.name.replace("split_", "额外拆分"), spec.seed, spec.dimension,
    spec.cells, raw.test_animals, raw.assignments, ...primaryNames.map(m => c[m])];
});
table(sensitivity, 6, ["设定", "随机种子", "基因数", "每鼠细胞", "测试小鼠", "标签分配", "MCA-kDET", "DISCO", "MMD", "MST"], sensitivityRows);
sensitivity.getRange("A6:A23").format.columnWidth = 28;
sensitivity.getRange("B7:F23").setNumberFormat("0");
sensitivity.getRange("G7:J23").setNumberFormat("0.000000");
sensitivity.getRange("G7:J23").conditionalFormats.add("cellIs", { operator: "lessThanOrEqual", formula: .05,
  format: { fill: "#E1EEE9", font: { color: "#155E49" } } });
sensitivity.freezePanes.freezeRows(6);

const animals = sheet("小鼠样本", "独立动物及原始采样结构");
const animalsData = tableObjects(data.animals);
table(animals, 6, ["小鼠编号", "分组", "可用细胞", "真实周数", "性别", "肿瘤数", "板数", "拆分", "选中细胞"], animalsData.map(r =>
  [r.verified_animal_block, ["12周", "18–20周", "30周"][r.group_number - 1], r.available_cells,
    r.weeks, r.sex, r.tumors, r.plates, r.split === "fit" ? "训练" : "测试", r.selected_cells]));
animals.getRange("A6:A20").format.columnWidth = 29;
animals.getRange("A4").values = [["肿瘤数按GEO字段计数；不将同鼠肿瘤当作独立动物。"]];

const aucSheet = sheet("判别AUC", "留出动物的一对其余AUC");
aucSheet.getRange("A4").values = [["点估计描述当前留出样本的判别能力；各鼠AUC共享负例，不是独立置信区间。"]];
const aucRows = tableObjects(data.auc).map(r => [short[`MCA_${r.learner}`], ["12周", "18–20周", "30周"][r.target - 1], r.auc, r.standardized_T]);
table(aucSheet, 6, ["学习器", "目标组", "测试AUC", "标准化统计量"], aucRows);
aucSheet.getRange("C7:D18").setNumberFormat("0.000000");
const animalAucRows = tableObjects(data.animal_auc).map(r => [short[`MCA_${r.learner}`], r.target, r.animal, r.auc]);
table(aucSheet, 22, ["学习器", "目标组编号", "目标组测试鼠", "该鼠对其余组AUC"], animalAucRows);
aucSheet.getRange("C22:C54").format.columnWidth = 29;
aucSheet.getRange("D23:D54").setNumberFormat("0.000000");
aucSheet.getRange("D6:D54").format.columnWidth = 24;

const genes = sheet("训练筛选基因", "训练细胞筛选的前200个基因");
genes.getRange("A4").values = [["主分析取前100；50/200维敏感性使用同一排序。未使用测试标签筛基因。"]];
const geneRows = tableObjects(data.genes).map(r => [r.variance_rank, r.geneID, r.ensgID, r.geneSymbol, r.fit_mean, r.fit_sd, r.fit_detection]);
table(genes, 6, ["排名", "geneID", "Ensembl", "基因符号", "训练均值", "训练标准差", "训练检出比例"], geneRows);
genes.getRange("A7:A206").setNumberFormat("0");
genes.getRange("B6:C206").format.columnWidth = 27;
genes.getRange("D6:D206").format.columnWidth = 22;
genes.getRange("E7:G206").setNumberFormat("0.000000");
genes.freezePanes.freezeRows(6);

const methods = sheet("方法说明", "分析定义与结果解释");
const explanations = [
  ["数据", "GSE154989；Marjanovic等（2020）；KP基因型时间序列"],
  ["表达单位", "作者提供log2(transcripts/100K+1)，不重复取log或归一化"],
  ["抽样与拆分", "每鼠等量无放回抽样；训练/测试按整鼠隔离"],
  ["筛选", "训练检出≥10%，排除mt-/Rpl/Rps，按训练log表达方差排序"],
  ["MCA-kDET", "L1、Quad、kNN、RF的一对其余标准化AUC最大值"],
  ["Sum / Comb", "Sum为各标准化AUC的有符号和；Comb为最大型与和型的联合排名"],
  ["Two-sided / Pair", "Two-sided对标准化AUC取绝对最大值；Pair为成对分类器最大值"],
  ["DISCO / MMD / MST", "同一测试矩阵，训练尺度；欧氏距离、Gaussian核、最小生成树"],
  ["整鼠置换", "固定测试组动物数，枚举560种分配；20周单组检查为210种"],
  ["原假设范围", "假设动物独立，抽样细胞块在原假设下可交换，强于单细胞边际同分布"],
  ["理论", "不直接声称现有iid细胞定理已覆盖动物块校准"],
  ["来源混杂", "18周雌性、20周雄性；取材、板和肿瘤结构也随阶段变化"],
  ["稳定性", "10次额外拆分报告全部p；共享动物，不是独立模拟或功效"],
  ["来源链接", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989"],
  ["原论文", "https://doi.org/10.1016/j.ccell.2020.06.012"],
];
table(methods, 6, ["项目", "定义"], explanations);
methods.getRange("A6:A21").format.columnWidth = 25;
methods.getRange("B6:B21").format.columnWidth = 112;

workbook.recalculate();
for (let i = 0; i < mainRows.length; i++) {
  const computed = main.getRange(`D${7 + i}`).values[0][0];
  const original = pvals.find(row => row.method === order[i]).p_value;
  if (typeof computed !== "number" || Math.abs(computed - original) > 1e-12) throw new Error("Workbook p-value formula mismatch");
}
console.log((await workbook.inspect({ kind: "table", range: "主结果!A6:E18", include: "values,formulas", tableMaxRows: 13, tableMaxCols: 5, maxChars: 3500 })).ndjson);
const preview = await workbook.render({ sheetName: "主结果", range: "A2:H18", scale: 1.5, format: "png" });
const previewDir = process.argv[2] || runtime;
await fs.mkdir(previewDir, { recursive: true });
await fs.writeFile(path.join(previewDir, "single_cell_workbook.png"), new Uint8Array(await preview.arrayBuffer()));
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(path.join(root, "results/单细胞实证结果.xlsx"));
console.log("Exported six-sheet research workbook.");
