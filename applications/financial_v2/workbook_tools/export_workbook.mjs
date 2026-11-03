import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const here=path.dirname(fileURLToPath(import.meta.url));
const input=JSON.parse(await fs.readFile(path.join(here,'workbook_input.json'),'utf8'));
const root=path.dirname(here);
const out=path.join(root,'outputs','01a17651-f660-7773-9b6e-77219fc0a1cf');
const qa=path.join(here,'qa');
await fs.mkdir(out,{recursive:true});await fs.mkdir(qa,{recursive:true});
const wb=Workbook.create();
const names=['日频特征试验','日频样本','事件风险筛查','风险后果','联合结构试验','真实检验','抽样稳定性','结果概览','股票与滤波','距离明细','原收益校准','创新校准','抽样置换校准'];
const sheets=Object.fromEntries(names.map(name=>[name,wb.worksheets.add(name)]));
const navy='#173A53',muted='#506273',light='#EDF3F7';
function col(n){let s='';for(let k=n;k;k=Math.floor((k-1)/26))s=String.fromCharCode(65+(k-1)%26)+s;return s;}
function table(name,title,note,headers,rows,widths){
  const s=sheets[name],last=rows.length+5,end=col(headers.length);
  s.showGridLines=false;s.tabColor=navy;
  s.getRange(`A1:${end}${Math.max(last,5)}`).format.font={name:'Arial',size:10,color:'#243746'};
  s.getRange(`A1:${end}${Math.max(last,5)}`).format.verticalAlignment='center';
  s.getRange(`A1:${end}${Math.max(last,5)}`).format.rowHeightPx=23;
  s.getRange('A1').values=[[title]];s.getRange('A1').format.font={name:'Arial',size:14,bold:true,color:navy};
  s.getRange('A1').format.rowHeightPx=30;
  s.getRange('A3').values=[[note]];s.getRange('A3').format.font={name:'Arial',size:10,color:muted};
  s.getRange(`A5:${end}5`).values=[headers];
  if(rows.length){
    s.getRange(`A6:${end}${last}`).values=rows;
    s.tables.add(`A5:${end}${last}`,true,`Research_${names.indexOf(name)}`);
    for(let row=6;row<=last;row+=2)s.getRange(`A${row}:${end}${row}`).format.fill=light;
    s.getRange(`A6:${end}${last}`).format.horizontalAlignment='left';
  }
  s.getRange(`A5:${end}5`).format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},
      horizontalAlignment:'center',verticalAlignment:'center',wrapText:true};
  s.getRange(`A5:${end}5`).format.rowHeightPx=42;
  widths.forEach((width,i)=>s.getRange(`${col(i+1)}1:${col(i+1)}${Math.max(last,5)}`).format.columnWidthPx=width);
  if(name!=='结果概览')s.freezePanes.freezeRows(5);
  return s;
}
const dailyTrial=table('日频特征试验','一天一条多特征观测：行业内小试','2023拟合、2024前20交易日缓冲后逐日测试。B999；以下置换数值仅是独立样本参考，尚不能解释为正式显著性。',
  input.daily_trial_headers,input.daily_trial,[90,150,170,170,110,110,110,110,140,270]);
dailyTrial.getRange('E6:I11').setNumberFormat('0.0000');dailyTrial.getRange('E6:I11').format.horizontalAlignment='right';
dailyTrial.getRange('A14').values=[['主五维：隔夜涨跌、日内涨跌、高低振幅、收盘位置、当日换手率；不用GARCH，不构造五日块，不稀疏抽样。']];
dailyTrial.getRange('A15').values=[['参考值0.001是B999的分辨率下限，不是p=0；同日跨股票和时间相关未获校准，不标显著性星号。']];
dailyTrial.getRange('A16').values=[['RF为随机森林；宏平均AUC是五股各自一对其余AUC的等权平均。以下完整展示四种学习器，差值不是提升显著性检验。']];
dailyTrial.getRange('A14:J16').format.font={name:'Arial',size:10,color:muted};
dailyTrial.getRange('A19').values=[['换手率提供多少额外区分信息：全部学习器']];dailyTrial.getRange('A19').format.font={name:'Arial',size:12,bold:true,color:navy};
dailyTrial.getRange('A21:E21').values=[input.daily_auc_headers];dailyTrial.getRange('A22:E33').values=input.daily_auc;
dailyTrial.tables.add('A21:E33',true,'Daily_AUC_diagnostic');
dailyTrial.getRange('C22:E33').setNumberFormat('0.0000');dailyTrial.getRange('C22:E33').format.horizontalAlignment='right';
for(let r=22;r<=33;r++)dailyTrial.getRange(`E${r}`).formulas=[[`=C${r}-D${r}`]];
dailyTrial.getRange('A36').values=[['实际样本数：每只股票各一张日频矩阵']];dailyTrial.getRange('A36').format.font={name:'Arial',size:12,bold:true,color:navy};
dailyTrial.getRange('A38:E38').values=[input.daily_count_headers];dailyTrial.getRange('A39:E53').values=input.daily_counts;
dailyTrial.tables.add('A38:E53',true,'Daily_stock_counts');dailyTrial.getRange('B39:B53').setNumberFormat('000000');dailyTrial.getRange('D39:E53').setNumberFormat('0');
for(const range of ['A21:E21','A38:E38']){
  dailyTrial.getRange(range).format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true};dailyTrial.getRange(range).format.rowHeightPx=38;
}
dailyTrial.getRange('A21:E53').format.font={name:'Arial',size:10};dailyTrial.getRange('A22:E53').format.rowHeightPx=24;
const dailySamples=table('日频样本','实际输入数据：一只股票的一天为一行','保留全部有效拟合/测试日。价格统一后复权；前3个特征为100×对数，收盘位置0至1，换手率保留米筐today原单位。',
  input.daily_sample_headers,input.daily_samples,[115,100,110,85,75,145,145,130,130,150]);
dailySamples.getRange(`B6:B${input.daily_samples.length+5}`).setNumberFormat('000000');
dailySamples.getRange(`F6:J${input.daily_samples.length+5}`).setNumberFormat('0.0000');dailySamples.getRange(`F6:J${input.daily_samples.length+5}`).format.horizontalAlignment='right';
const eventScreen=table('事件风险筛查','市场事件：哪些固定风险模块发生分布变化','主检验为共用GARCH过滤后的向量；每方法跨全部15项作Holm。均为探索性数值，不构成因果或韧性结论。',
  input.event_headers,input.event_screen,[120,130,70,130,110,120,120,120,120,240]);
eventScreen.getRange('C6:C20').setNumberFormat('0');eventScreen.getRange('E6:I20').setNumberFormat('0.0000');
eventScreen.getRange('E6:I20').format.horizontalAlignment='right';
eventScreen.getRange('E6:I20').conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#FFF0D6',font:{color:'#8A4F00'}}});
eventScreen.getRange('B23').values=[['2015：上涨/高位见顶/回撤；9·24：较早前期/临近政策前/政策发布后；美股2020：较早前期/冲击累积/政策响应后。']];
eventScreen.getRange('B24').values=[['半数越界比例为原收益描述，不是检验p值；两维食品饮料模块中包括单个成分越界，不据此跨模块排名。']];
eventScreen.getRange('B25').values=[['上表四方法统一B4999；下方全部15项MCA统一加密到B49999，原结果保留。精度增加不代表功效增加。']];
eventScreen.getRange('A23:J25').format.font={name:'Arial',size:10,color:muted};
eventScreen.getRange('A28').values=[['MCA置换加密：B4999 → B49999（全部15项；窗口、数据、学习器及Holm家族不变）']];
eventScreen.getRange('A28').format.font={name:'Arial',size:13,bold:true,color:navy};
eventScreen.getRange('A30:I30').values=[input.event_precision_headers];
eventScreen.getRange('A31:I45').values=input.event_precision;
eventScreen.tables.add('A30:I45',true,'Event_MC_precision');
eventScreen.getRange('A30:I45').format.font={name:'Arial',size:10,color:'#243746'};
eventScreen.getRange('A30:I45').format.rowHeightPx=26;
eventScreen.getRange('A30:I30').format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true};
eventScreen.getRange('A30:I30').format.rowHeightPx=42;
eventScreen.getRange('C31:H45').setNumberFormat('0.00000');eventScreen.getRange('I31:I45').setNumberFormat('0');
eventScreen.getRange('C31:I45').format.horizontalAlignment='right';
eventScreen.getRange('A47').values=[['区间仅为固定数据/分数下置换尾概率的逐格95%蒙特卡洛区间；不是经济效应区间或同时区间。']];
eventScreen.getRange('A48').values=[['仍属探索性；时间依赖及过滤估计误差的有效性尚待验证。增加置换次数不会补足每阶段21–40条测试观测。']];
eventScreen.getRange('A47:J48').format.font={name:'Arial',size:10,color:muted};
const eventRisk=table('风险后果','原收益的阶段风险后果：完整描述','所有阶段均保留。组合收益为日等权简单收益，尾部日数较少；无置信区间，不作独立筛选后验证或投资排名。',
  input.event_risk_headers,input.event_risk,[120,130,175,100,125,160,100,150,150]);
eventRisk.getRange('D6:D50').setNumberFormat('0');eventRisk.getRange('G6:G50').setNumberFormat('0');eventRisk.getRange('E6:E50').setNumberFormat('0.0%');
eventRisk.getRange('F6:F50').setNumberFormat('0.000');eventRisk.getRange('H6:I50').setNumberFormat('0.000');eventRisk.getRange('D6:I50').format.horizontalAlignment='right';
const joint=table('联合结构试验','新探索：行业内个股、行业间周横截面','行业内保留原始p；行业间改用五股周残差排序向量，四问题Holm重新计算。固定draw1、间隔4。',
  input.real_primary_headers,input.joint_primary,[180,125,80,110,115,120,120,105,205]);
joint.getRange('C6:D21').setNumberFormat('0');joint.getRange('E6:G21').setNumberFormat('0.0000');joint.getRange('H6:H21').setNumberFormat('0');
joint.getRange('C6:H21').format.horizontalAlignment='right';
joint.getRange('E6:F21').conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#FFF0D6',font:{color:'#8A4F00'}}});
joint.getRange('A23:I42').format.font={name:'Arial',size:10,color:muted};
joint.getRange('B24').values=[['行业内仍比较每股五日残差分布；行业间比较排序后的周横截面联合结构，未跨股票平均。']];
joint.getRange('B25').values=[['新MCA行业间原始探索p=0.0478，Holm=0.1912；不作为已确认差异。新表示未获水平保证。']];
joint.getRange('A28').values=[['行业间固定间隔4：全部十次抽样（同一数据，非独立实证）']];
joint.getRange('A28').format.font={name:'Arial',size:11,bold:true,color:navy};
joint.getRange('A29:G29').values=[input.joint_draw_headers];joint.getRange('A30:G39').values=input.joint_draws;
joint.tables.add('A29:G39',true,'JointSamplingStability');
joint.getRange('A29:G29').format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center'};
joint.getRange('A29:G29').format.rowHeightPx=38;joint.getRange('A30:G39').format.rowHeightPx=23;
joint.getRange('A30:B39').setNumberFormat('0');joint.getRange('C30:G39').setNumberFormat('0.0000');joint.getRange('A30:G39').format.horizontalAlignment='right';
for(let r=30;r<=39;r+=2)joint.getRange(`A${r}:G${r}`).format.fill=light;
joint.getRange('B41').values=[['此前个股混合分布的全部结果保留在原页；不能按新p值改变股票、日期或主抽样。']];
const primary=table('真实检验','原表示：行业间个股混合分布','原表示结果完整保留。固定第1次抽样、间隔4；校准未通过，不当作确认性发现。',
  input.real_primary_headers,input.real_primary,[170,125,80,110,115,120,120,105,205]);
primary.getRange('C6:D21').setNumberFormat('0');primary.getRange('E6:G21').setNumberFormat('0.0000');primary.getRange('H6:H21').setNumberFormat('0');
primary.getRange('C6:H21').format.horizontalAlignment='right';
primary.getRange('E6:F21').conditionalFormats.add('cellIs',{operator:'lessThan',formula:.05,format:{fill:'#FFF0D6',font:{color:'#8A4F00'}}});
primary.getRange('A24').values=[['读表说明']];primary.getRange('B24').values=[['每种方法内调整四项问题；不作四方法共16项的联合错误控制。']];
primary.getRange('B25').values=[['MCA行业间Holm探索p=0.5016；MST=0.0232。MST仅1/10次抽样低于0.05。']];
primary.getRange('B26').values=[['三项行业内未拒绝不能证明同分布；完整十次抽样及三种间隔见下一页。']];
primary.getRange('A24:I26').format.font={name:'Arial',size:10,color:muted};
const stability=table('抽样稳定性','同一真实数据的十次抽样与三个间隔','各方法数值均为探索性。完整展示，固定主结果始终保留；不取最小p、不平均p，也不把十次当作独立实证。',
  input.real_stability_headers,input.real_stability,[100,105,160,110,115,115,115,115,115,210]);
stability.getRange('A6:B125').setNumberFormat('0');stability.getRange('D6:D125').setNumberFormat('0');stability.getRange('E6:I125').setNumberFormat('0.0000');
stability.getRange('D6:I125').format.horizontalAlignment='right';
const overview=table('结果概览','股票双检验：描述结果与校准结论','真实探索结果见前两页；九场景各1000次确认显示候选仍超水平，不能确认性解释。',
  input.overview_headers,input.summaries.map(s=>[s.representation,s.d,s.n,null,null,null]),[150,80,115,145,145,135]);
const distance=table('距离明细','统一核下的个股配对距离','非平方经验MMD点估计。各表示采用各自训练尺度和带宽，不直接比较距离下降。',
  input.pair_headers,input.pairs,[150,75,100,100,100,100,90,135,110]);
distance.getRange(`B6:B${input.pairs.length+5}`).setNumberFormat('0');
distance.getRange(`C6:D${input.pairs.length+5}`).setNumberFormat('000000');
distance.getRange(`H6:H${input.pairs.length+5}`).setNumberFormat('0.000000');
distance.getRange(`I6:I${input.pairs.length+5}`).setNumberFormat('0');
distance.getRange(`H6:I${input.pairs.length+5}`).format.horizontalAlignment='right';
const last=input.pairs.length+5;
for(let i=0;i<input.summaries.length;i++){
  const r=i+6;
  overview.getRange(`D${r}`).formulas=[[`=AVERAGEIFS('距离明细'!$H$6:$H$${last},'距离明细'!$A$6:$A$${last},$A${r},'距离明细'!$B$6:$B$${last},$B${r},'距离明细'!$G$6:$G$${last},"行业内")`]];
  overview.getRange(`E${r}`).formulas=[[`=AVERAGEIFS('距离明细'!$H$6:$H$${last},'距离明细'!$A$6:$A$${last},$A${r},'距离明细'!$B$6:$B$${last},$B${r},'距离明细'!$G$6:$G$${last},"行业间")`]];
  overview.getRange(`F${r}`).formulas=[[`=E${r}-D${r}`]];
}
overview.getRange('B6:C9').setNumberFormat('0');overview.getRange('D6:F9').setNumberFormat('0.000000');
overview.getRange('B6:F9').format.horizontalAlignment='right';
const notes=[['研究目标','GARCH调整后的创新分布'],['原收益推断','校准未通过，保留描述'],['创新推断','已计算探索p；压力校准未通过'],
  ['当前候选','完整训练面板＋测试块随机单股抽样'],['固定面板','银行、医药、机械各五股'],['训练与测试','训练2015–2017，测试2018–2024'],
  ['距离解释','点估计不等于显著性或等效性'],['跨表示比较','各表示的尺度和带宽不同'],['校准区间','逐格95%二项蒙特卡洛区间'],
  ['校准次数','九场景各1000次独立确认已完成'],['固定主间隔','每隔4个五日块抽样；仍有误报偏高'],
  ['测试样本量','每股15–16条、每行业25–26条；功效有限'],['保留不利结果','GARCH(0.999)家族误报0.092，目标0.05']];
overview.getRange('A13:B25').values=notes;
overview.getRange('A13:F25').format.font={name:'Arial',size:10,color:'#243746'};
overview.getRange('A13:A25').format.font={name:'Arial',size:10,bold:true,color:navy};
overview.getRange('A13:F25').format.rowHeightPx=25;

const model=table('股票与滤波','训练期GARCH参数与测试期诊断','均值及参数仅由2015–2017估计；边界标记=1表示持续性或系数接近约束。',
  input.model_headers,input.models,[100,125,85,115,110,90,90,100,130,130,95,115]);
model.getRange('D6:J20').setNumberFormat('0.0000');model.getRange('K6:L20').setNumberFormat('0');
model.getRange('A6:A20').setNumberFormat('000000');
model.getRange('D6:L20').format.horizontalAlignment='right';
model.getRange('K6:K20').conditionalFormats.add('cellIs',{operator:'greaterThan',formula:0,format:{fill:'#FFF0D6',font:{bold:true,color:'#8A4F00'}}});
model.getRange('A24').values=[['数据来源']];model.getRange('B24').values=[[input.source]];
model.getRange('A25').values=[['数据接口']];model.getRange('B25').values=[[input.source_url]];
model.getRange('A24:L25').format.font={name:'Arial',size:10,color:muted};

for(const [name,rows,title,note]of[
  ['原收益校准',input.raw_calibration,'原收益候选的模拟校准','α=0.05；保留初轮、长块和环形块全部结果。不同场景不能合并重复次数。'],
  ['创新校准',input.innovation_calibration,'真创新与估计残差的配对校准','α=0.05；同一模拟面板分别使用真创新及重新估计的过滤器。固定残差重抽样尚未控制全部误报。'],
  ['抽样置换校准',input.sampling_calibration,'随机单股抽样：1000次确认及100次试跑','α=0.05，B=999；固定主间隔4，候选未通过。标注功效的行为备择；家族仅统计真原假设。']]){
  const s=table(name,title,note,name==='抽样置换校准'?input.sampling_headers:input.calibration_headers,rows,[240,145,160,120,95,95,100,100,100]);
  if(rows.length){
    s.getRange(`D6:F${rows.length+5}`).setNumberFormat('0');
    s.getRange(`G6:I${rows.length+5}`).setNumberFormat('0.000');
    s.getRange(`D6:I${rows.length+5}`).format.horizontalAlignment='right';
    s.getRange(`G6:I${rows.length+5}`).conditionalFormats.addCustom('AND($C6="Holm家族误报",$H6>0.05)',{fill:'#FBE7E4',font:{bold:true,color:'#9D302A'}});
  }
}
wb.recalculate();
const actual=overview.getRange('D6:F9').values;
for(let i=0;i<4;i++)for(let j=0;j<3;j++){
  const expected=[input.summaries[i].W,input.summaries[i].B,input.summaries[i].Delta][j];
  if(typeof actual[i][j]!=='number'||Math.abs(actual[i][j]-expected)>2e-12)throw new Error(`Summary formula mismatch ${i},${j}: ${actual[i][j]} vs ${expected}`);
}
for(let i=0;i<input.daily_auc.length;i++){
  const actual=dailyTrial.getRange(`E${i+22}`).values[0][0];
  if(typeof actual!=='number'||Math.abs(actual-input.daily_auc[i][4])>1e-14)throw new Error('Daily AUC difference formula mismatch');
}
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:50},summary:'final formula error scan'});
await fs.writeFile(path.join(qa,'formula_error_scan.ndjson'),errors.ndjson);
const qaArgument=process.argv.find(a=>a.startsWith('--qa-sheets='));
const qaNames=qaArgument?qaArgument.slice('--qa-sheets='.length).split(','):names;
for(const name of qaNames){
  if(!names.includes(name))throw new Error(`Unknown QA sheet ${name}`);
  const range=name==='日频特征试验'?'A1:J16':name==='日频样本'?'A1:J19':name==='事件风险筛查'?'A1:J26':name==='风险后果'?'A1:I21':name==='联合结构试验'?'A1:I42':name==='结果概览'?'A1:F26':name==='股票与滤波'?'A1:L21':name==='真实检验'?'A1:I27':name==='抽样稳定性'?'A1:J18':'A1:I18';
  const preview=await wb.render({sheetName:name,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(qa,`${name}.png`),new Uint8Array(await preview.arrayBuffer()));
  if(name==='事件风险筛查'){
    const precisionPreview=await wb.render({sheetName:name,range:'A28:J49',scale:1.3,format:'png'});
    await fs.writeFile(path.join(qa,'事件风险筛查_置换精度.png'),new Uint8Array(await precisionPreview.arrayBuffer()));
  }
  if(name==='日频特征试验'){
    const diagnostic=await wb.render({sheetName:name,range:'A19:F54',scale:1.3,format:'png'});
    await fs.writeFile(path.join(qa,'日频特征试验_样本与AUC.png'),new Uint8Array(await diagnostic.arrayBuffer()));
  }
}
const file=await SpreadsheetFile.exportXlsx(wb);
const filename=path.join(out,'股票双检验汇总.xlsx');await file.save(filename);
try{await fs.rename(`${filename}.inspect.ndjson`,path.join(qa,'saved_workbook_inspect.ndjson'));}catch(e){if(e.code!=='ENOENT')throw e;}
await fs.writeFile(path.join(qa,'verification.json'),JSON.stringify({summary_formulas_match_independent_results:true,sheets:names,
  imported_rows:{event_screen:input.event_screen.length,event_risk:input.event_risk.length,joint_primary:input.joint_primary.length,joint_draws:input.joint_draws.length,real_primary:input.real_primary.length,real_stability:input.real_stability.length,models:input.models.length,pairs:input.pairs.length,raw_calibration:input.raw_calibration.length,innovation_calibration:input.innovation_calibration.length,sampling_calibration:input.sampling_calibration.length},
  actual_summary_values:actual},null,2));
console.log(JSON.stringify({output:filename,sheets:names.length,summary_formulas_verified:true}));
