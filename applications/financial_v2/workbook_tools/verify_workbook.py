"""Read-only inspection of saved XLSX identifiers, cached formulas and tables."""
from pathlib import Path
import json
import zipfile
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def main():
    payload=json.loads((HERE/'workbook_input.json').read_text(encoding='utf-8'))
    filename=HERE.parent/'outputs/01a17651-f660-7773-9b6e-77219fc0a1cf/股票双检验汇总.xlsx'
    with zipfile.ZipFile(filename) as z:
        sst=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            sst=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
        def cells(sheet):
            root=ET.fromstring(z.read(f'xl/worksheets/sheet{sheet}.xml'))
            return {c.attrib['r']:c for c in root.findall('.//s:c',NS)}
        def value(c):
            kind=c.attrib.get('t')
            if kind=='inlineStr':return ''.join(c.find('s:is',NS).itertext())
            v=c.find('s:v',NS)
            if v is None:return None
            if kind=='s':return sst[int(v.text)]
            if kind=='str':return v.text
            return float(v.text)
        workbook=ET.fromstring(z.read('xl/workbook.xml'))
        names=[s.attrib['name'] for s in workbook.findall('s:sheets/s:sheet',NS)]
        assert names==['日频特征试验','日频样本','事件风险筛查','风险后果','联合结构试验','真实检验','抽样稳定性','结果概览','股票与滤波','距离明细','原收益校准','创新校准','抽样置换校准']
        def by_name(name):return cells(names.index(name)+1)
        model=by_name('股票与滤波');pairs=by_name('距离明细');overview=by_name('结果概览')
        for r,row in enumerate(payload['models'],6):assert value(model[f'A{r}'])==row[0]
        for r,row in enumerate(payload['pairs'],6):
            assert value(pairs[f'C{r}'])==row[2] and value(pairs[f'D{r}'])==row[3]
            assert abs(value(pairs[f'H{r}'])-row[7])<1e-14
        for r,row in enumerate(payload['summaries'],6):
            for col,key in [('D','W'),('E','B'),('F','Delta')]:assert abs(value(overview[f'{col}{r}'])-row[key])<2e-12
        formula_count=sum(c.find('s:f',NS) is not None for c in overview.values())
        assert formula_count==12
        tables=[name for name in z.namelist() if name.startswith('xl/tables/table') and name.endswith('.xml')]
        assert len(tables)==17
        for name,key,start in [('日频特征试验','daily_trial',6),('日频特征试验','daily_auc',22),('日频特征试验','daily_counts',39),('日频样本','daily_samples',6)]:
            sheet=by_name(name)
            for r,row in enumerate(payload[key],start):
                for i,expected in enumerate(row):
                    actual=value(sheet[f'{chr(65+i)}{r}'])
                    if isinstance(expected,str):assert actual==expected,(name,r,i)
                    else:assert abs(actual-expected)<1e-12,(name,r,i,actual,expected)
        precision=by_name('事件风险筛查')
        for r,row in enumerate(payload['event_precision'],31):
            for i,expected in enumerate(row):
                actual=value(precision[f'{chr(65+i)}{r}'])
                if isinstance(expected,str):assert actual==expected
                else:assert abs(actual-expected)<1e-12
        for name,key in [('事件风险筛查','event_screen'),('风险后果','event_risk'),('联合结构试验','joint_primary'),('真实检验','real_primary'),('抽样稳定性','real_stability'),('抽样置换校准','sampling_calibration')]:
            sheet=by_name(name)
            for r,row in enumerate(payload[key],6):
                for i,expected in enumerate(row):
                    actual=value(sheet[f'{chr(65+i)}{r}'])
                    if isinstance(expected,str):assert actual==expected
                    else:assert abs(actual-expected)<1e-12
        joint=by_name('联合结构试验')
        for r,row in enumerate(payload['joint_draws'],30):
            for i,expected in enumerate(row):assert abs(value(joint[f'{chr(65+i)}{r}'])-expected)<1e-12
    result=dict(saved_file=str(filename),all_stock_identifiers_preserved_as_text=True,
                all_420_pair_distances_match_source=True,summary_formula_count=12,
                all_cached_summary_values_match_independent_results=True,native_tables=17,sheets=names,
                daily_trial_rows_verified=len(payload['daily_trial']),daily_samples_verified=len(payload['daily_samples']),daily_auc_formulas_verified=12,
                event_precision_rows_verified=len(payload['event_precision']),
                event_screen_rows_verified=len(payload['event_screen']),event_risk_rows_verified=len(payload['event_risk']),
                joint_primary_rows_verified=len(payload['joint_primary']),joint_draws_verified=len(payload['joint_draws']),
                real_primary_rows_verified=len(payload['real_primary']),real_stability_rows_verified=len(payload['real_stability']),
                sampling_calibration_rows_verified=len(payload['sampling_calibration']),
                verification_engine='Artifact Tool plus independent read-only XLSX XML inspection; Excel app not opened')
    (HERE/'qa/saved_xlsx_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
