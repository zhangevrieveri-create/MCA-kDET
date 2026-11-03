"""Small fixed-cohort event application: source retrieval and calendar audit."""
from pathlib import Path
import csv
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from datetime import datetime,timezone
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"market_events";RAW=OUT/"raw"
URL="https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/49_Industry_Portfolios_daily_CSV.zip"


def french():
    path=RAW/"49_Industry_Portfolios_daily_CSV.zip"
    if not path.exists():
        request=urllib.request.Request(URL,headers={"User-Agent":"Academic research data preparation"})
        with urllib.request.urlopen(request,timeout=60)as response:data=response.read()
        if len(data)>30_000_000:raise ValueError("Unexpectedly large archive")
        zipfile.ZipFile(io.BytesIO(data)).testzip()
        path.write_bytes(data)
    with zipfile.ZipFile(path)as z:
        assert z.testzip() is None
        names=[n for n in z.namelist()if n.lower().endswith('.csv')];assert len(names)==1
        text=z.read(names[0]).decode('utf-8-sig')
    lines=text.splitlines();start=next(i for i,line in enumerate(lines)if 'Average Value Weighted Returns -- Daily' in line)
    header=next(i for i in range(start+1,start+10)if lines[i].startswith(','))
    cols=next(csv.reader([lines[header]]));cols=[c.strip()for c in cols];cols[0]='date'
    values=[]
    for line in lines[header+1:]:
        if not re.match(r'^\d{8},',line):break
        values.append(next(csv.reader([line])))
    frame=pd.DataFrame(values,columns=cols).set_index('date').astype(float)
    frame.index=pd.to_datetime(frame.index,format='%Y%m%d')
    frame=frame.mask(frame.isin([-99.99,-999]))
    frame=100*np.log1p(frame/100)
    assert len(frame)>10000 and frame.index.is_monotonic_increasing and len(frame.columns)==49
    frame.loc['2017':'2020'].to_csv(RAW/'us_industry_log_returns.csv.gz',index_label='date',compression='gzip')
    return frame,dict(url=URL,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),archive_bytes=path.stat().st_size,
      source_start=str(frame.index.min().date()),source_end=str(frame.index.max().date()),columns=frame.columns.tolist(),
      table='Daily value-weighted industry returns',return_transform='100*log(1+r_percent/100)',
      no_vix_in_China=True)


def main():
    RAW.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((OUT/'protocol.json').read_text(encoding='utf-8'))
    us,source=french()
    old=ROOT.parent/'financial/raw'
    prices=pd.read_csv(old/'daily_prices.csv.gz',parse_dates=['date'])
    close=prices.pivot(index='date',columns='order_book_id',values='close')
    volume=prices.pivot(index='date',columns='order_book_id',values='volume')
    calendar=close.index[close['000300.XSHG'].notna()]
    close=close.reindex(calendar);volume=volume.reindex(calendar)
    cn=100*np.log(close/close.shift(1));cn=cn.where(volume.gt(0))
    meta=pd.read_csv(old/'instruments.csv').merge(pd.read_csv(ROOT/'raw/hs300_industry_20141231.csv'),on='order_book_id')
    meta=meta.merge(pd.read_csv(old/'historical_weights.csv'),on='order_book_id').sort_values(['weight','order_book_id'],ascending=[False,True])
    modules=[]
    for sector,count in protocol['China_modules'].items():
        subset=meta[meta.first_industry_name==sector].head(count);assert len(subset)==count
        modules.append(dict(market='China',module=sector,assets=subset.order_book_id.tolist(),names=subset.symbol.tolist()))
    for name,assets in protocol['US_modules'].items():
        assert set(assets)<=set(us.columns)
        modules.append(dict(market='US',module=name,assets=assets,names=assets))
    prepared=[];phase_rows=[];feasibility=[]
    for event in ('China2015','China2024','US2020'):
        spec=protocol[event];market='US'if event.startswith('US')else'China';data=us if market=='US'else cn
        if 'phases'in spec:
            periods=[data.index[(data.index>=a)&(data.index<=b)]for a,b in spec['phases']]
        else:
            anchor=pd.Timestamp(spec['anchor']);pos=data.index.get_loc(anchor)
            periods=[data.index[pos+a:pos+b+1]for a,b in spec['trading_offsets']]
        first=min(d.min()for d in periods);last=max(d.max()for d in periods)
        for phase,days in enumerate(periods,1):
            fit_count=int(np.floor(.4*len(days)));parts=['classifier_fit']*fit_count+['buffer']*5+['test']*(len(days)-fit_count-5)
            assert len(parts)==len(days)
            phase_rows.extend(dict(event=event,phase=phase,label=spec['labels'][phase-1],date=str(day.date()),split=part)
                for day,part in zip(days,parts))
        for module in [m for m in modules if m['market']==market]:
            members=module['assets'];frame=data.loc[spec['baseline_fit'][0]:last,members].copy()
            prefix=f"{event}_m{len([a for a in prepared if a['event']==event])+1:02d}"
            path=OUT/f'{prefix}_returns.csv';frame.to_csv(path,index_label='date',float_format='%.17g')
            prepared.append(dict(event=event,module=module['module'],assets=members,names=module['names'],dimension=len(members),
              returns_file=path.name,baseline_fit=spec['baseline_fit'],phase_first=str(first.date()),phase_last=str(last.date())))
            for phase,days in enumerate(periods,1):
                fit_count=int(np.floor(.4*len(days)));valid=frame.reindex(days).notna().all(axis=1).to_numpy()
                nf=int(valid[:fit_count].sum());nt=int(valid[fit_count+5:].sum())
                feasibility.append(dict(event=event,module=module['module'],phase=phase,label=spec['labels'][phase-1],
                  calendar_days=len(days),complete_days=int(valid.sum()),classifier_fit_n=nf,test_n=nt,
                  first=str(days.min().date()),last=str(days.max().date()),estimable=(nf>=12 and nt>=12)))
    pd.DataFrame(phase_rows).to_csv(OUT/'phase_calendar.csv',index=False)
    pd.DataFrame(feasibility).to_csv(OUT/'feasibility.csv',index=False)
    (OUT/'modules.json').write_text(json.dumps(prepared,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    index=[]
    for i,m in enumerate(prepared,1):
        index.append(dict(job=i,event=m['event'],module=m['module'],dimension=m['dimension'],returns_file=m['returns_file'],
          baseline_start=m['baseline_fit'][0],baseline_end=m['baseline_fit'][1],assets=';'.join(m['assets'])))
    pd.DataFrame(index).to_csv(OUT/'modules_index.csv',index=False)
    sources=dict(downloaded_utc=datetime.now(timezone.utc).isoformat(),US=source,China=dict(provider='Previously authorized RQData cache',
      raw_file=str((old/'daily_prices.csv.gz').relative_to(ROOT.parent)),sha256=hashlib.sha256((old/'daily_prices.csv.gz').read_bytes()).hexdigest(),
      cohort='Fixed 2014 historical selection; no event-outcome replacement',role='Selected historical stock baskets, not every firm in each current sector'),
      protocol_sha256=hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),
      event_sources=dict(China2015='https://tjj.sh.gov.cn/tjnj/2016tjnj/C1713.htm',China2024='https://m.thepaper.cn/newsDetail_forward_28840070',
                        US2020='https://www.federalreserve.gov/newsevents/pressreleases/monetary20200323b.htm'))
    (OUT/'source_audit.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    f=pd.DataFrame(feasibility)
    print(f[['event','module','phase','calendar_days','classifier_fit_n','test_n','estimable']].to_string(index=False))
    print('No p-values computed. Modules:',len(prepared),'Public archive bytes:',source['archive_bytes'])


if __name__=='__main__':main()
