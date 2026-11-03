"""Fixed daily OHLC/turnover trial; one stock-date is one five-feature observation."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'daily_features_v1'
FEATURES = ['overnight_log_return', 'intraday_log_return', 'log_high_low_range', 'close_position', 'turnover_today']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def protocol():
    panel = pd.read_csv(ROOT/'design_two_tests/main_panel.csv')
    return dict(status='Exploratory feasibility trial; dependent-data significance not validated',
        author_request='One stock-day per observation; daily price shape and activity features. No GARCH in main analysis.',
        cohort=panel[['order_book_id','symbol','first_industry_name']].to_dict('records'),
        selection='Reuse fixed 2017 historical cohort, not current sector universe; no replacement by observed results.',
        fetch_start='2022-12-15', fetch_end='2024-12-31', fitting=['2023-01-01','2023-12-31'],
        testing='2024; first 20 market trading dates reserved as buffer, all remaining valid daily observations',
        features=dict(zip(FEATURES,['100*log(adjusted open / previous market trading-day adjusted close)',
            '100*log(adjusted close / adjusted open)', '100*log(adjusted high / adjusted low)',
            '(close-low)/(high-low); 0.5 if high==low on a traded day',
            'RQData get_turnover_rate today, vendor units retained, no rolling mean'])),
        scaling='Pooled fitting-only scale, shared across stocks. No per-stock or per-period centering/scaling.',
        missing='No price/feature interpolation. Current and previous date must be traded and unsuspended. Keep original calendar.',
        hypotheses='Within each fixed sector: equality of the five stocks daily five-feature marginal distributions.',
        diagnostics='Price-only first four coordinates as prespecified descriptive/learner diagnostic, no independent formal discoveries.',
        methods=['MCA-kDET','DISCO','MMD','MST'], B=999,
        interpretation='Ordinary label permutations are an IID reference diagnostic only, not valid real-market p-values under dependence. Never claim significance or equivalence from this trial.',
        frozen_sources={str(p.relative_to(ROOT.parents[1])):sha(p) for p in [ROOT.parents[1]/'R'/n for n in ['models.R','calib_params.R','calib_params_v2.R','tests.R','learners.R','learners_lib.R']]})


def fetch(secret_stdin):
    import rqdatac as rq
    secret = os.environ.get('RQDATA_LICENSE')
    if not secret and secret_stdin:
        print('Waiting for process-only API credential on stdin; it will not be logged or saved.', flush=True)
        secret = sys.stdin.readline().strip()
    if not secret:
        raise SystemExit('Process credential unavailable; existing files kept.')
    try:
        rq.init(username='license', password=secret)
    except Exception as e:
        raise SystemExit('RQData login failed: '+type(e).__name__) from None
    del secret
    p = protocol(); ids = [r['order_book_id'] for r in p['cohort']]
    raw = OUT/'raw'; raw.mkdir(parents=True, exist_ok=True)
    for name, call in [
        ('prices_post.csv.gz', lambda:rq.get_price(ids+['000300.XSHG'],start_date=p['fetch_start'],end_date=p['fetch_end'],
            frequency='1d',fields=['open','close','high','low','volume','total_turnover'],adjust_type='post',skip_suspended=False,expect_df=True)),
        ('turnover_today.csv.gz', lambda:rq.get_turnover_rate(ids,start_date=p['fetch_start'],end_date=p['fetch_end'],fields=['today'],expect_df=True)),
        ('suspended.csv.gz', lambda:rq.is_suspended(ids,start_date=p['fetch_start'],end_date=p['fetch_end']))]:
        path = raw/name
        if not path.exists():
            try:
                data = call()
                if data is None or len(data)==0:raise ValueError('Empty result')
                data.to_csv(path,float_format='%.17g',compression={'method':'gzip','mtime':0})
            except Exception as e:
                raise SystemExit('Fetch failed for '+name+': '+type(e).__name__) from None
        print('Cached '+name, flush=True)
    audit = dict(provider='RQData; existing user authorization; no credential persisted',
        retrieved_utc=datetime.now(timezone.utc).isoformat(),
        docs=['https://www.ricequant.com/doc/rqdata/python/generic-api','https://www.ricequant.com/doc/rqdata/python/stock-mod'],
        request=p, hashes={q.name:sha(q) for q in raw.glob('*.gz')})
    (OUT/'source_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def prepare():
    p = protocol(); panel = pd.DataFrame(p['cohort']); ids=panel.order_book_id.tolist()
    data=pd.read_csv(OUT/'raw/prices_post.csv.gz',parse_dates=['date'])
    turn=pd.read_csv(OUT/'raw/turnover_today.csv.gz');turn=turn.rename(columns={'tradedate':'date'})
    turn['date']=pd.to_datetime(turn.date)
    assert not data.duplicated(['date','order_book_id']).any()
    assert not turn.duplicated(['date','order_book_id']).any()
    days=pd.DatetimeIndex(sorted(data.loc[data.order_book_id.eq('000300.XSHG'),'date'].unique()))
    susp=pd.read_csv(OUT/'raw/suspended.csv.gz',index_col=0,parse_dates=True).reindex(index=days,columns=ids)
    fields={f:data.pivot(index='date',columns='order_book_id',values=f).reindex(index=days,columns=ids) for f in ['open','close','high','low','volume']}
    T=turn.pivot(index='date',columns='order_book_id',values='today').reindex(index=days,columns=ids)
    traded=fields['volume'].gt(0)&susp.eq(False)
    valid=traded&traded.shift(1,fill_value=False)&T.notna()&T.ge(0)
    for f in ['open','close','high','low']:valid &= fields[f].gt(0)&fields[f].notna()
    valid &= fields['close'].shift(1).gt(0)
    eps=1e-10
    bad=(fields['low']>fields['high']+eps)|(fields['open']<fields['low']-eps)|(fields['open']>fields['high']+eps)|(fields['close']<fields['low']-eps)|(fields['close']>fields['high']+eps)
    assert not (bad&valid).any().any(), 'OHLC ordering requires source review'
    O,C,H,L=[fields[f] for f in ['open','close','high','low']]
    spread=H-L
    position=(C-L)/spread.where(spread.ne(0));position=position.mask(spread.eq(0),.5)
    values=[100*np.log(O/C.shift(1)),100*np.log(C/O),100*np.log(H/L),position,T]
    assert np.allclose((values[0]+values[1]).where(valid), (100*np.log(C/C.shift(1))).where(valid),equal_nan=True,atol=1e-12)
    test_days=days[days.year==2024];buffer=set(test_days[:20]);test=set(test_days[20:])
    rows=[]
    for j,stock in enumerate(ids):
        meta=panel.iloc[j]
        for t,date in enumerate(days):
            split='fit' if date.year==2023 else 'buffer' if date in buffer else 'test' if date in test else 'warmup'
            r=dict(date=str(date.date()),stock=stock,name=meta.symbol,sector=meta.first_industry_name,split=split,valid=bool(valid.iloc[t,j]),
                equal_high_low=bool(spread.iloc[t,j]==0),traded_today=bool(traded.iloc[t,j]),traded_previous=bool(traded.shift(1,fill_value=False).iloc[t,j]))
            r.update({f:float(v.iloc[t,j]) if pd.notna(v.iloc[t,j]) else None for f,v in zip(FEATURES,values)})
            rows.append(r)
    frame=pd.DataFrame(rows);assert not frame.duplicated(['date','stock']).any()
    frame.to_csv(OUT/'daily_features.csv.gz',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    kept=frame[frame.valid&frame.split.isin(['fit','test'])]
    kept.to_csv(OUT/'complete_daily_features.csv.gz',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    counts=frame[frame.split.isin(['fit','test'])].groupby(['sector','stock','name','split']).agg(calendar_days=('date','size'),valid_days=('valid','sum'),equal_high_low_days=('equal_high_low','sum')).reset_index()
    counts.to_csv(OUT/'sample_counts.csv',index=False)
    desc=kept.groupby(['sector','stock','name','split'])[FEATURES].agg(['mean','std','median']).reset_index()
    desc.columns=['_'.join(c).rstrip('_') if isinstance(c,tuple) else c for c in desc.columns]
    desc.to_csv(OUT/'feature_descriptions.csv',index=False)
    kept.groupby('stock',sort=False).head(3).to_csv(OUT/'sample_preview.csv',index=False)
    quality=dict(no_duplicate_stock_dates=True,number_features=5,valid_rows=len(kept),invalid_rows=int((~frame.valid).sum()),
        first_test_date=str(test_days[20].date()),last_test_date=str(test_days[-1].date()),buffer_trading_days=20,
        return_decomposition_verified=True,ohlc_order_verified=True,turnover_units='RQData today original vendor units, not manually relabelled as a fraction',
        close_position_flat_day_rule=.5,source_hashes={q.name:sha(q) for q in (OUT/'raw').glob('*.gz')})
    (OUT/'data_quality.json').write_text(json.dumps(quality,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(quality),flush=True)
    print(counts[['sector','name','split','valid_days']].to_string(index=False),flush=True)


def main():
    a=argparse.ArgumentParser();a.add_argument('--fetch',action='store_true');a.add_argument('--secret-stdin',action='store_true');args=a.parse_args()
    OUT.mkdir(exist_ok=True)
    path=OUT/'protocol.json';p=protocol()
    if path.exists():assert json.loads(path.read_text(encoding='utf-8'))==p
    else:path.write_text(json.dumps(p,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if args.fetch:fetch(args.secret_stdin)
    prepare()


if __name__=='__main__':main()
