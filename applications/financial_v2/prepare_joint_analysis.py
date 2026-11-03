"""Build sector order-statistic vectors on the original uncompressed calendar."""
from pathlib import Path
import sys
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/"tests/out/financial_rdata_runtime"))
import rdata


def main():
    out=ROOT/"joint_analysis";out.mkdir(exist_ok=True)
    panel=pd.read_csv(ROOT/"design_two_tests/main_panel.csv")
    cal=pd.read_csv(ROOT/"design_two_tests/block_calendar.csv")
    cal=cal[cal.dimension==5]
    vectors=pd.read_csv(ROOT/"innovations/innovation_vectors.csv.gz",float_precision="round_trip")
    vectors=vectors[vectors.dimension==5]
    rows=[];components=[];checks=[]
    for split in ("fit","test"):
        calendar=cal[cal.split==split];blocks=calendar.block.to_numpy()
        for g,sector in enumerate(panel.first_industry_name.drop_duplicates(),1):
            ids=panel.loc[panel.first_industry_name==sector,"order_book_id"].tolist();assert len(ids)==5
            X=np.stack([vectors[(vectors.split==split)&(vectors.stock==stock)].set_index("block").loc[blocks,[f"r{i}" for i in range(1,6)]].to_numpy() for stock in ids],axis=1)
            weekly=X.sum(axis=2)/np.sqrt(5)
            ordered=np.sort(weekly,axis=1)
            valid=np.isfinite(X).all(axis=(1,2))
            complete=calendar.complete_all_main_stocks.astype(str).str.lower().eq("true").to_numpy()
            assert valid[complete].all() and np.all(np.diff(ordered[complete],axis=1)>=0)
            for i,block in enumerate(blocks):
                rows.append(dict(split=split,block=int(block),time=i+1,sector=g,sector_name=sector,
                    complete_all_stocks=bool(complete[i]),**{f"x{j+1}":ordered[i,j] for j in range(5)}))
                for j,stock in enumerate(ids):
                    components.append(dict(split=split,block=int(block),sector=g,stock=stock,weekly_score=weekly[i,j],complete_all_stocks=bool(complete[i])))
            checks.append(dict(split=split,sector=g,calendar_blocks=len(blocks),common_complete_blocks=int(complete.sum()),
                max_sum_preservation_error=float(np.max(np.abs(ordered[complete].sum(axis=1)-weekly[complete].sum(axis=1))))))
    pd.DataFrame(rows).to_csv(out/"sorted_sector_vectors.csv",index=False)
    pd.DataFrame(components).to_csv(out/"weekly_stock_scores.csv",index=False)
    schedules=[]
    for draw in range(1,11):
        job=rdata.read_rds(ROOT/f"real_analysis/jobs/draw{draw:02d}_q4.rds")
        for stride,audit in job["audits"].items():
            rec=audit["records"][["time","population"]].copy()
            assert not rec.time.duplicated().any() and len(set(rec.population.value_counts()))==1
            for r in rec.itertuples():schedules.append(dict(draw=draw,stride=int(stride),time=int(r.time),population=int(r.population)))
    pd.DataFrame(schedules).to_csv(out/"fixed_schedules.csv",index=False)
    sources=[ROOT/"innovations/innovation_vectors.csv.gz",ROOT/"design_two_tests/block_calendar.csv",ROOT/"design_two_tests/main_panel.csv",out/"protocol.json"]
    result=dict(representation="sorted weekly cross-stock residual scores",scope="exploratory",checks=checks,
        primary_fixed_before_new_pvalues=True,schedules_reused_without_outcome_selection=True,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (out/"preparation_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(rows=len(rows),weekly_stock_values=len(components),schedule_rows=len(schedules),checks=checks),ensure_ascii=False))


if __name__=="__main__":main()
