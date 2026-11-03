"""Audit forward GARCH recursion and retain the original vector calendar."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    src,out = ROOT/"design_two_tests",ROOT/"innovations"
    panel = pd.read_csv(src/"main_panel.csv")
    ids = panel.order_book_id.tolist()
    raw = pd.read_csv(src/"daily_returns.csv.gz",index_col="date",parse_dates=True,float_precision="round_trip").loc["2015":,ids]
    Z = pd.read_csv(out/"daily_innovations.csv.gz",index_col="date",parse_dates=True,float_precision="round_trip")
    H = pd.read_csv(out/"conditional_variance.csv.gz",index_col="date",parse_dates=True,float_precision="round_trip")
    parameters = pd.read_csv(out/"garch_parameters.csv",float_precision="round_trip").set_index("stock")
    assert raw.index.equals(Z.index) and raw.index.equals(H.index)
    assert Z.columns.tolist()==ids and H.columns.tolist()==ids
    np.testing.assert_array_equal(raw.isna(),Z.isna())
    fit = raw.index.year<=2017
    audits=[]
    for stock in ids:
        p=parameters.loc[stock]
        mu = raw.loc[fit,stock].mean();v0=((raw.loc[fit,stock]-mu)**2).mean()
        np.testing.assert_allclose([mu,v0],[p['mean'],p.initial_variance],rtol=1e-12,atol=1e-12)
        h=[];last_h=last_x2=v0
        for x in raw[stock].to_numpy()-p['mean']:
            current=p.omega+p.alpha*last_x2+p.beta*last_h
            h.append(current);last_h=current;last_x2=x*x if np.isfinite(x) else current
        h=np.array(h)
        np.testing.assert_allclose(h,H[stock],rtol=2e-10,atol=2e-10)
        z=(raw[stock].to_numpy()-p['mean'])/np.sqrt(h)
        np.testing.assert_allclose(z,Z[stock],rtol=2e-10,atol=2e-10,equal_nan=True)
        known=fit & raw[stock].notna().to_numpy()
        nll=.5*np.sum(np.log(h[known])+(raw[stock].to_numpy()[known]-p['mean'])**2/h[known])
        assert abs(nll-p.negative_log_likelihood)<1e-7
        audits.append(dict(stock=stock,max_variance_relative_error=float(np.max(np.abs(h-H[stock])/H[stock])),
                           fitting_observations=int(known.sum()),near_boundary=int(p.near_boundary)))
    calendar=pd.read_csv(src/"block_calendar.csv")
    rows=[]
    for block in calendar.itertuples():
        values=Z.loc[block.start:block.end,ids]
        assert len(values)==block.dimension
        common=bool(values.notna().all().all())
        assert common==block.complete_all_main_stocks
        for stock in ids:
            rows.append(dict(dimension=block.dimension,split=block.split,block=block.block,stock=stock,
                             industry=panel.set_index("order_book_id").loc[stock,"first_industry_name"],
                             complete=bool(values[stock].notna().all()),complete_all_main_stocks=common,
                             **{f"r{i+1}":v for i,v in enumerate(values[stock].to_numpy())}))
    pd.DataFrame(rows).to_csv(out/"innovation_vectors.csv.gz",index=False,float_format="%.17g",compression={"method":"gzip","mtime":0})
    diagnostics=pd.read_csv(out/"autocorrelation.csv")
    summary={}
    for rep in ("raw","filtered"):
        x=diagnostics[(diagnostics.representation==rep)&(diagnostics.power==2)&(diagnostics.lag==1)].autocorrelation
        summary[rep]=dict(median_squared_acf1=float(x.median()),max_abs_squared_acf1=float(x.abs().max()))
    correlations=pd.read_csv(out/"daily_correlation_filtered.csv",index_col=0).to_numpy()
    pairs=correlations[np.triu_indices(15,1)]
    result=dict(user_authorized_target="GARCH-adjusted innovation distributions; raw returns descriptive only",
                real_classifiers_fitted=False,real_pvalues_computed=False,calibration_validated=False,
                independent_recursion_audit=audits,calendar_and_missingness_match_original=True,
                fitting_days=int(fit.sum()),test_days=int((raw.index>=pd.Timestamp("2018-01-30")).sum()),
                daily_residual_cross_correlation=dict(median=float(np.median(pairs)),max_abs=float(np.abs(pairs).max())),
                squared_acf1=summary,
                sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [src/"daily_returns.csv.gz",src/"block_calendar.csv",out/"daily_innovations.csv.gz",out/"innovation_vectors.csv.gz",out/"garch_parameters.csv"]})
    (out/"preparation_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("fitting_days","test_days","squared_acf1","daily_residual_cross_correlation")},ensure_ascii=False))
    print("Independent GARCH recurrence, QML objective, fitting means and complete vector calendar verified for 15 stocks.")


if __name__=="__main__":
    main()
