"""Read original RDS in pure Python and independently recompute observed tests."""
from pathlib import Path
import argparse
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/"tests/out/financial_rdata_runtime"))
import rdata


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--directory",type=Path,default=ROOT/"real_analysis")
    parser.add_argument("--csv",default="results.csv")
    parser.add_argument("--expected-rows",type=int,default=480)
    parser.add_argument("--expected-jobs",type=int,default=40)
    parser.add_argument("--expected-primary-cells",type=int,default=12)
    args=parser.parse_args()
    out=args.directory;csv=pd.read_csv(out/args.csv)
    rows=[];checked_p=0;feature_cells=0
    features=None
    if (out/"sorted_sector_vectors.csv").exists():
        features=pd.read_csv(out/"sorted_sector_vectors.csv")
        sources=[ROOT/name for name in ("run_joint_analysis.R","joint_analysis/protocol.json","joint_analysis/sorted_sector_vectors.csv",
            "joint_analysis/fixed_schedules.csv","panel_calibration.R","native/panel_auc.dll")]
        sources += [ROOT.parents[1]/"R"/name for name in ("learners.R","learners_lib.R","tests.R")]
        expected=np.asarray(rdata.read_rds(out/"fingerprint.rds")["hashes"]).tolist()
        assert [hashlib.md5(p.read_bytes()).hexdigest() for p in sources]==expected
    for file in sorted((out/"jobs").glob("*.rds")):
        obj=rdata.read_rds(file)
        draw=int(obj["draw"][0]);q=int(obj["q"][0])
        for stride,a in obj["audits"].items():
            stats=a["statistics"]
            for method,values in stats.items():
                p=(1+np.count_nonzero(values[1:]>=values[0]-1e-10))/len(values)
                r=csv[(csv.draw==draw)&(csv.question==obj["rows"].question.iloc[0])&(csv.stride==int(stride))&(csv.method==method)].iloc[0]
                assert abs(p-r.p_exploratory)<1e-14
                checked_p+=1
            if draw!=1:continue
            X=np.asarray(a["X"]);S=np.asarray(a["S"]);L=np.asarray(a["Lall"])
            if features is not None:
                fit=features[(features.split=="fit")&features.complete_all_stocks.astype(str).str.lower().eq("true")][[f"x{i}" for i in range(1,6)]].to_numpy()
                st=a["coordinate_scale"];mu=np.asarray(st["mu"]);sd=np.asarray(st["sd"])
                np.testing.assert_allclose(fit.mean(axis=0),mu,atol=1e-13)
                np.testing.assert_allclose(fit.std(axis=0,ddof=1),sd,atol=1e-13)
                lookup=features[features.split=="test"].set_index(["time","sector"])
                rec=a["records"]
                raw=lookup.loc[list(zip(rec.time,rec.population)),[f"x{i}" for i in range(1,6)]].to_numpy()
                np.testing.assert_allclose((raw-mu)/sd,X,atol=1e-12)
                former=rdata.read_rds(ROOT/"real_analysis/jobs/draw01_q4.rds")["audits"][stride]
                assert np.array_equal(former["Lall"],L)
                feature_cells+=1
            targets=np.asarray(a["targets"]);N=len(X);k=int(L.max())
            n=np.array([(L[:,0]==i).sum() for i in range(1,k+1)])
            indicators=[(L==i).astype(float) for i in range(1,k+1)]
            ranks=rankdata(S,axis=0,method="average");Tm=np.zeros((S.shape[1],L.shape[1]))
            D=np.sqrt(np.maximum(cdist(X,X,"sqeuclidean"),0))
            K=np.exp(-D**2/(2*float(a["h2"][0])))
            disco=np.zeros(L.shape[1]);within=[];products=[]
            for i,I in enumerate(indicators):
                cols=np.where(targets==i+1)[0]
                auc=(ranks[:,cols].T@I-n[i]*(n[i]+1)/2)/(n[i]*(N-n[i]))
                Tm[cols]=(auc-.5)/np.sqrt((N+1)/(12*n[i]*(N-n[i])))
                disco-=(I*(D@I)).sum(axis=0)/n[i]
                products.append(K@I)
                within.append(((I*products[-1]).sum(axis=0)-n[i])/(n[i]*(n[i]-1)))
            mmd=np.zeros(L.shape[1]);pairs=0
            for i in range(k):
                for j in range(i+1,k):
                    mmd+=within[i]+within[j]-2*(indicators[i]*products[j]).sum(axis=0)/(n[i]*n[j]);pairs+=1
            mmd/=pairs
            # SciPy's graph implementation is independent of the R Prim routine.
            tree=minimum_spanning_tree(D).tocoo()
            assert tree.nnz==N-1
            mst=(L[tree.row]==L[tree.col]).sum(axis=0)
            values={"MCA-kDET":Tm.max(axis=0),"DISCO":disco,"MMD":mmd,"MST":mst}
            errors={m:float(np.max(np.abs(v-stats[m]))) for m,v in values.items()}
            assert max(errors.values())<1e-9,errors
            for method,v in values.items():
                p=(1+np.count_nonzero(v[1:]>=v[0]-1e-10))/len(v)
                original=(1+np.count_nonzero(stats[method][1:]>=stats[method][0]-1e-10))/len(v)
                assert p==original,(q,stride,method,p,original)
            rows.append(dict(question=q,stride=int(stride),permutation_labelings=L.shape[1],**errors))
    assert checked_p==args.expected_rows and len(rows)==args.expected_primary_cells
    assert len(list((out/"jobs").glob("*.rds")))==args.expected_jobs
    pd.DataFrame(rows).to_csv(out/"independent_python_method_audit.csv",index=False)
    result=dict(engine="Pure Python rdata parser, SciPy ranks/distances/MST, NumPy statistics; no R runtime",
                saved_job_files=args.expected_jobs,recomputed_pvalues=checked_p,all_pvalues_recomputed=True,
                all_480_pvalues_recomputed=(checked_p==480),primary_cells=len(rows),
                labelings_per_cell=5000,all_four_methods_recomputed_on_all_saved_primary_labelings=True,
                max_statistic_error=max(max(row[m] for m in ("MCA-kDET","DISCO","MMD","MST")) for row in rows),
                all_primary_pvalues_identical=True,
                independently_reconstructed_feature_cells=feature_cells,
                R_audit_runtime_note="The project's R runtime has intermittent decompression/internal errors; this audit uses a pure Python parser and independent numerical libraries.")
    (out/"independent_python_verification.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))


if __name__=="__main__":main()
