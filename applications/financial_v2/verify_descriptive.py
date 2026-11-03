"""Independent full pooled-Gram reconstruction of descriptive stock distances."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--representation",choices=["raw","innovations"],default="raw")
    args=parser.parse_args()
    src, out = ROOT/"design_two_tests", ROOT/("innovations/descriptive" if args.representation=="innovations" else "descriptive")
    panel = pd.read_csv(src/"main_panel.csv")
    path=ROOT/"innovations/innovation_vectors.csv.gz" if args.representation=="innovations" else src/"return_vectors.csv.gz"
    raw = pd.read_csv(path, float_precision="round_trip")
    findings = []
    for d in (5,20):
        info = json.loads((out/f"distance_summary_d{d}.json").read_text(encoding="utf-8"))
        for name, expected in info["source_sha256"].items():
            source=ROOT/info["source_paths"][name] if "source_paths" in info else src/name
            assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
        data = raw[(raw.dimension == d) & raw.complete_all_main_stocks]
        stocks = panel.order_book_id.tolist()
        arrays = {}
        for split in ("fit","test"):
            s = data[data.split == split]
            arrays[split] = np.array([s[s.stock == stock].sort_values("block")[[f"r{j+1}" for j in range(d)]].to_numpy() for stock in stocks])
        train = arrays["fit"].reshape(-1,d)
        np.testing.assert_array_equal(train.mean(axis=0),info["fitting_mean"])
        np.testing.assert_array_equal(train.std(axis=0,ddof=1),info["fitting_sd"])
        X = ((arrays["test"]-train.mean(axis=0))/train.std(axis=0,ddof=1)).reshape(-1,d)
        J, N = arrays["test"].shape[:2]
        # Independent implementation: one full pooled kernel, weighted by stock IDs.
        norm2 = np.sum(X*X,axis=1)
        squared = np.maximum(norm2[:,None]+norm2[None,:]-2*(X @ X.T),0)
        K = np.exp(-squared/(2*info["fitting_median_bandwidth"]**2))
        P = np.kron(np.eye(J),np.ones((1,N))/N)
        embedding = P @ K @ P.T
        recovered = np.sqrt(np.maximum(np.diag(embedding)[:,None]+np.diag(embedding)[None,:]-2*embedding,0))
        observed = pd.read_csv(out/f"distance_matrix_d{d}.csv",index_col=0).to_numpy()
        # Tiny roundoff on exact self-distances is separately set to its known zero.
        np.fill_diagonal(recovered,0)
        error = float(np.max(np.abs(recovered-observed)))
        assert error < 2e-12
        values_within,values_between = [],[]
        for j in range(J):
            for k in range(j+1,J):
                (values_within if panel.first_industry_name.iloc[j]==panel.first_industry_name.iloc[k] else values_between).append(recovered[j,k])
        assert len(values_within)==30 and len(values_between)==75
        np.testing.assert_allclose([np.mean(values_within),np.mean(values_between)], [info["W"],info["B"]],atol=2e-13,rtol=0)
        findings.append(dict(d=d,common_fit_blocks=arrays["fit"].shape[1],common_test_blocks=N,
                             all_source_hashes_match=True,shared_fitting_scale_exact=True,
                             independent_pooled_kernel_max_distance_error=error,
                             within_between_pair_counts_verified=True))
    (out/"independent_verification.json").write_text(json.dumps(findings,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(findings))


if __name__ == "__main__":
    main()
