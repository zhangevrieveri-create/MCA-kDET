"""Descriptive common-kernel distances. No tests or confidence statements."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, cdist

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def describe(d, out, vectors_path=None, representation="raw"):
    data = ROOT / "design_two_tests"
    panel = pd.read_csv(data / "main_panel.csv")
    vectors_path = vectors_path or data/"return_vectors.csv.gz"
    vectors = pd.read_csv(vectors_path, float_precision="round_trip")
    calendar = pd.read_csv(data / "block_calendar.csv")
    stocks = panel.order_book_id.tolist()
    sectors = panel.first_industry_name.drop_duplicates().tolist()
    arrays = {}
    for split in ("fit", "test"):
        cal = calendar[(calendar.dimension == d) & (calendar.split == split)]
        valid = cal.complete_all_main_stocks.astype(str).str.lower().eq("true")
        blocks = cal.loc[valid, "block"].tolist()
        subset = vectors[(vectors.dimension == d) & (vectors.split == split)]
        arrays[split] = np.stack([
            subset[subset.stock == stock].set_index("block").loc[blocks, [f"r{i+1}" for i in range(d)]].to_numpy()
            for stock in stocks])
    fit = arrays["fit"].reshape(-1, d)
    mu, sd = fit.mean(axis=0), fit.std(axis=0, ddof=1)
    if not np.isfinite(sd).all() or not (sd > 0).all():
        raise ValueError("Invalid shared fitting scale")
    sigma = float(np.median(pdist((fit - mu) / sd)))
    X = (arrays["test"] - mu) / sd
    J, N, _ = X.shape
    gram = np.empty((J, J))
    for j in range(J):
        for k in range(j, J):
            gram[j, k] = gram[k, j] = np.exp(-cdist(X[j], X[k], "sqeuclidean") / (2*sigma*sigma)).mean()
    distance2 = np.diag(gram)[:, None] + np.diag(gram)[None, :] - 2*gram
    if distance2.min() < -1e-12:
        raise ValueError("Negative squared kernel distance")
    D = np.sqrt(np.maximum(distance2, 0))
    np.testing.assert_allclose(D, D.T, atol=1e-14)
    pairs = []
    for j in range(J):
        for k in range(j+1, J):
            pairs.append(dict(d=d, stock_a=stocks[j], stock_b=stocks[k],
                              name_a=panel.symbol.iloc[j], name_b=panel.symbol.iloc[k],
                              sector_a=panel.first_industry_name.iloc[j], sector_b=panel.first_industry_name.iloc[k],
                              type="within" if panel.first_industry_name.iloc[j] == panel.first_industry_name.iloc[k] else "between",
                              mmd_unsquared=D[j, k], complete_test_calendar_blocks=N))
    pairs = pd.DataFrame(pairs)
    assert len(pairs) == 105 and (pairs.type == "within").sum() == 30
    weights = np.array([panel.first_industry_name.eq(s).to_numpy()/5 for s in sectors])
    sector_gram = weights @ gram @ weights.T
    sectorD = np.sqrt(np.maximum(np.diag(sector_gram)[:, None] + np.diag(sector_gram)[None, :] - 2*sector_gram, 0))
    # Independently check all three mixture distances using signed stock weights.
    for a in range(3):
        for b in range(a+1, 3):
            w = weights[a]-weights[b]
            np.testing.assert_allclose(sectorD[a,b]**2, w @ gram @ w, atol=2e-15)
    W = float(pairs.loc[pairs.type == "within", "mmd_unsquared"].mean())
    B = float(pairs.loc[pairs.type == "between", "mmd_unsquared"].mean())
    info = dict(d=d, representation=representation,scope="Descriptive point distances only; no calibrated uncertainty or p-values",
                kernel="exp(-||z-z'||^2/(2*sigma^2)); z=(return-shared_fit_mean)/shared_fit_sd",
                estimator="Norm of empirical kernel mean difference, including diagonal terms (V statistic)",
                common_fit_blocks=arrays["fit"].shape[1], common_test_blocks=N,
                fitting_mean=mu.tolist(), fitting_sd=sd.tolist(), fitting_median_bandwidth=sigma,
                within_stock_pairs=30, between_stock_pairs=75,
                W=W, B=B, Delta=B-W,
                within_by_sector={s:float(pairs.loc[(pairs.type == "within") & (pairs.sector_a == s), "mmd_unsquared"].mean()) for s in sectors},
                between_by_sector_pair=[dict(a=sectors[a],b=sectors[b],
                    mean_stock_pair_distance=float(pairs.loc[(pairs.sector_a==sectors[a]) & (pairs.sector_b==sectors[b]),"mmd_unsquared"].mean()))
                    for a in range(3) for b in range(a+1,3)],
                mixture_distances=[dict(a=sectors[a], b=sectors[b], distance=float(sectorD[a,b])) for a in range(3) for b in range(a+1,3)],
                source_sha256={p.name:sha(p) for p in (vectors_path, data/"main_panel.csv", data/"block_calendar.csv")},
                source_paths={p.name:str(p.relative_to(ROOT)) for p in (vectors_path,data/"main_panel.csv",data/"block_calendar.csv")})
    pairs.to_csv(out/f"stock_pair_distances_d{d}.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(D, index=stocks, columns=stocks).to_csv(out/f"distance_matrix_d{d}.csv", encoding="utf-8-sig")
    (out/f"distance_summary_d{d}.json").write_text(json.dumps(info, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key:info[key] for key in ("d", "common_test_blocks", "W", "B", "Delta", "within_by_sector", "mixture_distances")}, ensure_ascii=False))
    return info


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--representation",choices=["raw","innovations"],default="raw")
    args=parser.parse_args()
    innovations=args.representation=="innovations"
    out = ROOT/("innovations/descriptive" if innovations else "descriptive")
    out.mkdir(exist_ok=True,parents=True)
    path=ROOT/"innovations/innovation_vectors.csv.gz" if innovations else None
    results = [describe(d, out,path,args.representation) for d in (5,20)]
    lines = ["# 固定股票面板的描述性距离"+("（GARCH标准化残差）" if innovations else "（原始收益）"), "",
             "这里只汇报点估计，尚无经验证的置信区间或正式p值。不能把点估计的大小当成显著性或等效性。",
             "所有股票用同一个拟合段均值、标准差及高斯核带宽；没有逐股去均值/波动标准化。股票名单和顺序在看测试收益之前固定。",
             "每一维为单日复权对数收益，主分析为不重叠五交易日向量。共同缺失日期全部股票一起排除，完整日历仍在原数据保留。",
             "经验MMD为非平方核均值距离；有限样本下，即使真分布相同，点距离一般也不为零。横截面相关可能影响这些有限样本距离的比较。", ""]
    if innovations:
        lines[3]="数据先按各股训练期估计的GARCH去除均值与条件波动，再统一使用创新拟合段的公共坐标尺度/核带宽。它与原收益图采用不同表示和核尺度，不能把距离下降直接解释为改善幅度。"
        lines[4]="每一维为估计的单日标准化残差，主分析为原日历上非重叠五交易日向量；原缺失标记完全保留。估计残差不等于已知真创新，模型和第一阶段误差仍待校准。"
    for r in results:
        lines += [f"## d={r['d']}（共同训练/测试段 {r['common_fit_blocks']}/{r['common_test_blocks']}）", "",
                  f"行业内30对平均 W={r['W']:.6f}；行业间75对平均 B={r['B']:.6f}；差值 B−W={r['Delta']:.6f}。", "",
                  "行业内均值："+"；".join(f"{k} {v:.6f}" for k,v in r["within_by_sector"].items())+"。", "",
                  "各行业对的25个个股距离均值："+"；".join(f"{z['a']}—{z['b']} {z['mean_stock_pair_distance']:.6f}" for z in r["between_by_sector_pair"])+"。", "",
                  "行业混合分布距离（另一个量，不能替代75对个股距离均值）："+
                  "；".join(f"{z['a']}—{z['b']} {z['distance']:.6f}" for z in r["mixture_distances"])+"。", ""]
    lines += ["主目标只涵盖选定15只股票在共同可观测日期上的分布，不代表全部行业或停牌期间。",
              "股票对不是105个独立重复，不对它们直接作普通t检验；不因不符合理想结构而删股或改顺序。"]
    (out/"描述性距离.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
