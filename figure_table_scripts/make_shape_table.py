"""Main-text table of the shape alternatives (d = 20, oracle AUC 0.70/0.80) from
sim/motivated_v4_20261102/full_B200/report_tables/all_method_rates_and_intervals.csv."""
import csv, sys
from pathlib import Path
SIM = Path(__file__).resolve().parents[3] / "sim" if False else Path(__file__).resolve().parents[2].parent
src = SIM / "motivated_v4_20261102/full_B200/report_tables/all_method_rates_and_intervals.csv"
rows = list(csv.DictReader(open(src, encoding="utf-8")))
R = {}
for r in rows: R[(r["model"], r["target"], r["d"], r["method"])] = float(r["rate"])
OTHERS = ["maxmean", "clx_mean", "clx_cov", "maxmean_ovr", "clx_cov_ovr", "disco", "mmd", "mst",
          "schott", "schott_perm", "mmd_max", "kmd", "songchen", "songchen_A"]
NAMES = {"maxmean": "max-mean", "clx_mean": "CLX mean", "clx_cov": "max-cov", "maxmean_ovr": "max-mean",
         "clx_cov_ovr": "max-cov", "disco": "DISCO", "mmd": "MMD", "mst": "MST", "schott": "Schott",
         "schott_perm": "Schott", "mmd_max": "MMD-max", "kmd": "KMD", "songchen": "SC-S", "songchen_A": "SC-SA"}
lab = {"N1_bimodal": "(N1)", "N2_scale_mix": "(N2)", "N3_condvar": "(N3)", "N4_mixed": "(N4)",
       "N5_skewmix": "(N5)", "N6_ring": "(N6)"}
out = []
for m in lab:
    for t in ["0.7", "0.8"]:
        g = lambda meth: R[(m, t, "20", meth)]
        oth = [(g(o), o) for o in OTHERS if (m, t, "20", o) in R]
        best = max(oth)
        out.append((lab[m], t, g("mca_set_l1+quad+knn+rf"), g("mca_rf"), g("songchen"), g("songchen_A"),
                    g("kmd"), g("mmd_max"), g("disco"), best[0], NAMES[best[1]]))
lines = []
for o in out:
    lines.append("%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f (%s) \\\\" %
                 (o[0], float(o[1]), *o[2:10], o[10]))
Path(__file__).with_name("shape_table_rows.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
# averages over the 24 alternatives
alts = {(r["model"], r["target"], r["d"]) for r in rows if float(r["target"]) > 0.5}
for meth in ["mca_set_l1+quad+knn+rf", "mca_rf"]:
    print(meth, round(sum(R[a + (meth,)] for a in alts) / len(alts), 4), len(alts))
