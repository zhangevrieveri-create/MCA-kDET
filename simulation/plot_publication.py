"""Publication-sized power figures from existing MCA-kDET simulation chunks.

Run from sim/: python plot_publication.py --out results_v2_quick
Change --out to the actual full-run directory when available. No model fitting,
simulation, source-data mutation, or manuscript editing is performed.

Layout references: the user's Yu et al. (2022), Jiang et al. (2022), and Guo
et al. reference PDFs; typography/line conventions also informed by SciencePlots
(https://github.com/garrettj403/SciencePlots). SciencePlots is not a dependency.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, norm


ROOT = Path(__file__).resolve().parent
ALPHA = 0.05
MODELS = {"A": "M1", "B": "M2", "C": "M3", "C2": "M4", "D": "M5", "E": "M6"}
DESCRIPTIONS = {
    "A": "Sparse mean shift", "B": "Variance shift", "C": "Pair dependence",
    "C2": "Asymmetric dependence", "D": "Sparse correlated block", "E": "Distributed mean shifts",
}
# A method keeps its colour, marker and line style wherever it appears.
STYLE = {
    "mca": ("MCA-kDET", "#171717", "o", "-"),
    "disco": ("DISCO", "#0072B2", "^", (0, (5, 2))),
    "mmd": ("MMD", "#E69F00", "s", (0, (1.5, 1.5))),
    "mst": ("MST", "#008B67", "D", (0, (5, 2, 1, 2))),
    "schott": ("Schott", "#AA4499", "v", (0, (7, 2))),
    "macc": ("Accuracy", "#7E7E7E", "x", (0, (2, 2))),
    "mca_l1": ("L1", "#0072B2", "^", (0, (5, 2))),
    "mca_quad": ("Quadratic", "#E69F00", "s", (0, (1.5, 1.5))),
    "mca_knn": ("kNN", "#56B4E9", "D", (0, (5, 2, 1, 2))),
    "mca_rf": ("Random forest", "#882255", "v", (0, (7, 2))),
    "pair": ("Pairwise", "#D55E00", "s", (0, (5, 2))),
    "shared_nn": ("Shared network", "#AA4499", "D", (0, (1.5, 1.5))),
    "mca_sum": ("Sum", "#0072B2", "^", (0, (5, 2))),
    "mca_comb": ("Combined", "#008B67", "D", (0, (5, 2, 1, 2))),
    "mca_2s": ("Two-sided", "#AA4499", "v", (0, (1.5, 1.5))),
    "mca_ms_med": ("Five splits: median", "#D55E00", "s", (0, (7, 2))),
    "mca_ms_mean": ("Five splits: mean", "#AA4499", "D", (0, (1.5, 1.5))),
}
CORE = ["mca", "disco", "mmd", "mst", "schott", "macc"]
LEARNERS = ["mca", "mca_l1", "mca_quad", "mca_knn", "mca_rf"]
VARIANTS = ["mca", "mca_sum", "mca_comb", "mca_2s", "pair", "shared_nn"]


def configure_style():
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 9.0, "axes.labelsize": 10.0, "axes.titlesize": 9.6,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.linewidth": 0.65, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "legend.fontsize": 8.7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.unicode_minus": False,
    })


def read_experiment(exp: str, directory: Path) -> pd.DataFrame:
    grid = pd.read_csv(directory / f"{exp}_grid.csv")
    cells = set(grid["cell"])
    chunks = sorted(p for p in (directory / exp).glob("*.csv") if p.name.split("__r")[0] in cells)
    if not chunks:
        raise ValueError(f"No chunks for {exp} in {directory}")
    df = pd.concat([pd.read_csv(p, low_memory=False) for p in chunks], ignore_index=True)
    if df.duplicated(["cell", "rep"]).any():
        raise ValueError(f"Duplicate repetition identifiers in {exp}")
    if "err" in df and df["err"].fillna("").astype(str).str.strip().ne("").any():
        raise ValueError(f"Failed replications in {exp}; inspect the original output")
    for row in grid.itertuples(index=False):
        subset = df[df.cell == row.cell]
        expected = int(row.reps)
        if len(subset) != expected or set(subset.rep) != set(range(1, expected + 1)):
            raise ValueError(f"Incomplete cell {row.cell}: {len(subset)} vs {expected}")
    print(f"{exp}: {len(grid)} design cells; {len(df)} valid replications")
    return df


def aggregate_cells(experiments):
    records = []
    for exp, df in experiments.items():
        for cell, subset in df.groupby("cell", sort=False):
            first = subset.iloc[0]
            for method in STYLE:
                if method not in subset:
                    continue
                values = pd.to_numeric(subset[method], errors="coerce").dropna()
                if not len(values):
                    continue  # An uncomputed method is missing, never zero.
                if not values.between(0, 1).all():
                    raise ValueError(f"Invalid p-values in {cell}, {method}")
                rejects = int((values <= ALPHA).sum())
                record = {"record_type": "simulation_cell", "experiment": exp, "cell": cell,
                          "model": first.model, "method": method, "n_reps": len(values),
                          "n_reject": rejects, "power": rejects / len(values)}
                for key in ["beta", "k", "n", "d", "tau", "nsplit", "sC"]:
                    record[key] = first[key] if key in first else np.nan
                records.append(record)
    return pd.DataFrame.from_records(records)


def curve_style(method):
    label, colour, marker, linestyle = STYLE[method]
    return dict(label=label, color=colour, marker=marker, linestyle=linestyle,
                linewidth=1.9 if method == "mca" else 1.05,
                markersize=4.2 if method == "mca" else 3.7,
                markeredgewidth=0.75,
                markerfacecolor=colour if method == "mca" else ("none" if marker == "x" else "white"),
                zorder=5 if method == "mca" else 3)


def make_canvas(nrows, ncols, methods, footer, height=None):
    height = height or (7.1 if nrows == 3 else 3.15)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5, height), squeeze=False, sharey=True)
    fig.subplots_adjust(left=.085, right=.982, bottom=.095 if nrows == 3 else .23,
                        top=.88 if nrows == 3 else .755, wspace=.22, hspace=.54)
    handles = [Line2D([], [], **{k: v for k, v in curve_style(m).items() if k != "zorder"}) for m in methods]
    fig.legend(handles=handles, labels=[STYLE[m][0] for m in methods], loc="upper center",
               bbox_to_anchor=(.52, .995), ncol=3, frameon=False, handlelength=2.8,
               columnspacing=1.8, handletextpad=.55, labelspacing=.6)
    fig.text(.085, .018 if nrows == 3 else .025, footer, fontsize=8.0, color="#656565")
    for r in range(nrows):
        for c in range(ncols):
            ax = axes[r, c]
            ax.set_ylim(-.018, 1.025)
            ax.set_yticks([0, .25, .5, .75, 1])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.2f}"))
            ax.axhline(ALPHA, color="#BBBBBB", linestyle=(0, (2, 3)), linewidth=.6, zorder=0)
            ax.grid(axis="y", color="#ECECEC", linewidth=.45, zorder=0)
            ax.set_axisbelow(True)
            ax.tick_params(direction="out", pad=3)
            ax.spines["left"].set_color("#5A5A5A")
            ax.spines["bottom"].set_color("#5A5A5A")
            if c == 0:
                ax.set_ylabel("Power", labelpad=6)
    return fig, axes


def draw_curves(ax, data, xkey, methods, title, xlabel, logx=False):
    xvals = np.sort(data[xkey].unique())
    for method in reversed(methods):
        subset = data[data.method == method].sort_values(xkey)
        if subset.empty:
            continue
        if subset.duplicated(xkey).any():
            raise ValueError(f"More than one design cell per point: {title}, {method}")
        ax.plot(subset[xkey], subset.power, **curve_style(method))
    if logx:
        ax.set_xscale("log")
        ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(FixedLocator(xvals))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.margins(x=.055)
    ax.set_xlabel(xlabel, labelpad=4)
    ax.set_title(title, loc="left", pad=8.5)


def uncertainty_interval(rejects, count):
    lo = 0.0 if rejects == 0 else beta_dist.ppf(.025, rejects, count - rejects + 1)
    hi = 1.0 if rejects == count else beta_dist.ppf(.975, rejects + 1, count - rejects)
    return float(lo), float(hi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results_v2_quick", help="Existing result directory")
    parser.add_argument("--dest", default="写作/latex/figures/power_review", help="Final vector figures")
    parser.add_argument("--preview-dir", required=True, help="PNG previews outside the project")
    args = parser.parse_args()
    directory = Path(args.out)
    if not directory.is_absolute(): directory = ROOT / directory
    destination = Path(args.dest)
    if not destination.is_absolute(): destination = ROOT / destination
    previews = Path(args.preview_dir)
    destination.mkdir(parents=True, exist_ok=True)
    previews.mkdir(parents=True, exist_ok=True)
    configure_style()
    experiments = {e: read_experiment(e, directory) for e in ["E2", "E3", "E4", "E5", "E6"]}
    agg = aggregate_cells(experiments)
    alternative = agg[agg.beta > 0]
    nmin, nmax = int(alternative.n_reps.min()), int(alternative.n_reps.max())
    if "quick" in directory.name.lower():
        profile = "Quick run"
    elif nmin >= 1000:
        profile = "Full run"
    else:
        profile = "Simulation results"
    count_label = str(nmin) if nmin == nmax else f"{nmin}--{nmax}"
    footer = f"{profile}: {count_label} replications per point; nominal level 0.05."
    figure_records = []

    def save(fig, name, caption):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for text in fig.findobj(matplotlib.text.Text):
            if not text.get_visible() or not text.get_text(): continue
            box = text.get_window_extent(renderer)
            bounds = fig.bbox
            if box.width > 0 and box.height > 0 and (box.x0 < -1 or box.y0 < -1 or box.x1 > bounds.x1 + 1 or box.y1 > bounds.y1 + 1):
                raise ValueError(f"Clipped text in {name}: {text.get_text()!r}")
        pdf = destination / f"{name}.pdf"
        fig.savefig(pdf, metadata={"Title": caption, "Author": "MCA-kDET simulation study"})
        fig.savefig(previews / f"{name}.png", dpi=220)
        bundle.savefig(fig)
        figure_records.append((name, caption))
        plt.close(fig)

    with PdfPages(destination / "power_figures_review.pdf", metadata={"Title": "MCA-kDET power figures: publication-sized preview"}) as bundle:
        signal = alternative[alternative.experiment.isin(["E2", "E4"])].copy()
        signal["signal"] = np.where(signal.model == "D", np.ceil(100 * signal.beta - 1e-9), signal.beta)

        fig, axes = make_canvas(3, 2, CORE, footer + " Schott was not computed for M6.")
        for i, model in enumerate(MODELS):
            data = signal[signal.model == model]
            draw_curves(axes.flat[i], data, "signal", CORE,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}",
                        r"Signal strength $\beta$" if model != "D" else r"Dependent coordinates $s$")
        save(fig, "fig1_power_signal", "Power against signal strength in all six models; E2 and E4. Six comparison methods; Schott is missing, not zero, in M6.")

        fig, axes = make_canvas(1, 2, CORE, footer + " Horizontal axes use a logarithmic scale.")
        for i, model in enumerate(["A", "C"]):
            data = alternative[(alternative.experiment == "E6") & (alternative.model == model)]
            draw_curves(axes.flat[i], data, "d", CORE,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}", r"Dimension $d$", logx=True)
        save(fig, "fig2_power_dimension", "Power against dimension in M1 and M3, beta=0.8 and k=3; E6. Logarithmic dimension axes; only simulated dimensions are plotted.")

        kmethods = ["mca", "pair", "mca_l1", "mca_rf", "shared_nn", "macc"]
        fig, axes = make_canvas(1, 2, kmethods, footer)
        for i, model in enumerate(["A", "C"]):
            data = alternative[(alternative.experiment == "E3") & (alternative.model == model)]
            draw_curves(axes.flat[i], data, "k", kmethods,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}", r"Number of populations $k$")
        save(fig, "fig3_power_k", "Power against the number of populations in M1 and M3, beta=0.8; E3. Every point represents a separate experiment.")

        fig, axes = make_canvas(3, 2, LEARNERS, footer)
        for i, model in enumerate(MODELS):
            draw_curves(axes.flat[i], signal[signal.model == model], "signal", LEARNERS,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}",
                        r"Signal strength $\beta$" if model != "D" else r"Dependent coordinates $s$")
        save(fig, "figS1_power_learners", "Library versus the four individual learners across all six models; E2 and E4.")

        fig, axes = make_canvas(3, 2, VARIANTS, footer + " Pairwise was not computed for M6.")
        for i, model in enumerate(MODELS):
            draw_curves(axes.flat[i], signal[signal.model == model], "signal", VARIANTS,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}",
                        r"Signal strength $\beta$" if model != "D" else r"Dependent coordinates $s$")
        save(fig, "figS2_power_variants", "Maximum, sum, combined, two-sided, pairwise and shared-network tests across signal strengths; E2 and E4. Missing pairwise values in M6 are not drawn.")

        tmethods = ["mca", "mca_sum", "mca_comb", "mca_ms_med", "mca_ms_mean"]
        fig, axes = make_canvas(1, 2, tmethods, footer + " Five-split rules use twice the median/mean p-value.")
        for i, model in enumerate(["A", "C"]):
            data = alternative[(alternative.experiment == "E5") & (alternative.model == model)]
            draw_curves(axes.flat[i], data, "tau", tmethods,
                        f"({chr(97+i)}) {MODELS[model]}: {DESCRIPTIONS[model]}", r"Fitting fraction $\tau$")
        save(fig, "figS3_power_split", "Power against the fitting fraction in M1 and M3, beta=0.8; E5. Median/mean rules combine five splits as defined in the manuscript.")

        # This is a diagnostic conditional on the observed test-set AUC, not an
        # additional unconditional-power comparison. Preserve that distinction.
        aucc = [f"auc_{m}_dev" for m in ["l1", "quad", "knn", "rf"]]
        pool = pd.concat([experiments["E2"].query("beta > 0"),
                          experiments["E3"].query("k == 3"), experiments["E6"]], ignore_index=True)
        best = pool[aucc].max(axis=1)
        if best.isna().any():
            raise ValueError("Missing observed AUCs in the diagnostic")
        lower_edge = math.floor(float(best.min()) * 100) / 100
        edges = np.arange(lower_edge, 1.0101, .01)
        diagnostic = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = pool[(best >= lo) & (best < hi)]
            if not len(sub): continue
            hits = int((sub.mca <= ALPHA).sum()); ci = uncertainty_interval(hits, len(sub))
            diagnostic.append(dict(record_type="auc_diagnostic_bin", experiment="E2/E3(k=3)/E6",
                                   model="pooled", method="mca", auc_lower=lo, auc_upper=hi,
                                   auc_midpoint=(lo+hi)/2, n_reps=len(sub), n_reject=hits,
                                   power=hits/len(sub), ci_lower=ci[0], ci_upper=ci[1]))
        diag = pd.DataFrame(diagnostic)
        if int(diag.n_reps.sum()) != len(pool):
            raise ValueError("The AUC bins do not cover every diagnostic observation")
        visible = diag[diag.n_reps >= 10]
        fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.15))
        fig.subplots_adjust(left=.09, right=.98, top=.84, bottom=.24, wspace=.38)
        ax = axes[0]
        ax.errorbar(visible.auc_midpoint, visible.power,
                    yerr=np.vstack([visible.power-visible.ci_lower, visible.ci_upper-visible.power]),
                    fmt="o-", color="#171717", markersize=3.3, linewidth=1.1,
                    ecolor="#858585", elinewidth=.65, capsize=1.7)
        threshold = .5 + norm.ppf(1-ALPHA/12)*math.sqrt(451/(12*150*300))
        ax.axvline(threshold, color="#0072B2", linestyle=(0, (4, 2)), linewidth=.8)
        ax.set_ylabel("Rejection frequency within bin")
        ax.set_ylim(-.02, 1.025)
        ax.set_title("(a) Observed-AUC diagnostic", loc="left")
        axes[1].bar(diag.auc_midpoint, diag.n_reps, width=.008, color="#899DAA", linewidth=0)
        axes[1].set_ylabel("Number of replications")
        axes[1].set_title("(b) Replications supporting each bin", loc="left")
        for ax in axes:
            ax.set_xlabel("Largest observed test-set AUC")
            ax.set_xlim(lower_edge - .005, 1.005)
            ax.set_xticks([.5, .6, .7, .8, .9, 1])
        fig.text(.09, .035, f"{profile}; {len(pool)} pooled replications. Exact 95% binomial intervals; bins with n < 10 omitted in (a).", fontsize=7.6, color="#656565")
        save(fig, "figS4_auc_diagnostic", "Observed-AUC diagnostic, not an unconditional-power curve. Panel (a): conditional rejection frequencies and exact binomial intervals for bins with at least 10 observations. Panel (b): bin counts. Dashed line is a Bonferroni approximation, not an observed permutation critical value.")

    combined = pd.concat([agg, diag], ignore_index=True)
    combined.to_csv(destination / "plot_data.csv", index=False, float_format="%.10g")
    notes = ["# Power figure preview", "", f"Source: {directory}",
             f"Alternative replications per point: {count_label} ({profile}).", "",
             "The PDF pages are 6.5 inches wide; use them at full text width.",
             "All curves connect observed simulation points; there is no smoothing or invented grid.",
             "This plotting script does not edit raw results, original figures, manuscript sources, or hand-entered full-run Tables 1/2.",
             "The dimension symbol is d, consistent with the manuscript (p in many other papers).", "",
             "Current main-text figures: signal (Figure 1) and dimension (Figure 2).",
             "Current supplement: learners (S1), variants (S2), k (S3), fitting fraction (S4), AUC diagnostic (S5).",
             "The exact k results remain in main-text Table 4; E3 has no DISCO/MMD/MST/Schott results.", "",
             "The AUC diagnostic is explicitly conditional on observed test AUC and pools experiments. It is optional and does not replace setting-wise power comparisons.", "",
             "## Figures", ""]
    notes.extend(f"- {name}.pdf: {caption}" for name, caption in figure_records)
    (destination / "README.md").write_text("\n".join(notes)+"\n", encoding="utf-8")
    print(f"Wrote {len(figure_records)} individual PDFs and one review PDF to {destination}")
    print(f"PNG previews: {previews}")
    print(f"Auditable plotting data: {len(combined)} records; {len(pool)} diagnostic observations")


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
