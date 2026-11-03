"""Prepare a prespecified KP-stage analysis without fitting or testing models.

All random choices depend on IDs and a fixed seed, never expression/test results.
Gene selection and scaling summaries depend only on the selected fitting cells.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from prepare_data import ROOT, PREFIX, EXPRESSION, sha256


STAGES = {"06_KP_12w_ND": "KP_12w", "07_KP_20w_ND": "KP_18to20w", "08_KP_30w_ND": "KP_30w"}
FIT_COUNTS = {"KP_12w": 3, "KP_18to20w": 2, "KP_30w": 1}


def load_expression(genes: int, cells: int) -> sparse.csr_matrix:
    """Decode the author's MATLAB one-based sparse gene-by-cell coordinates."""
    with h5py.File(ROOT / "raw" / (PREFIX + EXPRESSION), "r") as handle:
        rows = handle["j"][0, :].astype(np.int32) - 1
        columns = handle["i"][0, :].astype(np.int32) - 1
        values = handle["v"][0, :].astype(np.float32)
    if rows.min() < 0 or rows.max() >= cells or columns.min() < 0 or columns.max() >= genes:
        raise ValueError("Expression coordinates outside companion metadata dimensions")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Invalid normalized expression values")
    matrix = sparse.coo_matrix((values, (rows, columns)), shape=(cells, genes)).tocsr()
    if matrix.nnz != len(values):
        raise ValueError("Duplicate sparse coordinates in source H5")
    return matrix


def choose_samples(metadata: pd.DataFrame, seed: int, cap: int, middle_20w_only: bool = False,
                   single_tumor: bool = False) -> pd.DataFrame:
    cohort = metadata.loc[metadata.author_stage.isin(STAGES)].copy()
    cohort["group"] = cohort.author_stage.map(STAGES)
    fit_counts = FIT_COUNTS.copy()
    if middle_20w_only:
        cohort = cohort.loc[~(cohort.group.eq("KP_18to20w") & cohort["geo_sack_timepoint (weeks)"].astype(int).eq(18))].copy()
        fit_counts["KP_18to20w"] = 1
    if single_tumor:
        tumor_rng = np.random.default_rng(seed + 179)
        for animal in sorted(cohort.loc[cohort.group.eq("KP_30w"), "verified_animal_block"].unique()):
            tumors = sorted(cohort.loc[cohort.verified_animal_block.eq(animal), "geo_tumor"].unique())
            chosen_tumor = tumor_rng.choice(tumors)
            cohort = cohort.loc[~cohort.verified_animal_block.eq(animal) | cohort.geo_tumor.eq(chosen_tumor)].copy()
    expected_strain = "KrasLSL-G12D; Trp53flox/flox; Rosa26LSL-tdTomato"
    if not cohort.geo_strain.eq(expected_strain).all():
        raise ValueError("Primary cohort contains more than one recorded genotype")
    weeks = cohort["geo_sack_timepoint (weeks)"].astype(int)
    expected_weeks = {"KP_12w": {12}, "KP_18to20w": {18, 20}, "KP_30w": {30}}
    for group, allowed in expected_weeks.items():
        if not set(weeks[cohort.group.eq(group)]).issubset(allowed):
            raise ValueError(f"Unexpected original collection time in {group}")
    cohort["split"] = "test"
    cohort["within_animal_random_rank"] = 0
    rng = np.random.default_rng(seed)
    for group in STAGES.values():
        animals = sorted(cohort.loc[cohort.group.eq(group), "verified_animal_block"].unique())
        chosen = set(rng.choice(animals, fit_counts[group], replace=False))
        cohort.loc[cohort.verified_animal_block.isin(chosen), "split"] = "fit"
        for animal in animals:
            indices = cohort.loc[cohort.verified_animal_block.eq(animal)].sort_values("sampleID").index.to_numpy()
            if len(indices) < cap:
                raise ValueError(f"{animal} has fewer than {cap} eligible cells")
            permuted = rng.permutation(indices)
            cohort.loc[permuted, "within_animal_random_rank"] = np.arange(1, len(indices) + 1)
    selected = cohort.loc[cohort.within_animal_random_rank.le(cap)].copy()
    if selected.groupby("verified_animal_block").split.nunique().max() != 1:
        raise ValueError("Animal leakage across fitting and test data")
    selected["group_number"] = selected.group.map({name: i + 1 for i, name in enumerate(STAGES.values())})
    selected = selected.sort_values(["split", "group_number", "verified_animal_block", "within_animal_random_rank"])
    if not selected.groupby("verified_animal_block").size().eq(cap).all():
        raise ValueError("Unequal selected cell counts per animal")
    return selected


def enumerate_animal_labels(test: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    """Enumerate label assignments only, without computing any test statistic."""
    units = test[["verified_animal_block", "group_number"]].drop_duplicates().sort_values("verified_animal_block")
    animals = units.verified_animal_block.tolist()
    counts = units.group_number.value_counts().to_dict()
    labelings = []
    for first in combinations(range(len(animals)), counts[1]):
        remaining = [index for index in range(len(animals)) if index not in first]
        for second in combinations(remaining, counts[2]):
            labels = np.full(len(animals), 3, dtype=np.int8)
            labels[list(first)] = 1
            labels[list(second)] = 2
            labelings.append(labels)
    labels = np.stack(labelings)
    expected = math.factorial(len(animals)) // math.prod(math.factorial(counts[g]) for g in [1, 2, 3])
    if len(labels) != expected or len(np.unique(labels, axis=0)) != expected:
        raise ValueError("Incorrect animal-label enumeration")
    observed = units.group_number.to_numpy()
    if np.sum(np.all(labels == observed, axis=1)) != 1:
        raise ValueError("Observed assignment must occur once in complete permutation space")
    return animals, labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--cells-per-mouse", type=int, default=30)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "prepared" / "analysis")
    parser.add_argument("--middle-20w-only", action="store_true")
    parser.add_argument("--single-tumor", action="store_true")
    args = parser.parse_args()
    metadata = pd.read_csv(ROOT / "prepared" / "cell_metadata.csv")
    genes = pd.read_csv(ROOT / "raw" / (PREFIX + "geneTable.csv.gz"))
    genes["matrix_column_0based"] = np.arange(len(genes))
    for column in ["geneID", "ensgID", "geneSymbol"]:
        genes[column] = genes[column].fillna("").astype(str).str.strip()
    if genes.geneID.duplicated().any():
        raise ValueError("Duplicate gene IDs")
    chosen = choose_samples(metadata, args.seed, args.cells_per_mouse, args.middle_20w_only, args.single_tumor)
    sample_indices = pd.Series(np.arange(len(metadata)), index=metadata.sampleID)
    chosen_indices = sample_indices.loc[chosen.sampleID].to_numpy()
    all_expression = load_expression(len(genes), len(metadata))
    expression = all_expression[chosen_indices, :]
    fitting = expression[chosen.split.eq("fit").to_numpy(), :].astype(np.float64)
    mean = np.asarray(fitting.mean(axis=0)).ravel()
    variance = (np.asarray(fitting.power(2).mean(axis=0)).ravel() - mean**2)
    variance = np.maximum(variance, 0) * fitting.shape[0] / (fitting.shape[0] - 1)
    detection = np.asarray((fitting > 0).mean(axis=0)).ravel()
    # This is explicitly log-expression variance ranking, not a claim to use
    # Seurat/scran's mean-adjusted highly-variable-gene algorithm.
    excluded = genes.geneSymbol.str.match(r"(?i)^(mt-|rpl|rps)") | genes.geneSymbol.eq("")
    eligible = (~excluded).to_numpy() & (detection >= 0.10) & (variance > 1e-12)
    ranking = genes.assign(fit_mean=mean, fit_sd=np.sqrt(variance), fit_detection=detection)
    ranking = ranking.loc[eligible].sort_values(["fit_sd", "geneID"], ascending=[False, True]).head(200).copy()
    if len(ranking) != 200:
        raise ValueError("Fewer than 200 eligible genes in fitting data")
    ranking.insert(0, "variance_rank", np.arange(1, len(ranking) + 1))
    columns = ranking.matrix_column_0based.to_numpy()
    selected_matrix = expression[:, columns].toarray()
    if not np.isfinite(selected_matrix).all():
        raise ValueError("Nonfinite selected expression values")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.output_dir / "genes.csv", index=False, encoding="utf-8-sig")
    chosen.to_csv(args.output_dir / "samples.csv", index=False, encoding="utf-8-sig")
    matrix_frame = pd.DataFrame(selected_matrix, columns=ranking.geneID)
    matrix_frame.insert(0, "sampleID", chosen.sampleID.to_numpy())
    matrix_frame.to_csv(args.output_dir / "expression.csv.gz", index=False, float_format="%.9g")
    summary = chosen.groupby(["group", "split"], sort=False).agg(
        animals=("verified_animal_block", "nunique"), cells=("sampleID", "size")
    ).reset_index()
    summary.to_csv(args.output_dir / "sample_counts.csv", index=False, encoding="utf-8-sig")
    animal_order, assignments = enumerate_animal_labels(chosen.loc[chosen.split.eq("test")])
    pd.DataFrame(assignments, columns=animal_order).to_csv(args.output_dir / "animal_label_assignments.csv", index=False)
    design = {
        "version": "KP_stage_preparation_v1", "seed": args.seed,
        "middle_20w_only": args.middle_20w_only, "single_tumor_per_late_mouse": args.single_tumor,
        "groups": list(STAGES.values()), "cells_per_mouse": args.cells_per_mouse,
        "primary_dimension": 100, "planned_dimension_sensitivity": [50, 200],
        "gene_selection": "Top log2(TP100K+1) fitting-cell variance; detected in >=10% fitting cells; exclude mt-/Rpl/Rps and blank symbols.",
        "normalization": "Author's published log2(transcripts/100K+1); do not normalize or log-transform a second time.",
        "matrix": "Unscaled log-expression, first 100 gene columns for primary analysis; standardize using fitting mean/SD only.",
        "all_genes_ranked_before_test_statistics": True,
        "planned_methods": ["MCA-kDET", "DISCO", "MMD", "MST"],
        "planned_calibration": "Complete animal-block label enumeration, subject to joint animal-block exchangeability; theoretical interpretation pending author discussion.",
        "test_animal_order": animal_order, "complete_animal_labelings": len(assignments),
        "minimum_attainable_complete_permutation_p": 1 / len(assignments),
        "inferential_tests_run": False,
        "source_expression_sha256": sha256(ROOT / "raw" / (PREFIX + EXPRESSION)),
        "metadata_sha256": sha256(ROOT / "prepared" / "cell_metadata.csv"),
        "fitting_sample_ids_sha256": __import__("hashlib").sha256("\n".join(chosen.loc[chosen.split.eq("fit"), "sampleID"]).encode()).hexdigest(),
        "scope": "Cross-sectional differences of sampled joint gene-expression distributions; no causal progression or dependence-only claim.",
    }
    (args.output_dir / "design.json").write_text(json.dumps(design, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"Prepared {len(chosen)} cells, 200 ranked genes; primary dimension 100.")
    print(f"Enumerated {len(assignments)} whole-animal labelings; no test statistics or p-values computed.")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
