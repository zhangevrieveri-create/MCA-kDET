"""Acquire and audit the published GSE154989 data; never run hypothesis tests.

Run with --expression to also download the published normalized expression H5.
Original downloads are kept unchanged. Outputs use sampleID for every join.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl/"
PREFIX = "GSE154989_mmLungPlate_fQC_"
METADATA = [
    "smpTable.csv.gz",
    "dZ_annot_smpTable.csv.gz",
    "dZ_QCstat_smpTable.csv.gz",
    "geneTable.csv.gz",
]
EXPRESSION = "dSp_normTPM.h5"
SOFT_NAME = "GSE154989_family.soft.gz"
SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/soft/" + SOFT_NAME


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire(name: str, url: str | None = None) -> dict:
    """Atomic downloads. Existing original files are checked and reused."""
    path = ROOT / "raw" / name
    url = url or GEO + name
    if not path.exists():
        part = path.with_suffix(path.suffix + ".part")
        print(f"Downloading {name}", flush=True)
        request = urllib.request.Request(url, headers={"User-Agent": "MCA-kDET-data-audit/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response, part.open("wb") as output:
            expected = response.headers.get("Content-Length")
            total = 0
            last = time.monotonic()
            while block := response.read(1024 * 1024):
                output.write(block)
                total += len(block)
                if time.monotonic() - last > 10:
                    print(f"  {total / 1024**2:.1f} MiB", flush=True)
                    last = time.monotonic()
        if expected is not None and total != int(expected):
            raise ValueError(f"Incomplete download: {name}")
        with part.open("rb") as stream:
            magic = stream.read(8)
        if name.endswith(".gz") and magic[:2] != b"\x1f\x8b":
            raise ValueError(f"Not a gzip file: {name}")
        if name.endswith(".h5"):
            import h5py
            # GEO's H5 is a MATLAB 7.3 file: its HDF5 signature follows
            # a MATLAB user block rather than occurring at byte zero.
            if not h5py.is_hdf5(part):
                raise ValueError(f"Not an HDF5 file: {name}")
        part.rename(path)
    return {"file": name, "url": url, "bytes": path.stat().st_size, "sha256": sha256(path)}


def sample_table(suffix: str) -> pd.DataFrame:
    frame = pd.read_csv(ROOT / "raw" / (PREFIX + suffix))
    if frame["sampleID"].isna().any() or frame["sampleID"].duplicated().any():
        raise ValueError(f"Missing or duplicate sampleID: {suffix}")
    return frame


def geo_sample_annotations() -> pd.DataFrame:
    """GEO source records explicitly identify animal, tumor, age and condition."""
    with gzip.open(ROOT / "raw" / SOFT_NAME, "rt", encoding="utf-8") as stream:
        text = stream.read()
    rows = []
    for record in text.split("^SAMPLE = ")[1:]:
        row = {"geo_accession": record.splitlines()[0].strip()}
        for line in record.splitlines():
            if line.startswith("!Sample_title = "):
                title = line.split(" = ", 1)[1]
                row["sampleID"] = title.split(": ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                key, value = line.split(" = ", 1)[1].split(": ", 1)
                row["geo_" + key] = value
        rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.sampleID.duplicated().any():
        raise ValueError("Duplicate sample IDs in GEO SOFT records")
    return frame


def join_geo(cells: pd.DataFrame, geo: pd.DataFrame) -> pd.DataFrame:
    """Resolve published ID formatting/relabeling without fuzzy matching.

    Some processed cell indices omit zero padding. In the STb2 batch, the
    author's mouseID supplies the renamed GEO animal prefix. Tumor, plate and
    cell-index suffixes must match exactly and the mapping must be one-to-one.
    """
    def key(value: str) -> str:
        return re.sub(r"_S0*(\d+)$", lambda match: "_S" + str(int(match[1])), value)

    cells = cells.copy()
    geo = geo.rename(columns={"sampleID": "geo_sampleID"}).copy()
    geo["join_key"] = geo.geo_sampleID.map(key)
    if geo.join_key.duplicated().any():
        raise ValueError("GEO identifiers collide after zero-padding normalization")
    cells["join_key"] = cells.sampleID.map(key)
    cells["geo_join_rule"] = "cell_index_zero_padding"
    known = set(geo.join_key)
    for index in cells.index[~cells.join_key.isin(known)]:
        sample = cells.at[index, "sampleID"]
        suffix = re.search(r"_(T\d+_P\d+(?:_[NO]\d+)?_S\d+)$", sample)
        if not suffix:
            raise ValueError(f"Unrecognized cell suffix: {sample}")
        animal = re.sub(r"_T\d+$", "", cells.at[index, "mouseID"])
        candidate = key(animal + "_" + suffix[1])
        if candidate not in known:
            raise ValueError(f"No GEO counterpart after author animal relabeling: {sample}")
        cells.at[index, "join_key"] = candidate
        cells.at[index, "geo_join_rule"] = "author_mouseID_relabel_and_zero_padding"
    if cells.join_key.duplicated().any():
        raise ValueError("Multiple processed cells map to the same GEO record")
    result = cells.merge(geo, on="join_key", validate="one_to_one", sort=False)
    if list(result.sampleID) != list(cells.sampleID):
        raise ValueError("GEO join changed cell order")
    return result.drop(columns="join_key")


def audit_expression(path: Path, sample_ids: list[str], gene_count: int) -> dict:
    import h5py

    datasets = {}
    with h5py.File(path, "r") as handle:
        def record(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets[name] = {"shape": list(obj.shape), "dtype": str(obj.dtype)}
        handle.visititems(record)
        # Export the actual schema before interpreting any sparse conventions.
        sample_matches = []
        for name in datasets:
            dataset = handle[name]
            if dataset.ndim == 1 and dataset.size == len(sample_ids) and dataset.dtype.kind in "OSU":
                values = [x.decode() if isinstance(x, bytes) else str(x) for x in dataset[:]]
                if set(values) == set(sample_ids):
                    sample_matches.append({"dataset": name, "same_order": values == sample_ids})
        sparse_audit = {}
        if set(handle) == {"i", "j", "v"}:
            if not (handle["i"].shape == handle["j"].shape == handle["v"].shape):
                raise ValueError("Different sparse coordinate/value lengths")
            for name, limit in [("i", gene_count), ("j", len(sample_ids))]:
                minimum, maximum = float("inf"), -float("inf")
                for start in range(0, handle[name].size, 1_000_000):
                    values = handle[name][0, start:start + 1_000_000]
                    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
                        raise ValueError(f"Noninteger sparse coordinates: {name}")
                    minimum = min(minimum, float(values.min()))
                    maximum = max(maximum, float(values.max()))
                if minimum < 1 or maximum > limit:
                    raise ValueError(f"Sparse coordinates outside metadata dimension: {name}")
                sparse_audit[name + "_range"] = [int(minimum), int(maximum)]
            minimum, maximum = float("inf"), -float("inf")
            for start in range(0, handle["v"].size, 1_000_000):
                values = handle["v"][0, start:start + 1_000_000]
                if not np.isfinite(values).all() or (values < 0).any():
                    raise ValueError("Nonfinite/negative normalized expression")
                minimum = min(minimum, float(values.min()))
                maximum = max(maximum, float(values.max()))
            sparse_audit.update({"index_origin": 1, "nonzero_entries": handle["v"].size,
                                 "matrix_shape_genes_by_cells": [gene_count, len(sample_ids)],
                                 "value_range": [minimum, maximum],
                                 "axis_labels": "Published companion geneTable/smpTable row order; H5 has no embedded IDs."})
    return {"datasets": datasets, "metadata_cells": len(sample_ids),
            "gene_table_rows": gene_count, "sample_id_matches": sample_matches,
            "sparse_coordinates": sparse_audit}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", action="store_true")
    args = parser.parse_args()
    (ROOT / "raw").mkdir(parents=True, exist_ok=True)
    output = ROOT / "prepared"
    output.mkdir(exist_ok=True)
    manifest_path = ROOT / "raw" / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    for suffix in METADATA + ([EXPRESSION] if args.expression else []):
        entry = acquire(PREFIX + suffix)
        manifest[entry["file"]] = entry
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    entry = acquire(SOFT_NAME, SOFT_URL)
    manifest[SOFT_NAME] = entry
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    samples = sample_table("smpTable.csv.gz")
    annot = sample_table("dZ_annot_smpTable.csv.gz")
    qc = sample_table("dZ_QCstat_smpTable.csv.gz")
    for name, other in [("annotations", annot), ("QC", qc)]:
        if set(other.sampleID) != set(samples.sampleID):
            raise ValueError(f"Different sample IDs in {name}")
    cells = samples.merge(annot, on="sampleID", validate="one_to_one", sort=False)
    cells = cells.merge(qc, on="sampleID", validate="one_to_one", sort=False)
    geo = geo_sample_annotations()
    cells = join_geo(cells, geo)
    if list(cells.sampleID) != list(samples.sampleID):
        raise ValueError("Sample order changed in metadata joins")

    # Conservative leakage guard, not an assertion about verified animal identity.
    # Euthanized mice cannot recur at later time points. Keep the full condition
    # prefix and only remove a terminal tumor suffix within that prefix.
    cells["donor_block"] = cells.mouseID.str.replace(r"_T\d+$", "", regex=True)
    cells["tumor_suffix"] = cells.mouseID.str.extract(r"(_T\d+)$", expand=False).fillna("")
    cells["source_condition"] = cells.typeID
    cells["author_stage"] = cells.timesimple
    # Use the actual condition and animal recorded in GEO, not a bare m1/m2
    # identifier reused across different sacrificed cohorts.
    condition = cells["geo_timepoint"].str.replace(r"_[mf]\d+$", "", regex=True)
    cells["verified_animal_block"] = (condition + "_" +
                                       cells["geo_biological_replicate_animal"])
    verified_counts = cells.groupby(["author_stage", "source_condition", "geo_timepoint",
                                     "geo_sack_timepoint (weeks)", "geo_strain"], dropna=False).agg(
        cells=("sampleID", "size"), animals=("verified_animal_block", "nunique"),
        source_units=("mouseID", "nunique")
    ).reset_index()
    verified_counts.to_csv(output / "animal_counts.csv", index=False, encoding="utf-8-sig")
    cells.to_csv(output / "cell_metadata.csv", index=False, encoding="utf-8-sig")

    unit_counts = cells.groupby(["author_stage", "source_condition", "donor_block", "mouseID"], sort=True).agg(
        cells=("sampleID", "size"), plates=("plateID", "nunique")
    ).reset_index()
    unit_counts.to_csv(output / "source_units.csv", index=False, encoding="utf-8-sig")
    summary = cells.groupby(["author_stage", "source_condition"], sort=True).agg(
        cells=("sampleID", "size"), source_units=("mouseID", "nunique"),
        conservative_donor_blocks=("donor_block", "nunique"), plates=("plateID", "nunique")
    ).reset_index()
    summary.to_csv(output / "stage_counts.csv", index=False, encoding="utf-8-sig")
    genes = pd.read_csv(ROOT / "raw" / (PREFIX + "geneTable.csv.gz"))

    audit = {
        "series": "GSE154989", "cells": len(cells), "gene_table_rows": len(genes),
        "gene_table_columns": list(genes), "author_mouseID_values": cells.mouseID.nunique(),
        "conservative_donor_blocks": cells.donor_block.nunique(),
        "geo_animal_blocks": cells.verified_animal_block.nunique(),
        "geo_join_rules": cells.geo_join_rule.value_counts().to_dict(),
        "identity_status": "GEO timepoint + biological_replicate_animal; repeated tumor labels stay in the same animal block.",
        "type_time_discrepancies": summary.loc[summary.source_condition.eq("07_KP_18w_ND") |
                                                 summary.source_condition.eq("07_KP_20w_STb2")].to_dict("records"),
        "inferential_tests_run": False,
    }
    expression = ROOT / "raw" / (PREFIX + EXPRESSION)
    if expression.exists():
        audit["expression"] = audit_expression(expression, list(samples.sampleID), len(genes))
    (output / "data_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(summary.to_string(index=False), flush=True)
    print(verified_counts.to_string(index=False), flush=True)
    print(json.dumps(audit, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
