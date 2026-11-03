"""Freeze, prepare and execute the single-cell study; no result-dependent choices."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
RSCRIPT = Path("D:/R-4.6.1/bin/Rscript.exe")
PRIMARY_SEED = 20261027


def specifications() -> list[dict]:
    common = dict(seed=PRIMARY_SEED, cells=30, dimension=100, cohort="all", role="sensitivity")
    specs = [dict(common, name="primary", role="primary"),
             dict(common, name="genes_50", dimension=50),
             dict(common, name="genes_200", dimension=200),
             dict(common, name="cells_20", cells=20),
             dict(common, name="cells_40", cells=40),
             dict(common, name="middle_20w_only", cohort="middle_20w_only"),
             dict(common, name="one_tumor", cohort="single_tumor")]
    specs.extend(dict(common, name=f"split_{i:02d}", seed=PRIMARY_SEED + i,
                      role="split_stability") for i in range(1, 11))
    return specs


def prepared_path(spec: dict) -> Path:
    if spec["name"] in {"primary", "genes_50", "genes_200"}:
        return ROOT / "prepared" / "analysis"
    return ROOT / "prepared" / spec["name"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    specs = specifications()
    plan = dict(version="KP_stage_execution_v1", primary_seed=PRIMARY_SEED,
                alpha=0.05, specifications=specs,
                main_methods=["MCA", "DISCO", "MMD", "MST"],
                library=["l1", "quad", "knn", "rf"],
                calibration="Complete enumeration of whole-test-animal labels with fixed animal group counts; observational block exchangeability assumed.",
                inference_scope="Joint distributions of sampled animal cell blocks, not iid-cell marginal-distribution theorem.",
                sensitivity_policy="Report every prespecified configuration and split; no min-p selection or power claim.",
                one_tumor_rule="For each 30-week mouse choose one tumor uniformly before cell sampling using seed+179.",
                middle_20w_rule="Omit 18-week mice; middle group 1 fitting and 2 testing mice, all other stage counts unchanged.")
    plan_path = ROOT / "study_plan.json"
    if plan_path.exists() and json.loads(plan_path.read_text(encoding="utf-8")) != plan:
        raise ValueError("Study plan differs from the frozen plan; make a versioned plan")
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    for key in ("LC_ALL", "LC_CTYPE", "LANG"):
        env.pop(key, None)
    env["R_LIBS_USER"] = "C:/Users/zhang/AppData/Local/R/win-library/4.6"
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    logs = ROOT / "results" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    prepared = set()
    for spec in specs:
        destination = prepared_path(spec)
        if destination in prepared:
            continue
        prepared.add(destination)
        if destination == ROOT / "prepared" / "analysis":
            if not (destination / "design.json").exists():
                raise FileNotFoundError("Prespecified primary input missing")
            continue
        command = [sys.executable, str(ROOT / "prepare_analysis.py"),
                   f"--seed={spec['seed']}", f"--cells-per-mouse={spec['cells']}",
                   f"--output-dir={destination}"]
        if spec["cohort"] == "middle_20w_only":
            command.append("--middle-20w-only")
        if spec["cohort"] == "single_tumor":
            command.append("--single-tumor")
        # Preserve existing prepared inputs; fail on unexpected parameters.
        if (destination / "design.json").exists():
            design = json.loads((destination / "design.json").read_text(encoding="utf-8"))
            if design["seed"] != spec["seed"] or design["cells_per_mouse"] != spec["cells"]:
                raise ValueError(f"Prepared parameters differ: {destination}")
            continue
        result = subprocess.run(command, env=env, capture_output=True, timeout=300)
        (logs / f"{spec['name']}_prepare.log").write_bytes(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Preparation failed: {spec['name']}; inspect log")
        print(f"Prepared {spec['name']}", flush=True)

    def run(spec: dict) -> str:
        command = [str(RSCRIPT), str(ROOT / "run_analysis.R"),
                   f"--input={prepared_path(spec)}", f"--output={ROOT / 'results' / spec['name']}",
                   f"--dimension={spec['dimension']}", f"--seed={spec['seed']}"]
        # Transient R startup failures are known on this host. Same seed each attempt.
        for attempt in range(1, 4):
            result = subprocess.run(command, env=env, capture_output=True, timeout=1800)
            (logs / f"{spec['name']}_R_attempt{attempt}.log").write_bytes(result.stdout + result.stderr)
            if result.returncode == 0 and (ROOT / "results" / spec["name"] / "analysis.rds").exists():
                export = subprocess.run([str(RSCRIPT), str(ROOT / "export_exact_scores.R"),
                                         str(ROOT / "results" / spec["name"])], env=env, capture_output=True, timeout=60)
                if export.returncode:
                    raise RuntimeError(f"Full-precision score export failed: {spec['name']}")
                return spec["name"]
        raise RuntimeError(f"R analysis failed: {spec['name']}; inspect logs")

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        tasks = [pool.submit(run, spec) for spec in specs]
        for done, future in enumerate(as_completed(tasks), 1):
            print(f"Completed {done}/{len(specs)}: {future.result()}", flush=True)
    provenance = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in [ROOT / name for name in ("study_plan.json", "prepare_analysis.py", "run_analysis.R", "run_study.py", "export_exact_scores.R")]}
    (ROOT / "results" / "execution_manifest.json").write_text(
        json.dumps(dict(completed=len(specs), sha256=provenance), indent=2) + "\n", encoding="utf-8")
    print("All planned configurations completed.", flush=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
