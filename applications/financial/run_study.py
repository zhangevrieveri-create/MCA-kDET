"""Run all frozen financial configurations with deterministic, verified checkpoints."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
R = "D:/R-4.6.1/bin/Rscript.exe"


def main():
    plan = json.loads((ROOT / "study_plan.json").read_text(encoding="utf-8"))
    env = os.environ.copy()
    for key in ("LC_ALL", "LC_CTYPE", "LANG"):
        env.pop(key, None)
    env["R_LIBS_USER"] = "C:/Users/zhang/AppData/Local/R/win-library/4.6"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[key] = "1"
    logs = ROOT / "results" / "logs"; logs.mkdir(parents=True, exist_ok=True)

    def run(spec):
        for attempt in range(1, 4):
            completed = subprocess.run([R, str(ROOT / "run_analysis.R"), f"--config={spec['name']}"],
                                       env=env, capture_output=True, timeout=1800)
            (logs / f"{spec['name']}_attempt{attempt}.log").write_bytes(completed.stdout + completed.stderr)
            if completed.returncode == 0 and (ROOT / "results" / spec["name"] / "analysis.rds").exists():
                return spec["name"]
        raise RuntimeError(f"Failed {spec['name']}; inspect the log")

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run, spec) for spec in plan["specifications"]]
        for i, completed in enumerate(as_completed(futures), 1):
            print(f"Completed {i}/{len(futures)}: {completed.result()}", flush=True)
    files = [ROOT / name for name in ("study_plan.json", "prepare_data.py", "prepare_filter.R", "run_analysis.R", "run_study.py")]
    (ROOT / "results" / "execution_manifest.json").write_text(json.dumps(dict(
        completed=len(plan["specifications"]), sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
