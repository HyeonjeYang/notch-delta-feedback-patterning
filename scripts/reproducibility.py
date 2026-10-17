"""Record the environment and content hashes without promising bitwise ODE output."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    paths = [ROOT / "requirements.txt", ROOT / "requirements-lock.txt"]
    for directory in ("model", "control", "geometry", "simulation", "stability", "scripts", "figures", "tests"):
        paths.extend((ROOT / directory).glob("*.py"))
    return {p.relative_to(ROOT).as_posix(): sha256(p) for p in sorted(paths) if p.is_file()}


def new_manifest(steps: list[str]) -> dict:
    import numpy as np
    packages = {}
    for line in (ROOT / "requirements-lock.txt").read_text(encoding="utf-8").splitlines():
        if "==" in line:
            name, expected = line.split("==", 1)
            try:
                actual = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                actual = None
            packages[name] = {"expected": expected, "actual": actual}
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    config = getattr(np.__config__, "CONFIG", {})
    blas = config.get("Build Dependencies", {}).get("blas", {})
    return {
        "started_utc": datetime.now(timezone.utc).isoformat(), "status": "running", "steps": steps,
        "python": platform.python_version(), "python_implementation": platform.python_implementation(),
        "platform": platform.platform(), "machine": platform.machine(),
        "git_head": git.stdout.strip() if git.returncode == 0 else None,
        "source_sha256_at_start": source_hashes(), "packages": packages,
        "blas": {k: blas.get(k) for k in ("name", "version", "openblas configuration")},
        "thread_environment": {k: os.environ.get(k) for k in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "frozen_global_parameters_sha256": sha256(ROOT / "figures/data/global_sample_parameters.csv"),
        "random_protocols": {
            "global": {"seed": 20260817, "n": 4000, "sampler": "scipy.stats.qmc.LatinHypercube",
                       "api": "legacy seed=integer", "scramble": True, "strength": 1, "optimization": None},
            "mode_search": {"seed": 20261016, "generator": "Generator(PCG64(seed))",
                            "n_per_box": 20000, "box_order": [1, 10], "reset_between_boxes": False,
                            "sampling": "independent draws; rates log-uniform, Hill exponents uniform"},
            "nonlinear_reference": "deterministic eigenvector initial state; no random seed",
            "initial_conditions": {"seeds": list(range(2026101901,2026101909)),
                "generator": "Generator(PCG64(seed))", "deterministic_after_initialisation": True,
                "protocol": "figures/data/initial_conditions/protocol.json",
                "protocol_sha256": sha256(ROOT / "figures/data/initial_conditions/protocol.json")
                    if (ROOT / "figures/data/initial_conditions/protocol.json").exists() else None},
            "reviewer_extensions": {"protocol": "figures/data/reviewer_extensions/protocol.json",
                "geometry_seeds": [2026102001,2026102002,2026102003],
                "heterogeneity_stream": "PCG64(seed+1)",
                "pattern_recovery_seeds": list(range(2026102201,2026102213)),
                "protocol_sha256": sha256(ROOT / "figures/data/reviewer_extensions/protocol.json")
                    if (ROOT / "figures/data/reviewer_extensions/protocol.json").exists() else None}},
        "reproduction_contract": "Frozen inputs and numerical conclusions; floating-point outputs and PDF metadata need not be byte-identical."
    }


def write_manifest(manifest: dict, name: str) -> None:
    target = ROOT / "figures/data" / name
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(target)


def finish_manifest(manifest: dict, status: str) -> None:
    manifest["status"] = status
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["source_changed_during_run"] = source_hashes() != manifest["source_sha256_at_start"]
    manifest["artifact_sha256_at_end"] = {
        p.relative_to(ROOT).as_posix(): sha256(p)
        for folder in ("data", "main", "supplementary")
        for p in sorted((ROOT / "figures" / folder).rglob("*"))
        if p.is_file() and "manifest" not in p.name
    }
