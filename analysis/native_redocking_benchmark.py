#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path.cwd()
PIPELINE_PATH = ROOT / "analysis" / "corrected_resubmission_pipeline.py"
spec = importlib.util.spec_from_file_location("corrected_pipeline", PIPELINE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not import {PIPELINE_PATH}")
p = importlib.util.module_from_spec(spec)
sys.modules["corrected_pipeline"] = p
spec.loader.exec_module(p)

OUT = ROOT / "results" / "native_redocking_benchmark"
WORK = OUT / "work"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

raw_pdb = ROOT / "data" / "raw" / "proteins" / "4DVF.pdb"
if not raw_pdb.exists():
    raw_pdb.write_text(p.request_text("https://files.rcsb.org/download/4DVF.pdb"), encoding="utf-8")

chain_a = WORK / "4DVF_chain_A_protein_only.pdb"
receptor_pdb = WORK / "4DVF_chain_A_prepared.pdb"
receptor_pdbqt = WORK / "4DVF_chain_A_prepared.pdbqt"
p.extract_chain_a_protein(raw_pdb, chain_a)
p.prepare_receptor(chain_a, receptor_pdb, receptor_pdbqt)

native_chain = p.choose_native_chain(raw_pdb)
raw_native_sdf = p.get_native_sdf(raw_pdb, native_chain, WORK)
all_native = p.load_first_sdf(raw_native_sdf, remove_hs=False, sanitize=False)
native = p.largest_covalent_fragment(all_native)
native_sdf = WORK / "4DVF_native_largest_fragment.sdf"
w = Chem.SDWriter(str(native_sdf))
w.write(native)
w.close()
_, native_pdbqt = p.prepare_native_ligand(native_sdf, WORK)
_, xyz, _ = p.mol_heavy_coords(native)
center = xyz.mean(axis=0)
span = xyz.max(axis=0) - xyz.min(axis=0)

configs = [
    {"name": "current_22A_ex8", "box": 22.0, "exhaustiveness": 8, "num_modes": 20},
    {"name": "extent_plus4_28A_ex32", "box": 28.0, "exhaustiveness": 32, "num_modes": 20},
    {"name": "extent_plus6_32A_ex64", "box": 32.0, "exhaustiveness": 64, "num_modes": 20},
]


def parse_scores(text: str) -> list[float]:
    vals: list[float] = []
    for line in text.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+", line)
        if m:
            vals.append(float(m.group(2)))
    return vals


rows = []
for cfg in configs:
    stem = cfg["name"]
    pose_pdbqt = WORK / f"{stem}.pdbqt"
    log_file = WORK / f"{stem}.log"
    cmd = [
        "smina", "-r", str(receptor_pdbqt), "-l", str(native_pdbqt),
        "--center_x", f"{center[0]:.4f}", "--center_y", f"{center[1]:.4f}", "--center_z", f"{center[2]:.4f}",
        "--size_x", str(cfg["box"]), "--size_y", str(cfg["box"]), "--size_z", str(cfg["box"]),
        "--exhaustiveness", str(cfg["exhaustiveness"]), "--seed", "123",
        "--num_modes", str(cfg["num_modes"]), "--energy_range", "20",
        "--out", str(pose_pdbqt),
    ]
    proc = p.run(cmd, check=False, timeout=5400)
    log_file.write_text(proc.stdout or "", encoding="utf-8")
    if proc.returncode != 0:
        rows.append({"Configuration": stem, "Status": "smina failed", "Return_code": proc.returncode})
        continue
    scores = parse_scores(proc.stdout or "")
    pose_sdf = WORK / f"{stem}.sdf"
    p.run(["obabel", str(pose_pdbqt), "-O", str(pose_sdf)])
    poses = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=False) if m is not None]
    rmsds = []
    for i, mol in enumerate(poses):
        try:
            rmsd, _ = p.direct_symmetry_aware_rmsd(native, mol)
        except Exception as exc:
            print(f"RMSD failure {stem} mode {i+1}: {exc}", flush=True)
            rmsd = float("nan")
        rmsds.append(rmsd)
    if not poses:
        rows.append({"Configuration": stem, "Status": "No converted poses"})
        continue
    finite = [(i, x) for i, x in enumerate(rmsds) if math.isfinite(x)]
    min_i, min_rmsd = min(finite, key=lambda t: t[1]) if finite else (-1, float("nan"))
    rows.append({
        "Configuration": stem,
        "Status": "Completed",
        "Box_A": cfg["box"],
        "Exhaustiveness": cfg["exhaustiveness"],
        "Requested_modes": cfg["num_modes"],
        "Returned_modes": len(poses),
        "Top_score_kcal_mol": scores[0] if scores else float("nan"),
        "Top_pose_RMSD_A": rmsds[0] if rmsds else float("nan"),
        "Minimum_RMSD_A_among_modes": min_rmsd,
        "Minimum_RMSD_mode_rank": min_i + 1 if min_i >= 0 else "N/A",
        "Minimum_RMSD_mode_score_kcal_mol": scores[min_i] if min_i >= 0 and min_i < len(scores) else float("nan"),
    })

# Local refinement is diagnostic only and is not treated as self-docking validation.
local_pose = WORK / "local_only.pdbqt"
local_proc = p.run([
    "smina", "-r", str(receptor_pdbqt), "-l", str(native_pdbqt),
    "--local_only", "--seed", "123", "--out", str(local_pose),
], check=False, timeout=1800)
(WORK / "local_only.log").write_text(local_proc.stdout or "", encoding="utf-8")
if local_proc.returncode == 0 and local_pose.exists():
    local_sdf = WORK / "local_only.sdf"
    p.run(["obabel", str(local_pose), "-O", str(local_sdf)])
    lm = p.load_first_sdf(local_sdf, remove_hs=False, sanitize=False)
    local_rmsd, _ = p.direct_symmetry_aware_rmsd(native, lm)
    local_score = parse_scores(local_proc.stdout or "")
    rows.append({
        "Configuration": "local_only_diagnostic",
        "Status": "Completed; not global redocking",
        "Box_A": "N/A",
        "Exhaustiveness": "N/A",
        "Requested_modes": 1,
        "Returned_modes": 1,
        "Top_score_kcal_mol": local_score[0] if local_score else float("nan"),
        "Top_pose_RMSD_A": local_rmsd,
        "Minimum_RMSD_A_among_modes": local_rmsd,
        "Minimum_RMSD_mode_rank": 1,
        "Minimum_RMSD_mode_score_kcal_mol": local_score[0] if local_score else float("nan"),
    })

result = pd.DataFrame(rows)
result.to_csv(OUT / "Native_redocking_benchmark.csv", index=False)
metadata = {
    "native_chain": native_chain,
    "native_heavy_atom_count": native.GetNumHeavyAtoms(),
    "native_heavy_atom_span_A": [round(float(x), 3) for x in span],
    "grid_center_A": [round(float(x), 4) for x in center],
    "note": "Global top-pose RMSD is the prespecified validation criterion. Minimum RMSD among returned modes and local-only refinement are diagnostics, not substitutes for top-pose global self-docking.",
}
(OUT / "Benchmark_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(result.to_string(index=False), flush=True)
