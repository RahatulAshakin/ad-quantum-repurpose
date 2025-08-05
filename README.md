# AD Quantum Repurposing (Proposal 1)

Baseline docking of 10 CNS-approved drugs vs BACE1 (PDB 4DVF) with smina; quick BBB heuristics; composite ranking.  
Plus a minimal Qiskit Nature VQE demo.

## Repo structure
- `data/raw/` PDB (4DVF) and seed SMILES
- `data/processed/` prepared receptor & ligands
- `results/docking/` docking outputs + `docking_results.csv`
- `results/reports/` `ranking.md`, `ranking_table.csv`, `docking_bar.png`, `vqe_demo.md`
- `results/quantum/` `vqe_h2.json`
- `scripts/` helper scripts (incl. `06_vqe_workflow.py`)
