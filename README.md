# AD Quantum Repurposing (Proposal 1 Baseline)

**Goal:** classical docking + BBB heuristics + quantum demo to scaffold a quantum-enhanced repurposing pipeline for Alzheimer’s disease (AD).

## What’s in this repo
- `data/raw/`  – PDB **4DVF** and seed SMILES  
- `data/processed/` – prepared receptor & ligands (PDBQT/SDF)  
- `results/docking/` – docking logs/poses + `docking_results.csv`  
- `results/reports/` – `ranking.md`, `ranking_table.csv`, `docking_bar.png`, `vqe_demo.md`  
- `results/quantum/` – `vqe_h2.json`  
- `scripts/` – helper scripts (incl. `06_vqe_workflow.py`)  
- `envs/wsl_adquant_q.yml` – WSL environment for Qiskit/PySCF

## Methods (summary)
- **Protein prep:** PDBFixer/OpenMM on 4DVF; box center `[3.061, −0.128, 17.532]`, size `[22,22,22]` Å.
- **Ligands:** 10 CNS-approved drugs; RDKit for 3D + MW/PSA/cLogP; OpenBabel for PDBQT & Gasteiger charges.
- **Docking:** `smina` (Vina score); parsed best-pose energies; composite rank = binding + BBB features.
- **Quantum demo:** Qiskit Nature VQE (H₂) to validate the quantum stack; outputs JSON + report.

## Results (top-10)
| Rank | Ligand | Affinity (kcal/mol) | MW | PSA | cLogP | BBB | Composite |
|---:|---|---:|---:|---:|---:|:---:|---:|
| 1 | **donepezil** | **−8.2** | 392.21 | 42.01 | 3.856 | ✅ | **2.2470** |
| 2 | **galantamine** | **−7.6** | 278.09 | 39.44 | 4.245 | ✅ | **1.5164** |
| 3 | **nilvadipine** | **−7.0** | 433.05 | 78.51 | 4.464 | ✅ | **0.7858** |
| 4 | fluoxetine | −6.7 | 259.14 | 35.25 | 3.796 | ✅ | 0.4205 |
| 5 | sertraline | −6.6 | 283.09 | 12.03 | 4.467 | ✅ | 0.2987 |
| 6 | riluzole | −6.4 | 253.02 | 76.72 | 2.140 | ✅ | 0.0552 |
| 7 | rivastigmine | −6.1 | 221.11 | 38.77 | 2.036 | ✅ | −0.3101 |
| 8 | memantine | −5.9 | 181.18 | 26.02 | 3.084 | ✅ | −0.5537 |
| 9 | rasagiline | −5.7 | 191.11 | 12.03 | 2.175 | ✅ | −0.7972 |
| 10 | selegiline | −5.4 | 173.12 | 12.03 | 2.036 | ✅ | −1.1625 |

**Highlights**
- 10/10 ligands docked and parsed (**100%** success).  
- Best vs worst |ΔE| improvement: **~52%** (donepezil vs selegiline).  
- Best vs median |ΔE| improvement: **~26%** (donepezil vs median).  
- All pass BBB heuristics used here (CNS-plausible).

See `results/reports/ranking.md` and `results/reports/docking_bar.png`.

## Quantum check
VQE demo for H₂ succeeded (see `results/quantum/vqe_h2.json`, `results/reports/vqe_demo.md`). This validates Qiskit/PySCF on WSL and is a placeholder for target-specific quantum simulations we’ll add next.

## Limitations
Rigid receptor; approximate scoring; single target; simple composite score; no experimental validation yet.

## Next steps
- Rescoring (RFScore/MMGBSA), interaction analysis, ensemble docking.  
- Add a second AD target (e.g., tau) for polypharmacology.  
- Swap demo H₂ for ligand-fragment **active-space VQE**.  
- Plan organoid/iPSC validation with ≥30% reduction thresholds for Aβ/p-tau.

