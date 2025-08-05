# 06_vqe_workflow.py
# Minimal VQE demo compatible with Qiskit 1.4 / Nature 0.7.2
# Computes H2 ground-state energy at 0.735 Å (sto-3g) and writes JSON/MD

from pathlib import Path
import json
from datetime import datetime

import numpy as np

from qiskit.primitives import Estimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.algorithms import GroundStateEigensolver

# --- paths
PROJECT = Path("/mnt/c/Users/ashak/ad-quantum-repurpose").resolve()
OUTDIR = PROJECT / "results" / "quantum"
OUTDIR.mkdir(parents=True, exist_ok=True)

# --- define tiny test system (no Molecule API needed in v0.7)
driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",  # Å
    basis="sto-3g",
    charge=0,
    spin=0,
)

problem = driver.run()
num_spatial_orbitals = problem.num_spatial_orbitals
num_particles = problem.num_particles

# --- mapping & ansatz
mapper = JordanWignerMapper()
hf = HartreeFock(num_spatial_orbitals, num_particles, mapper)
ansatz = UCCSD(
    num_spatial_orbitals=num_spatial_orbitals,
    num_particles=num_particles,
    qubit_mapper=mapper,
    initial_state=hf,
)

# Start at HF point (zeros); avoids random starts
initial_point = np.zeros(ansatz.num_parameters)

# --- classical optimizer & quantum primitive
optimizer = SLSQP(maxiter=200)
estimator = Estimator()

vqe = VQE(estimator, ansatz=ansatz, optimizer=optimizer)
vqe.initial_point = initial_point

# --- ground state solver
solver = GroundStateEigensolver(mapper, vqe)
result = solver.solve(problem)

energy_ha = float(result.total_energies[0].real)
energy_ev = energy_ha * 27.211386245988  # Hartree -> eV

payload = {
    "system": "H2",
    "bond_length_A": 0.735,
    "basis": "sto-3g",
    "mapper": "Jordan-Wigner",
    "ansatz": "UCCSD",
    "optimizer": "SLSQP(maxiter=200)",
    "energy_hartree": energy_ha,
    "energy_eV": energy_ev,
    "timestamp": datetime.utcnow().isoformat() + "Z",
    "qiskit_versions": {
        "note": "verified in adquant_q: qiskit~=1.4, qiskit-nature~=0.7.2",
    },
}

# --- write artifacts
json_path = OUTDIR / "vqe_h2.json"
json_path.write_text(json.dumps(payload, indent=2))

md_path = PROJECT / "results" / "reports" / "vqe_demo.md"
md_path.parent.mkdir(parents=True, exist_ok=True)
md_path.write_text(
    "# Quantum demo (VQE) — H2\n"
    f"- Bond length: 0.735 Å  \n"
    f"- Basis: sto-3g  \n"
    f"- Mapper: Jordan–Wigner; Ansatz: UCCSD; Optimizer: SLSQP  \n"
    f"- **Ground-state energy**: {energy_ha:.8f} Ha ({energy_ev:.6f} eV)\n"
)

print("WROTE:", json_path)
print("WROTE:", md_path)
print(f"Energy: {energy_ha:.8f} Ha ({energy_ev:.6f} eV)")

