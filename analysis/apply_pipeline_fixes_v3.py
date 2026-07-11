#!/usr/bin/env python3
from pathlib import Path

path = Path("analysis/corrected_resubmission_pipeline.py")
text = path.read_text(encoding="utf-8")

text = text.replace(
    "from rdkit.Chem import AllChem, ChemicalFeatures, Crippen, Descriptors, Lipinski, rdMolDescriptors",
    "from rdkit.Chem import AllChem, ChemicalFeatures, Crippen, Descriptors, Lipinski, rdMolDescriptors, rdFMCS",
)

old = '''def direct_symmetry_aware_rmsd(reference: Chem.Mol, pose: Chem.Mol) -> tuple[float, list[int]]:
    ref = Chem.RemoveHs(reference, sanitize=False)
    mob = Chem.RemoveHs(pose, sanitize=False)
    if ref.GetNumAtoms() != mob.GetNumAtoms():
        raise RuntimeError(f"RMSD atom-count mismatch: reference={ref.GetNumAtoms()}, pose={mob.GetNumAtoms()}")
    ref_conf = ref.GetConformer()
    mob_conf = mob.GetConformer()
    matches = mob.GetSubstructMatches(ref, uniquify=False, useChirality=False, maxMatches=20000)
    if not matches:
        # Fallback to order only when element sequences agree.
        if [a.GetSymbol() for a in ref.GetAtoms()] != [a.GetSymbol() for a in mob.GetAtoms()]:
            raise RuntimeError("Could not establish an atom mapping for direct RMSD")
        matches = [tuple(range(ref.GetNumAtoms()))]
    best = (float("inf"), [])
    for match in matches:
        diffs = []
        for ref_idx, mob_idx in enumerate(match):
            rp = ref_conf.GetAtomPosition(ref_idx)
            mp = mob_conf.GetAtomPosition(mob_idx)
            diffs.append((rp.x-mp.x)**2 + (rp.y-mp.y)**2 + (rp.z-mp.z)**2)
        rmsd = math.sqrt(float(np.mean(diffs)))
        if rmsd < best[0]:
            best = (rmsd, list(match))
    return best
'''

new = '''def largest_covalent_fragment(mol: Chem.Mol) -> Chem.Mol:
    """Return the largest heavy-atom fragment without inventing connectivity.

    RCSB ModelServer may include crystallographic waters in the ligand SDF. Those
    disconnected water fragments must not contribute to the docking center or RMSD.
    """
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not fragments:
        raise RuntimeError("No fragments found in native-ligand structure")
    return max(fragments, key=lambda m: (m.GetNumHeavyAtoms(), m.GetNumAtoms()))


def direct_symmetry_aware_rmsd(reference: Chem.Mol, pose: Chem.Mol) -> tuple[float, list[int]]:
    """Direct (no superposition) heavy-atom RMSD with graph/symmetry mapping.

    PDBQT-to-SDF conversion can reorder atoms and alter bond orders. A full-molecule
    maximum common substructure therefore maps atoms by element and connectivity while
    allowing bond-order differences. All full-atom mappings are evaluated in the fixed
    receptor coordinate frame; the smallest direct RMSD is reported.
    """
    ref = Chem.RemoveHs(reference, sanitize=False)
    mob = Chem.RemoveHs(pose, sanitize=False)
    if ref.GetNumAtoms() != mob.GetNumAtoms():
        raise RuntimeError(f"RMSD atom-count mismatch: reference={ref.GetNumAtoms()}, pose={mob.GetNumAtoms()}")
    for mol in (ref, mob):
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            Chem.GetSymmSSSR(mol)

    mcs = rdFMCS.FindMCS(
        [ref, mob],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        matchValences=False,
        timeout=120,
    )
    if mcs.canceled or mcs.numAtoms != ref.GetNumAtoms():
        raise RuntimeError(
            f"Could not establish a full heavy-atom mapping for RMSD: "
            f"MCS={mcs.numAtoms}/{ref.GetNumAtoms()} atoms"
        )
    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        raise RuntimeError("Could not parse the full-molecule MCS query for RMSD")
    ref_matches = ref.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=2000)
    mob_matches = mob.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=2000)
    if not ref_matches or not mob_matches:
        raise RuntimeError("The full-molecule MCS did not yield RMSD atom mappings")

    ref_conf = ref.GetConformer()
    mob_conf = mob.GetConformer()
    best_rmsd = float("inf")
    best_mapping: list[int] = []
    for ref_match in ref_matches:
        for mob_match in mob_matches:
            sq = []
            mapping = [-1] * ref.GetNumAtoms()
            for ref_idx, mob_idx in zip(ref_match, mob_match):
                rp = ref_conf.GetAtomPosition(ref_idx)
                mp = mob_conf.GetAtomPosition(mob_idx)
                sq.append((rp.x-mp.x)**2 + (rp.y-mp.y)**2 + (rp.z-mp.z)**2)
                mapping[ref_idx] = mob_idx
            rmsd = math.sqrt(float(np.mean(sq)))
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_mapping = mapping
    if any(i < 0 for i in best_mapping):
        raise RuntimeError("Incomplete full-molecule RMSD atom mapping")
    return best_rmsd, best_mapping
'''

if old not in text:
    raise RuntimeError("Expected RMSD function block was not found")
text = text.replace(old, new)

old_native = '''    native_crystal_sdf = get_native_sdf(raw_pdb, native_chain, work)
    native_crystal_mol = load_first_sdf(native_crystal_sdf, remove_hs=False, sanitize=False)
    native_prepared_sdf, native_pdbqt = prepare_native_ligand(native_crystal_sdf, work)
    native_heavy_idx, native_xyz, _ = mol_heavy_coords(native_crystal_mol)
'''
new_native = '''    native_crystal_sdf_raw = get_native_sdf(raw_pdb, native_chain, work)
    native_crystal_all = load_first_sdf(native_crystal_sdf_raw, remove_hs=False, sanitize=False)
    native_crystal_mol = largest_covalent_fragment(native_crystal_all)
    native_crystal_sdf = work / "4DVF_native_inhibitor_crystal_largest_fragment.sdf"
    writer = Chem.SDWriter(str(native_crystal_sdf))
    writer.write(native_crystal_mol)
    writer.close()
    native_prepared_sdf, native_pdbqt = prepare_native_ligand(native_crystal_sdf, work)
    native_heavy_idx, native_xyz, _ = mol_heavy_coords(native_crystal_mol)
'''
if old_native not in text:
    raise RuntimeError("Expected native-ligand block was not found")
text = text.replace(old_native, new_native)

text = text.replace(
    '"grid_definition": "Centroid of all heavy atoms in the 4DVF native inhibitor paired with BACE1 chain A",',
    '"grid_definition": "Centroid of the largest covalently connected heavy-atom fragment of the 4DVF native inhibitor paired with BACE1 chain A; disconnected crystallographic waters excluded",',
)
text = text.replace(
    "The docking center was defined from the heavy-atom centroid of native inhibitor chain {native_chain}, which is the inhibitor copy paired with chain A.",
    "The docking center was defined from the heavy-atom centroid of the largest covalently connected fragment of native inhibitor chain {native_chain}, which is the inhibitor copy paired with chain A; three disconnected crystallographic waters returned in the ModelServer ligand SDF were excluded.",
)
text = text.replace(
    "Direct symmetry-aware heavy-atom RMSD was {rmsd:.2f} A",
    "Direct full-graph/symmetry-mapped heavy-atom RMSD was {rmsd:.2f} A",
)
text = text.replace(
    'f"Direct symmetry-aware heavy-atom RMSD = {rmsd:.2f} A ({status}); smina score = {score:.2f} kcal/mol",',
    'f"Direct full-graph/symmetry-mapped heavy-atom RMSD = {rmsd:.2f} A ({status}); smina score = {score:.2f} kcal/mol",',
)

path.write_text(text, encoding="utf-8")
print(f"Patched {path}")
