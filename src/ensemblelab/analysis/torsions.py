import numpy as np
from ase import Atoms
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from ensemblelab.generators import Conformer, Ensemble


def find_heavy_atom_torsions(mol: Ensemble) -> dict:
    mol = mol.get_mol()
    torsions = {}

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        if bond.IsInRing():
            continue

        atom_2 = bond.GetBeginAtom()
        atom_3 = bond.GetEndAtom()

        if atom_2.GetAtomicNum() == 1 or atom_3.GetAtomicNum() == 1:
            continue

        left_neighbors = [
            atom.GetIdx()
            for atom in atom_2.GetNeighbors()
            if atom.GetIdx() != atom_3.GetIdx() and atom.GetAtomicNum() > 1
        ]

        right_neighbors = [
            atom.GetIdx()
            for atom in atom_3.GetNeighbors()
            if atom.GetIdx() != atom_2.GetIdx() and atom.GetAtomicNum() > 1
        ]

        if not left_neighbors or not right_neighbors:
            continue

        atom_1 = left_neighbors[0]
        atom_4 = right_neighbors[0]

        name = f"torsion_{len(torsions) + 1}_{atom_1}_{atom_2.GetIdx()}_{atom_3.GetIdx()}_{atom_4}"
        torsions[name] = (atom_1, atom_2.GetIdx(), atom_3.GetIdx(), atom_4)

    return torsions


"""manually override this later with chemically meaningful names:
torsion_definitions = {
    "backbone_1": (0, 1, 2, 3),
    "backbone_2": (1, 2, 3, 4),
}

if "torsion_definitions" not in globals() or not torsion_definitions:
    torsion_definitions = find_heavy_atom_torsions(m    ol)"""



