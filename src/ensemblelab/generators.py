"""Conformer-ensemble generation backed by RDKit ETKDG."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from ase import Atoms
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


@dataclass(slots=True)
class Conformer:
    """One generated conformer and its conformer-level data.

    ``energy`` is ``None`` until an optimization backend assigns a value. It
    belongs to the conformer rather than a parallel ensemble-level array so a
    conformer's geometry and computed properties remain together.
    """

    id: int
    atoms: Atoms
    energy: float | None = None
    energy_unit: str | None = None
    optimization_method: str | None = None
    optimization_converged: bool | None = None

    def show(self) -> str:
        """Return a concise, human-readable conformer summary.

        This inspection method only renders stored conformer data; it does not
        perform an optimization or any other calculation.
        """
        from .display.summaries import conformer_summary

        return conformer_summary(self)
# generation history generator helper
def _generation_history(
    smiles: str,
    molecule: Chem.Mol,
    n_requested: int,
    n_generated: int,
) -> dict[str, Any]:
    """Create the initial provenance record for an ensemble."""

    return {
        "process": "generation",
        "method": "ETKDGv3",
        "requested_smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(
            Chem.RemoveHs(molecule),
            canonical=True,
        ),
        "n_requested": n_requested,
        "n_generated": n_generated,
        "random_seed": 42,
    }

@dataclass(slots=True)
class Ensemble:
    """A molecular conformational ensemble with aligned conformer-level data.

    Generation initializes structures and provenance only. Every conformer has
    ``energy=None`` until optimization assigns a physical energy.
    """

    smiles: str
    molecule: Chem.Mol
    conformers: tuple[Conformer, ...]
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        ids = [conformer.id for conformer in self.conformers]
        if len(ids) != len(set(ids)):
            raise ValueError("conformer IDs must be unique within an ensemble.")

    @classmethod
    def from_smiles(cls, smiles: str, n_confs: int = 25) -> Ensemble:
        """Create an unoptimized ensemble from SMILES.

        This is the object oriented equivalent of "generate" function and keeps
        the public workflow ready for subsequent "optimize" and "cluster"
        methods.
        """
        return generate(smiles, n_confs=n_confs)

    @property
    def conformer_ids(self) -> tuple[int, ...]:
        """RDKit conformer IDs in the alignment order used by ensemble data."""
        return tuple(conformer.id for conformer in self.conformers)

    def rdkit_conformer(self, conformer_id: int) -> Chem.Conformer:
        """Return an RDKit conformer by its stable ID."""
        if conformer_id not in self.conformer_ids:
            raise KeyError(f"Unknown conformer ID: {conformer_id}")
        return self.molecule.GetConformer(conformer_id)

    def show(
        self,
        *,
        history: bool = False,
        metadata: bool = False,
        conformers: bool = True,
    ) -> str:
        """Return a human-readable view of this ensemble.

        The default view contains the ensemble summary and conformer table.
        Set ``history`` or ``metadata`` to include the corresponding stored
        provenance sections. Rendering is read-only and never triggers
        expensive analysis or optimization.
        """
        from .display.summaries import ensemble_history, ensemble_metadata, ensemble_summary

        sections = [ensemble_summary(self, include_conformers=conformers)]
        if history:
            sections.append(ensemble_history(self))
        if metadata:
            sections.append(ensemble_metadata(self))
        return "\n\n".join(sections)


def generate(smiles: str, n_confs: int = 25) -> Ensemble:
    """Embed an unoptimized conformational ensemble from a SMILES string.

    Uses the notebook's RDKit ETKDGv3 setup (explicit seed 42) on a
    hydrogenated molecule. Optimization and energies intentionally belong to
    the next module.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty SMILES string.")
    if isinstance(n_confs, bool) or not isinstance(n_confs, int) or n_confs < 1:
        raise ValueError("n_confs must be a positive integer.")

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")
    molecule = Chem.AddHs(molecule)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    conformer_ids = tuple(
        int(identifier)
        for identifier in AllChem.EmbedMultipleConfs(
            molecule, numConfs=n_confs, params=params
        )
    )
    if not conformer_ids:
        raise ValueError("RDKit could not embed any conformers for this molecule.")

    symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
    conformers = tuple(
        Conformer(
            id=conformer_id,
            atoms=Atoms(
                symbols=symbols,
                positions=np.asarray(
                    molecule.GetConformer(conformer_id).GetPositions(), dtype=float
                ).copy(),
            ),
        )
        for conformer_id in conformer_ids
    )
    metadata = {
    "n_conformers": len(conformers),
    "optimization_status": "unoptimized",
    "energy_status": "uncomputed",
    "energy_unit": None,
    "rdkit_version": rdBase.rdkitVersion,
    "processing_history": [
        _generation_history(
            smiles,
            molecule,
            n_confs,
            len(conformers),
        )
    ],}
    return Ensemble(
        smiles=smiles,
        molecule=molecule,
        conformers=conformers,
        metadata=metadata,
    )
