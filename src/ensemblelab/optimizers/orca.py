"""ORCA conformer optimization through ASE."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from ase import Atoms
from ase.optimize import BFGS
from rdkit import Chem

from ..generators import Ensemble
from .base import BaseOptimizer

_EV_TO_KCAL_PER_MOL = 23.0605478306


@dataclass(frozen=True, slots=True)
class ORCAOptimizer(BaseOptimizer):
    """Optimize conformers using ASE's ORCA calculator."""

    orca_command: str | None = None
    orca_simple_input: str = "r2SCAN-3c"
    orca_blocks: str | None = None
    workdir: str | Path | None = None
    fmax: float = 0.05
    max_steps: int = 500
    charge: int | None = None
    multiplicity: int | None = None

    def __post_init__(self) -> None:
        self._validate_fmax(self.fmax)
        self._validate_max_steps(self.max_steps)

    @property
    def method(self) -> str:
        """Canonical optimization method name."""
        return "ORCA"

    def optimize(self, ensemble: Ensemble) -> Ensemble:
        """Optimize all conformers after resolving the ensemble's electronic state."""
        self._resolved_charge = self._resolve_charge(ensemble, self.charge)
        self._resolved_multiplicity = self._resolve_multiplicity(ensemble, self.multiplicity)
        try:
            return super().optimize(ensemble)
        finally:
            del self._resolved_charge
            del self._resolved_multiplicity

    def _optimize_conformer(
        self, molecule: Chem.Mol, conformer_id: int, atoms: Atoms
    ) -> tuple[float, bool]:
        del molecule
        if not self.orca_command:
            raise ValueError(
                "ORCA optimization requires orca_command, such as '/path/to/orca'."
            )
        try:
            from ase.calculators.orca import ORCA, OrcaProfile
        except ImportError as error:
            raise ImportError("ORCA optimization requires ASE's ORCA calculator support.") from error

        def run(directory: str) -> tuple[float, bool]:
            atoms.calc = ORCA(
                profile=OrcaProfile(command=self.orca_command),
                directory=directory,
                charge=self._resolved_charge,
                mult=self._resolved_multiplicity,
                orcasimpleinput=self.orca_simple_input,
                orcablocks=self.orca_blocks or "%pal nprocs 1 end",
            )
            optimizer = BFGS(atoms, logfile=None)
            converged = bool(optimizer.run(fmax=self.fmax, steps=self.max_steps))
            return float(atoms.get_potential_energy()) * _EV_TO_KCAL_PER_MOL, converged

        if self.workdir is None:
            with TemporaryDirectory(prefix=f"ensemblelab-orca-{conformer_id}-") as directory:
                return run(directory)

        directory = Path(self.workdir) / f"conformer-{conformer_id}"
        directory.mkdir(parents=True, exist_ok=True)
        return run(str(directory))

    def _history_details(self) -> dict[str, Any]:
        return {
        "max_steps": self.max_steps,
        "fmax": self.fmax,
        "charge": self._resolved_charge,
        "multiplicity": self._resolved_multiplicity,
        "orca_simple_input": self.orca_simple_input,
        "orca_blocks": self.orca_blocks,
        }
    # add n_cores and memory through orca? affect reproducibility and performance on other machines
