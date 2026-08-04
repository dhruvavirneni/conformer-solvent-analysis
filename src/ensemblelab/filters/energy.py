from ensemblelab.generators import Ensemble
from abc import ABC, abstractmethod
from ensemblelab.filters.base import BaseFilter
from dataclasses import dataclass
from copy import deepcopy

@dataclass(frozen=True, slots=True)
class EnergyFilter(BaseFilter):

    def __post_init__(self, window: float | None = None, top_n: int | None = None):
        if (window is None) == (top_n is None):
            raise ValueError("Specify exactly one of window or top_n.")
        self.window = window
        self.top_n = top_n
  
    def apply(self, ensemble: Ensemble) -> Ensemble:
        if self.window is not None:
            metadata = deepcopy(ensemble.metadata)
            ensemble_energies = [conformer.energy for conformer in ensemble.conformers if conformer.energy is not None]
            min_energy = min(ensemble_energies)
            filtered_conformers = [conformer for conformer in ensemble.conformers if conformer.energy is not None and conformer.energy - min_energy <= self.window]
            return Ensemble(
                smiles=ensemble.smiles,
                molecule=ensemble.molecule,
                conformers=tuple(filtered_conformers),
                metadata=metadata
)
        elif self.top_n is not None:
            metadata = deepcopy(ensemble.metadata)
            ensemble_energies = [conformer.energy for conformer in ensemble.conformers if conformer.energy is not None]
            min_energy = min(ensemble_energies)
            sorted_conformers = sorted(ensemble.conformers, key=lambda c: c.energy - min_energy)
            filtered_conformers = sorted_conformers[:self.top_n]
            return Ensemble(
                smiles=ensemble.smiles,
                molecule=ensemble.molecule,
                conformers=tuple(filtered_conformers),
                metadata=metadata
)