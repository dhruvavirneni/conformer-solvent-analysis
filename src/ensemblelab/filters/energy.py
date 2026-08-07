from ensemblelab.generators import Ensemble
from ensemblelab.filters.base import BaseFilter
from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class EnergyFilter(BaseFilter):

    def __post_init__(self, window: float | None = None, top_n: int | None = None):
        if (window is None) == (top_n is None):
            raise ValueError("Specify exactly one of window or top_n.")
        self.window = window
        self.top_n = top_n
  
    def apply(self, ensemble: Ensemble) -> Ensemble:
        self._validate_ensemble(ensemble)
        if self.window is not None:
            ensemble_energies = [conformer.energy for conformer in ensemble.conformers if conformer.energy is not None]
            min_energy = min(ensemble_energies)
            filtered_conformers = [conformer for conformer in ensemble.conformers if conformer.energy is not None and conformer.energy - min_energy <= self.window]
            return self._build_filtered_ensemble(ensemble, filtered_conformers, {"window": self.window}
)
        elif self.top_n is not None:
            ensemble_energies = [conformer.energy for conformer in ensemble.conformers if conformer.energy is not None]
            min_energy = min(ensemble_energies)
            sorted_conformers = sorted(ensemble.conformers, key=lambda c: c.energy)
            filtered_conformers = sorted_conformers[:self.top_n]
            return self._build_filtered_ensemble(ensemble, filtered_conformers, {"top_n": self.top_n})