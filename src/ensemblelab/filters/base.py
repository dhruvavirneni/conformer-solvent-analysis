from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Any, Sequence
from copy import deepcopy

from ensemblelab.generators import Conformer, Ensemble


class BaseFilter(ABC):
    """Abstract base class for ensemble filters."""
    @property
    @abstractmethod
    def method(self) -> str:
        """Returns filter method name"""
        ...

    @abstractmethod
    def _history_details(self) -> dict[str, Any]:
        """Return filter settings to store in ensemble provenance"""
        ...

    @staticmethod
    def _append_history(
        metadata: dict[str, Any],
        process: str,
        method: str,
        **details: Any,
    ) -> None:
        """Append a workflow event to the ensemble history"""
        history = list(metadata.get("history", []))
        history.append(
            {
                "process": process,
                "method": method,
                **details,
            })
        metadata["history"] = history

    @abstractmethod
    def apply(self, ensemble: Ensemble) -> Ensemble:
        """Apply the filter to an ensemble and return a new filtered ensemble."""
        ...
    @staticmethod
    def _validate_ensemble(ensemble: Ensemble) -> None:
        """Validate that the ensemble is suitable for filtering."""
        if not isinstance(ensemble, Ensemble):
            raise TypeError("Input must be an instance of Ensemble.")
        if not ensemble.conformers:
            raise ValueError("Ensemble must contain at least one conformer.")

    # reusable function for building new ensemble following filtration; filter_metadata is specific to filter type applied and is added to the ensemble metadata for provenance tracking
    def _build_filtered_ensemble(self, ensemble: Ensemble, conformers: Sequence[Conformer], filter_metadata: dict[str, Any]) -> Ensemble:
        """Construct a new Ensemble with filtered conformers and updated metadata."""
        metadata = deepcopy(ensemble.metadata)

        # update filter metadata history of ensemble
        self._append_history(
            metadata,
            process="filter",
            method=self.method,
            n_input_conformers=len(ensemble.conformers),
            n_output_conformers=len(conformers),
            **self._history_details(),
        )

        metadata["n_conformers"] = len(conformers)
        
        return Ensemble(
            smiles=ensemble.smiles,
            molecule=ensemble.molecule,
            conformers=tuple(conformers),
            metadata=metadata
        )