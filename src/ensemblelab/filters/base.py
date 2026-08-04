from ensemblelab.generators import Ensemble
from abc import ABC, abstractmethod

class BaseFilter(ABC):

    @abstractmethod
    def apply(self, ensemble: Ensemble) -> Ensemble:
        """Apply the filter to an ensemble and return a new filtered ensemble."""
        pass
