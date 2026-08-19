# Filter

## Overview
Filters select or remove conformers from an `Ensemble` and return a new `Ensemble`; the input is not modified. Each filter preserves the molecule and SMILES, recording the input/output conformer counts and filter settings in `metadata["history"]`. Filters that use energies require energies to be present on input Conformers (via optimization).

## Functions
- **EnergyFilter:** keeps conformers within an energy window of the lowest-energy conformer, or keeps the lowest-energy `top_n` conformers.
- **PopulationFilter:** keeps conformers by individual or cumulative Boltzmann population.
- **RMSDFilter:** removes geometrically redundant conformers using an RMSD cutoff, prioritizing lower-energy conformers.


### `EnergyFilter`
**Overview**
Filters conformers by computed energy. Specify exactly one of `window` or `top_n`. With `window`, all conformers whose energy is within the specified kcal/mol window of the minimum energy are retained. With `top_n`, conformers are sorted by energy and the lowest-energy entries are retained. Conformers without an energy are excluded.

Metadata recorded: `window` or `top_n`


**Usage**
```python
from ensemblelab.filters.energy import EnergyFilter

energy_filter = EnergyFilter(window=5.0)
filtered = energy_filter.apply(ensemble)
```


### `PopulationFilter`
**Overview**
Filters conformers using normalized Boltzmann populations. Specify exactly one of `population_cutoff` or `cumulative_cutoff`; both values are fractions between 0 and 1. `population_cutoff` retains every conformer whose individual population meets the cutoff. `cumulative_cutoff` ranks conformers by population and retains entries until the cumulative population reaches the cutoff. Every conformer must have an energy.

The default temperature is 298.15 K. Populations are calculated using the kcal/mol gas constant and the supplied temperature.

Metadata recorded: `temperature`, `population_cutoff`, `cumulative_cutoff`


**Usage**
```python
from ensemblelab.filters.population import PopulationFilter

population_filter = PopulationFilter(
	cumulative_cutoff=0.95,
	temperature=298.15,
)
filtered = population_filter.apply(ensemble)
```


### `RMSDFilter`
**Overview**
Removes conformers that are within the RMSD cutoff of an already approved conformer. Conformers are sorted from lowest to highest energy, so the lowest energy representative is considered first. The default cutoff is 0.5 Angstroms. Every conformer must have a computed energy.

Metadata recorded: `rmsd_cutoff_angstrom`


**Usage**
```python
from ensemblelab.filters.rmsd import RMSDFilter

rmsd_filter = RMSDFilter(cutoff=0.5)
filtered = rmsd_filter.apply(ensemble)
```
