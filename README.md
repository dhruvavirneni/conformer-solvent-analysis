# ensemblelab

`ensemblelab` is an open-source Python library for quantitative analysis of molecular conformational ensembles. It fills the gap between creating conformers and describing flexibility at the ensemble level.

## Install for development

```bash
pip install -e ".[dev]"
pytest
```

## Datatypes

- `Conformer`: EnsembleLab native object type storing ASE `Atoms`-based molecules; records conformer energy in kcal/mol and optimization type.
- `Ensemble`: EnsembleLab native object storing `Conformer` objects.
- [`Mol`](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol): RDKit-native molecule object.
- [`Atoms`](https://ase-lib.org/ase/atoms.html): ASE-native atomic structure object.

## Current API

The package currently provides conformer generation, single-stage optimization, and hierarchical optimization. Import the public functions directly from `ensemblelab`:

```python
from ensemblelab import generate, optimize, hierarchical_optimize
```

### Generate

`generate(smiles, n_confs=25)` builds explicit-hydrogen conformers with RDKit ETKDGv3. Generation uses a fixed random seed of `42`, so the chosen settings are reproducible and stored in `ensemble.metadata`.

```python
ensemble = generate("CCO", n_confs=20)

print(ensemble.conformer_ids)
print(ensemble.metadata["energy_status"])  # "uncomputed"
```

`Ensemble.from_smiles()` is the equivalent object-oriented constructor:

```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=20)
```

An `Ensemble` keeps conformer-level data aligned:

- `ensemble.molecule` is the hydrogenated RDKit molecule with all coordinates.
- `ensemble.conformers` is a tuple of `Conformer` records with ASE `Atoms`.
- `ensemble.conformer_ids` gives the stable RDKit conformer IDs.
- `ensemble.rdkit_conformer(conformer_id)` retrieves a conformer by its ID.

Generated conformers have `energy=None`. Optimization assigns energies in `kcal/mol`, along with the method and convergence result.

### Optimize

`optimize()` accepts a method name and returns a new `Ensemble`; it does not change the input ensemble. The returned RDKit and ASE coordinates remain synchronized.

```python
optimized = optimize(ensemble, method="MMFF", max_steps=500)

for conformer in optimized.conformers:
    print(conformer.id, conformer.energy, conformer.energy_unit)

print(optimized.metadata["optimization_history"][-1])
```

Available method names are `"MMFF"`, `"UFF"`, `"GFN2-xTB"`, and `"ORCA"`. `"xtb"`, `"gfn2-xtb"`, and `"gfn2_xtb"` are accepted aliases for `"GFN2-xTB"`.

- **MMFF / UFF:** RDKit force-field optimization. `max_steps` is the maximum force-field iteration count.
- **GFN2-xTB:** uses TBLite through ASE. Install it with `pip install -e ".[xtb]"`. `fmax` is in eV/angstrom; `solvent` selects an ALPB solvent model. `charge` defaults to RDKit formal charge, and `multiplicity` defaults from radical electrons.
- **ORCA:** uses ASE's ORCA calculator. Supply `orca_command`; optional `orca_simple_input`, `orca_blocks`, and `workdir` configure the calculation.

All backends return a new `Ensemble`. Optimized energies are stored as kcal/mol on each conformer, and the optimization method, convergence status, and settings are recorded in metadata. If a calculation reaches its maximum number of steps, the returned conformer is still retained with `optimization_converged=False`, and the conformer ID is listed in optimization metadata. Optional backends raise actionable errors when their dependency or executable is unavailable.

### Hierarchical optimization

Use `hierarchical_optimize()` to screen an ensemble with one to three optimizer instances. Non-final stages optimize all current conformers, rank them by energy, and retain the lowest-energy subset before the next stage. The final stage is also trimmed to `target_size` when needed.

```python
from ensemblelab.optimizers import GFN2xTBOptimizer, MMFFOptimizer

refined = hierarchical_optimize(
    ensemble,
    workflow=[
        MMFFOptimizer(max_steps=500),
        GFN2xTBOptimizer(fmax=0.05, max_steps=500),
    ],
    target_size=25,
    retention_rates=[0.2],
)
```

For every non-final stage, the retained count is:

```text
max(target_size, ceil(current_size * retention_rate))
```

`retention_rates` therefore needs one value for each non-final stage. For a single optimizer stage, pass an empty sequence:

```python
from ensemblelab.optimizers import MMFFOptimizer

screened = hierarchical_optimize(
    ensemble,
    workflow=[MMFFOptimizer()],
    target_size=25,
    retention_rates=[],
)
```

The returned ensemble records the workflow in `metadata["hierarchical_optimization"]`, including the workflow names, target size, retention rates, and per-stage input/output sizes, energy range, and retained energy cutoff. Each optimizer also appends its normal record to `metadata["optimization_history"]`.

## Modules

- **Generate:** reproducible, unoptimized conformer ensembles with aligned RDKit and ASE representations.
- **Optimize:** single-stage and hierarchical conformer refinement.
- **Descriptors:** in progress.
- **Filter:** in progress.
- **Clustering:** in progress.
- **Utils:** in progress.

## Roadmap

1. Generate, optimize, filter, and cluster conformer ensembles.
2. Define and validate an interpretable conformational fingerprint.
3. Benchmark across molecule classes and solvent models.
4. Investigate experimental-property correlations and molecular-ML embeddings.
