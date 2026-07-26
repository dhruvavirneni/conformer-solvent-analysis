# ensemblelab

`ensemblelab` is an open-source Python library for quantitative analysis of molecular conformational ensembles. It fills the gap between creating conformers and describing flexibility at the ensemble level.

## Current capability

Generate a reproducible, unoptimized conformer ensemble with aligned RDKit and ASE representations:

```python
from ensemblelab import generate

ensemble = generate("CCO", n_confs=20)
print(ensemble.metadata["n_conformers_generated"])
assert ensemble.conformers[0].energy is None  # assigned by optimize() later
```

Optimize an ensemble without modifying the original object:

```python
from ensemblelab import optimize

optimized = optimize(ensemble, method="MMFF")
print(optimized.conformers[0].energy)  # kcal/mol
```

## Optimization backends

- `MMFF` and `UFF` use RDKit and are available in the base install.
- `GFN2-xTB` uses TBLite with ASE BFGS optimization. Install it with `pip install -e ".[xtb]"`. The optional `solvent` argument selects TBLite ALPB solvation.
- `ORCA` uses ASE's ORCA calculator. Pass `orca_command` and optional `orca_simple_input`, `orca_blocks`, `charge`, `multiplicity`, and `workdir` settings to `optimize()`.

All backends return a new `Ensemble`. Optimized energies are stored as kcal/mol on each conformer, and the optimization method, convergence status, and settings are recorded in metadata.

The returned `Ensemble` stores the input SMILES, a hydrogenated RDKit molecule containing all conformers, ASE `Atoms` snapshots, and generation provenance. Each conformer owns its own placeholder `energy` value (`ensemble.conformers[0].energy`), which remains `None` until optimization.

## Install for development

```bash
pip install -e ".[dev]"
pytest
```

## Roadmap

1. Generate, optimize, filter, and cluster conformer ensembles.
2. Define and validate an interpretable conformational fingerprint.
3. Benchmark across molecule classes and solvent models.
4. Investigate experimental-property correlations and molecular-ML embeddings.
