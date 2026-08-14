# ensemblelab

`ensemblelab` is an open source Python library for quantitative analysis of molecular conformational ensembles. It provides functionality for ensemble analysis (basin visualization, RMSD heatmapping, etc.) with simple object-based workflows.

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




## Modules

- **Generate:** reproducible, unoptimized conformer ensembles with aligned RDKit and ASE representations.
- **Optimize:** single-stage and hierarchical conformer refinement.
- **Descriptors:** in progress.
- **Filter:** in progress.
- **Clustering:** in progress.
- **Utils:** in progress.