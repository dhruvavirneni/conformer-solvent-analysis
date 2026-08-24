# Generate

## Overview
`generate()` creates a new unoptimized `Ensemble` from a SMILES string using RDKit ETKDGv3 with a fixed random seed of 42. Each generated conformer is stored as a `Conformer` with an ASE `Atoms` geometry and optional energy metadata that is populated later by optimization steps. The returned ensemble retains the original SMILES, the RDKit molecule object, conformer IDs, and provenance information in `metadata["processing_history"]`.

## Functions
- **Conformer:** immutable per-conformer object with geometry, energy, and optimization metadata.
- **Ensemble:** aligned conformer collection with the original molecule and provenance.
- **generate():** generates a new ensemble from a SMILES string.
- **Ensemble.from_smiles():** convenience constructor equivalent to `generate()`.


### `Conformer`
**Overview:** 
A single generated conformer and its conformer-level data. The `energy` field is `None` until an optimization backend assigns a value. The geometry is stored as an ASE `Atoms` object, while `energy_unit`, `optimization_method`, and `optimization_converged` are attached to the same conformer rather than to a separate ensemble-level array.

Metadata recorded: `id`, `atoms`, `energy`, `energy_unit`, `optimization_method`, `optimization_converged`


**Usage**
```python
from ase import Atoms
from ensemblelab import Conformer

conformer = Conformer(
    id=0,
    atoms=Atoms("H2O", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.9], [0.0, 0.9, 0.0]]),
    energy=None,
    energy_unit=None,
    optimization_method=None,
    optimization_converged=None,
)
```


### `Ensemble`
**Overview:** 
An aligned conformational ensemble containing the input SMILES, the RDKit molecule, and a tuple of `Conformer` objects. The ensemble enforces unique conformer IDs and preserves one explicit conformer index for all downstream analyses. The molecular geometry and energies remain synchronized with the conformer-level metadata.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `processing_history`


**Usage**
```python
from ensemblelab import Ensemble

ensemble = Ensemble(
    smiles="CCO",
    molecule=molecule,
    conformers=(conformer_1, conformer_2),
    metadata={
        "n_conformers": 2,
        "processing_history": [],
    },
)
```


### `generate`
**Overview:** 
Creates an unoptimized ensemble by parsing a SMILES string, adding explicit hydrogens, and embedding multiple conformers with RDKit ETKDGv3. The default random seed is fixed at 42 for reproducibility, and the ensemble stores the canonicalized molecule and generation provenance in `metadata["processing_history"]`. Validation rejects empty or invalid SMILES strings and requires a positive integer value for `n_confs`.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `processing_history`


**Usage**
```python
from ensemblelab import generate

ensemble = generate("CCO", n_confs=25)
```


### `Ensemble.from_smiles`
**Overview:** 
Class method equivalent to calling `generate()`. It is the object-oriented entry point for creating an unoptimized `Ensemble` from a SMILES string and is intended to fit naturally into the broader `Ensemble` workflow for optimization, filtering, and analysis.

Metadata recorded: `n_conformers`, `optimization_status`, `energy_status`, `energy_unit`, `rdkit_version`, `processing_history`


**Usage**
```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=25)
```
