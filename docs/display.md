# Display

## Overview
`Ensemble.show()` and `Conformer.show()` return read-only, human-readable summaries of stored ensemble data. They never optimize conformers, compute descriptors, or otherwise perform analysis. The default ensemble view contains the molecular summary and a compact conformer table; optional sections expose workflow provenance and raw metadata for inspection and debugging.

## Functions

- **Ensemble.show():** returns a default ensemble summary and conformer table.
- **Conformer.show():** returns stored energy, optimization, convergence, and atom-count data for one conformer.

### `Ensemble.show`

**Overview**
The default display reports SMILES, explicit atom count, conformer count, whether energy values are available, the stored optimization method, and one row per conformer. When every conformer has an energy in one unit, the table shows relative energy (Delta E) from the lowest-energy conformer; this keeps the comparison compact while the absolute energies remain on each `Conformer`.

Set `history=True` to add a workflow history section. The display supports the current generation (`processing_history`), optimizer/filter (`history`), and functional optimizer (`optimization_history`) provenance keys without modifying metadata. Set `metadata=True` to add a raw key/value metadata table. Sections can be combined, and `conformers=False` suppresses the default conformer table.

**Usage**

```python
from ensemblelab import Ensemble

ensemble = Ensemble.from_smiles("CCO", n_confs=25)
print(ensemble.show())

print(ensemble.show(history=True, metadata=True, conformers=False))
```

### `Conformer.show`

**Overview**
Returns the conformer ID, stored energy and unit, optimization method, convergence state, and number of atoms. It does not include a coordinate table; molecular geometry is retained in `conformer.atoms` for programmatic use and future visualization tools.

**Usage**

```python
conformer = ensemble.conformers[0]
print(conformer.show())
```
