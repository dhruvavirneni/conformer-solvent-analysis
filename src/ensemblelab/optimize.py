"""Conformer geometry optimization with RDKit and optional ASE backends."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
import warnings

import numpy as np
from ase.optimize import BFGS
from rdkit import Chem
from rdkit.Chem import AllChem

from .generators import Conformer, Ensemble

_EV_TO_KCAL_PER_MOL = 23.0605478306
_Method = Literal["MMFF", "UFF", "GFN2-xTB", "ORCA"]


def optimize(
    ensemble: Ensemble,
    method: str = "GFN2-xTB",
    fmax: float = 0.05,
    solvent: str | None = None,
    *,
    max_steps: int = 500,
    charge: int | None = None,
    multiplicity: int | None = None,
    orca_command: str | None = None,
    orca_simple_input: str = "r2SCAN-3c",
    orca_blocks: str | None = None,
    workdir: str | Path | None = None,
) -> Ensemble:
    """Optimize every conformer and return a new ensemble.

    The input ensemble is never modified. The returned conformers hold the
    optimized ASE coordinates and a per-conformer energy in kcal/mol. Their
    corresponding coordinates in the returned RDKit molecule are updated too.

    Parameters
    ----------
    ensemble
        Ensemble to optimize. It must contain at least one conformer.
    method
        Optimization backend: ``"MMFF"``, ``"UFF"``, ``"GFN2-xTB"``, or
        ``"ORCA"``. ``"xtb"`` is accepted as an alias for ``"GFN2-xTB"``.
    fmax
        Maximum force in eV/angstrom for ASE-backed GFN2-xTB and ORCA
        optimizations. It is not used by RDKit MMFF or UFF.
    solvent
        ALPB solvent name for GFN2-xTB. It must be ``None`` for MMFF, UFF,
        and ORCA because their solvent models require backend-specific input.
    max_steps
        Maximum geometry-optimization steps or RDKit force-field iterations.
    charge, multiplicity
        Electronic state for GFN2-xTB and ORCA. Charge defaults to the RDKit
        formal charge. Multiplicity defaults to one plus radical electrons.
    orca_command
        Path or command for the ORCA executable. Required for ``"ORCA"``.
    orca_simple_input, orca_blocks
        ORCA input settings. Use these to select a method, basis, parallelism,
        and any backend-specific solvent treatment.
    workdir
        Parent directory for ORCA calculation folders. If omitted, temporary
        folders are used and removed after each conformer finishes.

    Returns
    -------
    Ensemble
        A new ensemble with optimized geometries, per-conformer energies, and
        optimization provenance in ``metadata``.

    Raises
    ------
    ImportError
        If an optional backend dependency is unavailable.
    RuntimeError
        If a backend cannot construct or run a calculation.
    ValueError
        If parameters are incompatible with the requested backend.
    """
    canonical_method = _normalise_method(method)
    _validate_inputs(ensemble, canonical_method, fmax, solvent, max_steps)

    optimized_molecule = Chem.Mol(ensemble.molecule)
    optimized_conformers: list[Conformer] = []
    convergence: dict[int, bool] = {}
    resolved_charge = _resolve_charge(ensemble, charge)
    resolved_multiplicity = _resolve_multiplicity(ensemble, multiplicity)

    for conformer in ensemble.conformers:
        atoms = conformer.atoms.copy()
        _set_rdkit_positions(optimized_molecule, conformer.id, atoms.positions)

        if canonical_method == "MMFF":
            energy, converged = _optimize_rdkit(
                optimized_molecule, conformer.id, method="MMFF", max_steps=max_steps
            )
        elif canonical_method == "UFF":
            energy, converged = _optimize_rdkit(
                optimized_molecule, conformer.id, method="UFF", max_steps=max_steps
            )
        elif canonical_method == "GFN2-xTB":
            energy, converged = _optimize_xtb(
                atoms,
                fmax=fmax,
                max_steps=max_steps,
                charge=resolved_charge,
                multiplicity=resolved_multiplicity,
                solvent=solvent,
            )
            _set_rdkit_positions(optimized_molecule, conformer.id, atoms.positions)
        else:
            energy, converged = _optimize_orca(
                atoms,
                fmax=fmax,
                max_steps=max_steps,
                charge=resolved_charge,
                multiplicity=resolved_multiplicity,
                command=orca_command,
                simple_input=orca_simple_input,
                blocks=orca_blocks,
                parent_workdir=workdir,
                conformer_id=conformer.id,
            )
            _set_rdkit_positions(optimized_molecule, conformer.id, atoms.positions)

        if canonical_method in {"MMFF", "UFF"}:
            atoms.positions[:] = np.asarray(
                optimized_molecule.GetConformer(conformer.id).GetPositions(), dtype=float
            )

        optimized_conformers.append(
            Conformer(
                id=conformer.id,
                atoms=atoms,
                energy=energy,
                energy_unit="kcal/mol",
                optimization_method=canonical_method,
                optimization_converged=converged,
            )
        )
        convergence[conformer.id] = converged

    unconverged_ids = [identifier for identifier, done in convergence.items() if not done]
    if unconverged_ids:
        warnings.warn(
            "Optimization reached max_steps before convergence for conformer IDs: "
            f"{unconverged_ids}.",
            RuntimeWarning,
            stacklevel=2,
        )

    metadata = deepcopy(ensemble.metadata)
    history = list(metadata.get("optimization_history", []))
    history.append(
        {
            "method": canonical_method,
            "fmax_eV_per_angstrom": fmax if canonical_method in {"GFN2-xTB", "ORCA"} else None,
            "max_steps": max_steps,
            "solvent": solvent,
            "charge": resolved_charge if canonical_method in {"GFN2-xTB", "ORCA"} else None,
            "multiplicity": resolved_multiplicity if canonical_method in {"GFN2-xTB", "ORCA"} else None,
            "orca_simple_input": orca_simple_input if canonical_method == "ORCA" else None,
            "orca_blocks": orca_blocks if canonical_method == "ORCA" else None,
            "energy_unit": "kcal/mol",
            "converged_conformer_ids": [identifier for identifier, done in convergence.items() if done],
            "unconverged_conformer_ids": unconverged_ids,
        }
    )
    metadata.update(
        {
            "optimization_status": "optimized" if not unconverged_ids else "not_fully_converged",
            "energy_status": "computed",
            "energy_unit": "kcal/mol",
            "optimization_history": history,
        }
    )
    return Ensemble(
        smiles=ensemble.smiles,
        molecule=optimized_molecule,
        conformers=tuple(optimized_conformers),
        metadata=metadata,
    )


def _normalise_method(method: str) -> _Method:
    if not isinstance(method, str):
        raise ValueError("method must be a string.")
    aliases: dict[str, _Method] = {
        "mmff": "MMFF",
        "uff": "UFF",
        "xtb": "GFN2-xTB",
        "gfn2-xtb": "GFN2-xTB",
        "gfn2_xtb": "GFN2-xTB",
        "orca": "ORCA",
    }
    try:
        return aliases[method.strip().lower()]
    except KeyError as error:
        choices = ", ".join(sorted(set(aliases.values())))
        raise ValueError(f"Unsupported optimization method {method!r}. Choose: {choices}.") from error


def _validate_inputs(
    ensemble: Ensemble,
    method: _Method,
    fmax: float,
    solvent: str | None,
    max_steps: int,
) -> None:
    if not isinstance(ensemble, Ensemble):
        raise TypeError("ensemble must be an Ensemble instance.")
    if not ensemble.conformers:
        raise ValueError("ensemble must contain at least one conformer.")
    if not isinstance(fmax, (float, int)) or isinstance(fmax, bool) or fmax <= 0:
        raise ValueError("fmax must be a positive number in eV/angstrom.")
    if isinstance(max_steps, bool) or not isinstance(max_steps, int) or max_steps < 1:
        raise ValueError("max_steps must be a positive integer.")
    if solvent is not None and (not isinstance(solvent, str) or not solvent.strip()):
        raise ValueError("solvent must be a non-empty string or None.")
    if solvent is not None and method != "GFN2-xTB":
        raise ValueError(
            "solvent is currently supported only for GFN2-xTB through ALPB. "
            "Configure ORCA solvent input with orca_simple_input or orca_blocks."
        )


def _resolve_charge(ensemble: Ensemble, charge: int | None) -> int:
    if charge is None:
        return int(Chem.GetFormalCharge(ensemble.molecule))
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError("charge must be an integer or None.")
    return charge


def _resolve_multiplicity(ensemble: Ensemble, multiplicity: int | None) -> int:
    if multiplicity is None:
        return 1 + sum(atom.GetNumRadicalElectrons() for atom in ensemble.molecule.GetAtoms())
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 1:
        raise ValueError("multiplicity must be a positive integer or None.")
    return multiplicity


def _optimize_rdkit(
    molecule: Chem.Mol, conformer_id: int, *, method: Literal["MMFF", "UFF"], max_steps: int
) -> tuple[float, bool]:
    """Optimize one RDKit conformer and return kcal/mol energy and convergence."""
    if method == "MMFF":
        properties = AllChem.MMFFGetMoleculeProperties(molecule)
        if properties is None:
            raise RuntimeError("MMFF parameters are unavailable for this molecule.")
        force_field = AllChem.MMFFGetMoleculeForceField(
            molecule, properties, confId=conformer_id
        )
    else:
        force_field = AllChem.UFFGetMoleculeForceField(molecule, confId=conformer_id)

    if force_field is None:
        raise RuntimeError(f"{method} could not construct a force field for conformer {conformer_id}.")
    status = force_field.Minimize(maxIts=max_steps)
    if status not in {0, 1}:
        raise RuntimeError(f"{method} optimization failed for conformer {conformer_id} (status {status}).")
    return float(force_field.CalcEnergy()), status == 0


def _optimize_xtb(
    atoms,
    *,
    fmax: float,
    max_steps: int,
    charge: int,
    multiplicity: int,
    solvent: str | None,
) -> tuple[float, bool]:
    """Optimize one conformer with TBLite and return energy in kcal/mol."""
    try:
        from tblite.ase import TBLite
    except ImportError as error:
        raise ImportError(
            "GFN2-xTB optimization requires the optional 'tblite' package. "
            "Install ensemblelab with its xTB extra."
        ) from error

    calculator_options = {"method": "GFN2-xTB", "charge": charge, "multiplicity": multiplicity}
    if solvent is not None:
        calculator_options["solvation"] = ("alpb", solvent)
    atoms.calc = TBLite(**calculator_options)
    optimizer = BFGS(atoms, logfile=None)
    converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
    return float(atoms.get_potential_energy()) * _EV_TO_KCAL_PER_MOL, converged


def _optimize_orca(
    atoms,
    *,
    fmax: float,
    max_steps: int,
    charge: int,
    multiplicity: int,
    command: str | None,
    simple_input: str,
    blocks: str | None,
    parent_workdir: str | Path | None,
    conformer_id: int,
) -> tuple[float, bool]:
    """Optimize one conformer with ASE's optional ORCA calculator."""
    if not command:
        raise ValueError("ORCA optimization requires orca_command, such as '/path/to/orca'.")
    try:
        from ase.calculators.orca import ORCA, OrcaProfile
    except ImportError as error:
        raise ImportError("ORCA optimization requires ASE's ORCA calculator support.") from error

    def run(directory: str) -> tuple[float, bool]:
        atoms.calc = ORCA(
            profile=OrcaProfile(command=command),
            directory=directory,
            charge=charge,
            mult=multiplicity,
            orcasimpleinput=simple_input,
            orcablocks=blocks or "%pal nprocs 1 end",
        )
        optimizer = BFGS(atoms, logfile=None)
        converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
        return float(atoms.get_potential_energy()) * _EV_TO_KCAL_PER_MOL, converged

    if parent_workdir is None:
        with TemporaryDirectory(prefix=f"ensemblelab-orca-{conformer_id}-") as directory:
            return run(directory)

    directory = Path(parent_workdir) / f"conformer-{conformer_id}"
    directory.mkdir(parents=True, exist_ok=True)
    return run(str(directory))


def _set_rdkit_positions(
    molecule: Chem.Mol, conformer_id: int, positions: np.ndarray
) -> None:
    """Copy optimized ASE coordinates into the corresponding RDKit conformer."""
    conformer = molecule.GetConformer(conformer_id)
    for atom_index, (x, y, z) in enumerate(np.asarray(positions, dtype=float)):
        conformer.SetAtomPosition(atom_index, (float(x), float(y), float(z)))
