"""ASCII tables used by the human-readable display functions."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from .formatters import format_energy, format_value


def conformer_table(conformers: Iterable[Any]) -> str:
    """Render conformer-level values as a compact ASCII table.

    Relative energies are calculated only when every conformer has an energy;
    this avoids implying a complete ranking for a partially optimized ensemble.
    """
    rows = list(conformers)
    energies = [conformer.energy for conformer in rows]
    complete_energies = bool(rows) and all(energy is not None for energy in energies)
    minimum_energy = min(energies) if complete_energies else None

    headers = ("ID", "Energy", "Delta E (kcal/mol)", "Method", "Converged", "Atoms")
    values = []
    for conformer in rows:
        delta_energy = (
            "N/A"
            if minimum_energy is None or conformer.energy is None
            else f"{conformer.energy - minimum_energy:.3f}"
        )
        values.append(
            (
                str(conformer.id),
                format_energy(conformer.energy, conformer.energy_unit),
                delta_energy,
                format_value(conformer.optimization_method),
                format_value(conformer.optimization_converged),
                str(len(conformer.atoms)),
            )
        )
    return _ascii_table(headers, values)


def metadata_table(metadata: Mapping[str, Any]) -> str:
    """Render metadata as a two-column ASCII table."""
    rows = [(str(key), format_value(value)) for key, value in metadata.items()]
    return _ascii_table(("Key", "Value"), rows)


def _ascii_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    """Build a bordered table while keeping every column aligned."""
    widths = [len(header) for header in headers]
    for row in rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    separator = "+" + "+".join("-" * (width + 2) for width in widths) + "+"
    lines = [
        separator,
        "| " + " | ".join(header.ljust(width) for header, width in zip(headers, widths)) + " |",
        separator,
    ]
    lines.extend(
        "| " + " | ".join(value.ljust(width) for value, width in zip(row, widths)) + " |"
        for row in rows
    )
    lines.append(separator)
    return "\n".join(lines)
