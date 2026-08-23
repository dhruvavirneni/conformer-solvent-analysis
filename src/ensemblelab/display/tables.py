"""Plain-text tables used by :mod:`ensemblelab.display.summaries`."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .formatters import format_bool, format_energy, format_value


def conformer_table(conformers: Sequence[Any]) -> str:
    """Return a table of conformer-level, stored optimization results."""
    rows = list(conformers)
    energies = [conformer.energy for conformer in rows]
    computed = [float(energy) for energy in energies if energy is not None]
    units = {conformer.energy_unit for conformer in rows if conformer.energy is not None}
    use_relative_energy = len(computed) == len(rows) and len(units) == 1
    minimum_energy = min(computed) if use_relative_energy else None
    unit = next(iter(units), None) if use_relative_energy else None
    energy_heading = f"Delta E ({unit})" if unit else "Energy"

    table_rows = [
        (
            str(conformer.id),
            format_energy(
                None if conformer.energy is None else float(conformer.energy) - minimum_energy,
                unit if use_relative_energy else conformer.energy_unit,
            ),
            format_value(conformer.optimization_method),
            format_bool(conformer.optimization_converged),
        )
        for conformer in rows
    ]
    return _table(("ID", energy_heading, "Method", "Converged"), table_rows)


def metadata_table(metadata: Mapping[str, Any]) -> str:
    """Return all metadata as a two-column table without interpretation."""
    return _table(
        ("Key", "Value"),
        [(str(key), format_value(value)) for key, value in metadata.items()],
    )


def _table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """Render an aligned ASCII table with stable, dependency-free formatting."""
    widths = [len(header) for header in headers]
    for row in rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    def render(row: Sequence[str]) -> str:
        return "  ".join(value.ljust(widths[index]) for index, value in enumerate(row)).rstrip()

    lines = [render(headers), render(tuple("-" * width for width in widths))]
    lines.extend(render(row) for row in rows)
    return "\n".join(lines)
