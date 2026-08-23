"""Small, shared formatting helpers for human-readable displays."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


def format_bool(value: bool | None) -> str:
    """Format an optional boolean for a compact text display."""
    if value is True:
        return "✓"
    if value is False:
        return "✗"
    return "N/A"


def format_energy(value: float | None, unit: str | None = None) -> str:
    """Format an optional energy value to three decimal places."""
    if value is None:
        return "N/A"
    result = f"{value:.3f}"
    return f"{result} {unit}" if unit else result


def format_value(value: Any) -> str:
    """Format metadata without requiring metadata to be JSON serializable."""
    if value is None:
        return "N/A"
    if isinstance(value, bool):
        return format_bool(value)
    if isinstance(value, Mapping):
        return ", ".join(f"{key}={format_value(item)}" for key, item in value.items())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return ", ".join(format_value(item) for item in value)
    return str(value)
