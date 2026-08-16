# -*- coding: utf-8 -*-
"""Reusable validation contracts that never mutate scientific data."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def require(self, condition: bool, message: str) -> None:
        if not condition: self.errors.append(message)

    def warn(self, condition: bool, message: str) -> None:
        if not condition: self.warnings.append(message)

    def raise_if_invalid(self, heading: str = "Validation failed") -> None:
        if self.errors:
            raise ValueError(heading + ":\n" + "\n".join(f"  - {x}" for x in self.errors))

    def lines(self) -> tuple[str, ...]:
        return tuple(self.messages + [f"WARNING: {x}" for x in self.warnings])

def require_files(paths: Iterable[Path], label: str = "Required files") -> tuple[Path, ...]:
    normalized=tuple(Path(p) for p in paths)
    missing=tuple(p for p in normalized if not p.is_file())
    if missing:
        raise FileNotFoundError(label + " are missing:\n" + "\n".join(f"  - {p}" for p in missing))
    return normalized

def ensure_directories(paths: Iterable[Path]) -> tuple[Path, ...]:
    normalized=tuple(Path(p) for p in paths)
    for path in normalized: path.mkdir(parents=True, exist_ok=True)
    return normalized

def require_unique(values: Sequence[object], label: str) -> None:
    if len(values) != len(set(values)):
        raise ValueError(f"{label} must be unique.")

def require_between(value: float, minimum: float, maximum: float, label: str, *, inclusive: bool=True) -> None:
    valid = minimum <= value <= maximum if inclusive else minimum < value < maximum
    if not valid:
        brackets = "[ ]" if inclusive else "( )"
        raise ValueError(f"{label}={value} must lie in {brackets[0]}{minimum}, {maximum}{brackets[1]}.")
