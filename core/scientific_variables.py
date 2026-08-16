# -*- coding: utf-8 -*-
"""Authoritative scientific-variable registry shared by PWP modules."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Final, Iterator

@dataclass(frozen=True)
class ScientificVariable:
    key: str
    column: str
    label: str
    units: str
    symbol: str
    circular: bool = False

VARIABLE_REGISTRY: Final[dict[str, ScientificVariable]] = {
    "lon_360": ScientificVariable("lon_360", "lon_360", "PWP centroid longitude", "degrees east", "λc", True),
    "lat": ScientificVariable("lat", "lat", "PWP centroid latitude", "degrees", "φc"),
    "area_km2": ScientificVariable("area_km2", "area_km2", "PWP total area", "km²", "A"),
    "mean_pwp_sst_c": ScientificVariable("mean_pwp_sst_c", "mean_pwp_sst_c", "Area-weighted mean PWP SST", "°C", "T̄PWP"),
}
PWP_VARIABLE_KEYS: Final[tuple[str, ...]] = tuple(VARIABLE_REGISTRY)
VARIABLE_LABELS: Final[dict[str, str]] = {k:v.label for k,v in VARIABLE_REGISTRY.items()}
VARIABLE_UNITS: Final[dict[str, str]] = {k:v.units for k,v in VARIABLE_REGISTRY.items()}

def get_variable(key: str) -> ScientificVariable:
    try: return VARIABLE_REGISTRY[key]
    except KeyError as error:
        raise KeyError(f"Unknown scientific variable {key!r}; available: {', '.join(VARIABLE_REGISTRY)}") from error

def iter_variables() -> Iterator[ScientificVariable]:
    return iter(VARIABLE_REGISTRY.values())
