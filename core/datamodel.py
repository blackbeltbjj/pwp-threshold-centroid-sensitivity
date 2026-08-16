# -*- coding: utf-8 -*-
"""Typed, validated scientific-product models for module interoperability."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping
import numpy as np
import pandas as pd

@dataclass(frozen=True)
class ProductIdentity:
    product_type: str
    program: str
    program_version: str
    threshold_c: float | None = None
    variable: str | None = None
    pair_key: str | None = None

@dataclass(frozen=True)
class ScientificProduct:
    identity: ProductIdentity
    dates: pd.DatetimeIndex
    metadata: Mapping[str, Any] = field(default_factory=dict)
    source_files: tuple[Path, ...] = field(default_factory=tuple)

    def validate(self) -> None:
        if self.dates.empty: raise ValueError("Scientific product has no dates.")
        if not self.dates.is_monotonic_increasing: raise ValueError("Dates must be monotonic increasing.")
        if self.dates.has_duplicates: raise ValueError("Dates must not contain duplicates.")

@dataclass(frozen=True)
class WaveletProduct(ScientificProduct):
    scales_days: np.ndarray = field(default_factory=lambda: np.empty(0))
    periods_days: np.ndarray = field(default_factory=lambda: np.empty(0))
    coi_days: np.ndarray = field(default_factory=lambda: np.empty(0))
    coefficients: np.ndarray | None = None
    power: np.ndarray | None = None
    real: np.ndarray | None = None
    imaginary: np.ndarray | None = None
    phase: np.ndarray | None = None

    def validate(self) -> None:
        super().validate()
        if self.scales_days.ndim != 1 or self.periods_days.ndim != 1: raise ValueError("Scales and periods must be one-dimensional.")
        if self.scales_days.size != self.periods_days.size: raise ValueError("Scale and period counts differ.")
        if self.coi_days.shape != (self.dates.size,): raise ValueError("COI length differs from date count.")
        expected=(self.periods_days.size, self.dates.size)
        for name, array in (("coefficients",self.coefficients),("power",self.power),("real",self.real),("imaginary",self.imaginary),("phase",self.phase)):
            if array is not None and array.shape != expected: raise ValueError(f"{name} shape {array.shape} differs from {expected}.")

@dataclass(frozen=True)
class CrossWaveletProduct(WaveletProduct):
    coherence: np.ndarray | None = None
    x_name: str = ""
    y_name: str = ""

    def validate(self) -> None:
        super().validate()
        expected=(self.periods_days.size, self.dates.size)
        if self.coherence is not None and self.coherence.shape != expected: raise ValueError(f"coherence shape {self.coherence.shape} differs from {expected}.")
