# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

FILE
    config/wavelet_config.py

VERSION
    4.0.0

PURPOSE
    Authoritative scientific and graphical configuration for Program 11.

SCIENTIFIC BASIS
    Torrence and Compo (1998), Domingues et al. (2005), Gu and Philander
    (1995), and Torrence and Webster (1999).

USER POLICY
    Edit spectral-band limits and wavelet-figure families in this file rather
    than editing Program 11.

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal

PeriodUnit = Literal["days", "years"]
BandNormalization = Literal["absolute", "mean", "standard_score", "maximum"]

DAYS_PER_YEAR: Final[float] = 365.2425


@dataclass(frozen=True)
class WaveletBand:
    """Definition of one user-selected Fourier-period band."""

    key: str
    label: str
    short_label: str
    minimum_period: float
    maximum_period: float
    period_unit: PeriodUnit = "days"
    enabled: bool = True
    order: int = 1
    show_on_scalogram: bool = True
    show_on_global_spectrum: bool = True

    @property
    def minimum_period_days(self) -> float:
        if self.period_unit == "years":
            return self.minimum_period * DAYS_PER_YEAR
        return self.minimum_period

    @property
    def maximum_period_days(self) -> float:
        if self.period_unit == "years":
            return self.maximum_period * DAYS_PER_YEAR
        return self.maximum_period

    @property
    def period_range_text(self) -> str:
        return (
            f"{self.minimum_period:g}–{self.maximum_period:g} "
            f"{self.period_unit}"
        )


# Between one and six bands may be enabled.
WAVELET_BANDS: Final[tuple[WaveletBand, ...]] = (
    WaveletBand(
        key="annual",
        label="Annual cycle",
        short_label="Annual",
        minimum_period=330.0,
        maximum_period=400.0,
        period_unit="days",
        enabled=True,
        order=1,
    ),
    WaveletBand(
        key="quasi_biennial",
        label="Quasi-biennial variability",
        short_label="QBO-like",
        minimum_period=1.5,
        maximum_period=3.0,
        period_unit="years",
        enabled=True,
        order=2,
    ),
    WaveletBand(
        key="enso",
        label="ENSO-scale variability",
        short_label="ENSO",
        minimum_period=3.0,
        maximum_period=8.0,
        period_unit="years",
        enabled=True,
        order=3,
    ),
    WaveletBand(
        key="interdecadal",
        label="Interdecadal variability",
        short_label="Interdecadal",
        minimum_period=8.0,
        maximum_period=16.0,
        period_unit="years",
        enabled=True,
        order=4,
    ),
)

MOTHER_WAVELET_NAME: Final[str] = "Complex Morlet"
MORLET_OMEGA0: Final[float] = 6.0
SAMPLING_INTERVAL_DAYS: Final[float] = 1.0
SCALE_SPACING_DJ: Final[float] = 1.0 / 12.0
SMALLEST_SCALE_DAYS: Final[float] = 2.0
MAXIMUM_PERIOD_YEARS: Final[float] = 16.0
SIGNIFICANCE_LEVEL: Final[float] = 0.95
AR1_MINIMUM: Final[float] = -0.99
AR1_MAXIMUM: Final[float] = 0.99
POWER_FLOOR: Final[float] = 1.0e-15
MINIMUM_SCALES_PER_BAND: Final[int] = 2
MINIMUM_RECORDS: Final[int] = 3 * 365
MAXIMUM_INTERPOLATION_GAP_DAYS: Final[int] = 7

REMOVE_LINEAR_TREND: Final[bool] = True
STANDARDIZE_SERIES: Final[bool] = True
ALLOW_SHORT_GAP_INTERPOLATION: Final[bool] = True
LONGITUDE_UNWRAP_BEFORE_ANALYSIS: Final[bool] = True
RETAIN_CHECK_QC_RECORDS: Final[bool] = True

GENERATE_WAVELET_FIGURES: Final[dict[str, bool]] = {
    "traditional_complete": True,
    "bands_multipanel": True,
    "individual_band_figures": True,
    "all_bands_overlay": True,
    "series_scalogram_bands_global": True,
    "series_and_scalogram": True,
    "global_spectrum_only": True,
    "real_part_only": True,
}

BAND_OVERLAY_NORMALIZATION: Final[BandNormalization] = "standard_score"
BAND_MULTIPANEL_SHARE_Y: Final[bool] = False
SHOW_BAND_SIGNIFICANCE: Final[bool] = True
SHADE_SIGNIFICANT_BAND_POWER: Final[bool] = True
SHOW_BAND_MEAN_LINE: Final[bool] = True
SHOW_BAND_PERCENTILE_LINES: Final[bool] = False
BAND_PERCENTILES: Final[tuple[float, ...]] = (0.90, 0.95)

SCALOGRAM_PERIOD_UNIT: Final[PeriodUnit] = "years"
SCALOGRAM_SHOW_BAND_LIMITS: Final[bool] = True
SCALOGRAM_SHOW_SIGNIFICANCE_CONTOUR: Final[bool] = True
SCALOGRAM_SHOW_CONE_OF_INFLUENCE: Final[bool] = True
SCALOGRAM_POWER_LEVELS: Final[int] = 24
REAL_PART_LEVELS: Final[int] = 25
REAL_PART_PERCENTILE_LIMIT: Final[float] = 99.0

GLOBAL_SPECTRUM_PERIOD_UNIT: Final[PeriodUnit] = "years"
GLOBAL_SPECTRUM_LOG2_X_AXIS: Final[bool] = True
GLOBAL_SPECTRUM_SHOW_BANDS: Final[bool] = True
GLOBAL_SPECTRUM_SHADE_BANDS: Final[bool] = True
GLOBAL_SPECTRUM_SHOW_BACKGROUND: Final[bool] = True
GLOBAL_SPECTRUM_SHOW_SIGNIFICANCE: Final[bool] = True
GLOBAL_SPECTRUM_SHOW_PEAKS: Final[bool] = True
GLOBAL_SPECTRUM_MAXIMUM_ANNOTATED_PEAKS: Final[int] = 5

EXPORT_PNG: Final[bool] = True
EXPORT_PDF: Final[bool] = True
EXPORT_SVG: Final[bool] = False
EXPORT_COMPLEX_COEFFICIENTS: Final[bool] = True
EXPORT_REAL_COEFFICIENTS: Final[bool] = True
EXPORT_IMAGINARY_COEFFICIENTS: Final[bool] = True
EXPORT_PHASE: Final[bool] = True
EXPORT_LOCAL_POWER: Final[bool] = True
EXPORT_FULL_BAND_TIME_SERIES: Final[bool] = True
EXPORT_GLOBAL_SPECTRUM_TABLE: Final[bool] = True
EXPORT_SUMMARY_TABLES: Final[bool] = True


def enabled_wavelet_bands() -> tuple[WaveletBand, ...]:
    """Return enabled bands in publication order."""
    return tuple(sorted(
        (band for band in WAVELET_BANDS if band.enabled),
        key=lambda band: (band.order, band.key),
    ))


def expected_figure_count_per_variable() -> int:
    """Return expected figure count for one variable and one output format."""
    count = 0
    bands = len(enabled_wavelet_bands())
    count += int(GENERATE_WAVELET_FIGURES["traditional_complete"])
    count += int(GENERATE_WAVELET_FIGURES["bands_multipanel"])
    count += (
        bands
        if GENERATE_WAVELET_FIGURES["individual_band_figures"]
        else 0
    )
    count += int(GENERATE_WAVELET_FIGURES["all_bands_overlay"])
    count += int(GENERATE_WAVELET_FIGURES["series_scalogram_bands_global"])
    count += int(GENERATE_WAVELET_FIGURES["series_and_scalogram"])
    count += int(GENERATE_WAVELET_FIGURES["global_spectrum_only"])
    count += int(GENERATE_WAVELET_FIGURES["real_part_only"])
    return count


def validate_wavelet_configuration() -> tuple[str, ...]:
    """Validate scientific and graphical wavelet configuration."""
    errors: list[str] = []
    warnings: list[str] = []
    bands = enabled_wavelet_bands()

    if not 1 <= len(bands) <= 6:
        errors.append("Between one and six wavelet bands must be enabled.")

    for attribute, label in (
        ("key", "keys"),
        ("short_label", "short labels"),
        ("order", "orders"),
    ):
        values = [getattr(band, attribute) for band in bands]
        if len(values) != len(set(values)):
            errors.append(f"Enabled band {label} must be unique.")

    for band in bands:
        if band.period_unit not in ("days", "years"):
            errors.append(f"{band.key}: period unit must be days or years.")
        if band.minimum_period <= 0:
            errors.append(f"{band.key}: minimum period must be positive.")
        if band.maximum_period <= band.minimum_period:
            errors.append(f"{band.key}: maximum must exceed minimum.")
        if band.maximum_period_days > MAXIMUM_PERIOD_YEARS * DAYS_PER_YEAR:
            errors.append(
                f"{band.key}: upper limit exceeds maximum CWT period."
            )

    for index, first in enumerate(bands):
        for second in bands[index + 1:]:
            overlap = (
                min(first.maximum_period_days, second.maximum_period_days)
                - max(first.minimum_period_days, second.minimum_period_days)
            )
            if overlap > 0:
                warnings.append(
                    f"Overlapping bands: {first.key} and {second.key}."
                )

    if SAMPLING_INTERVAL_DAYS <= 0:
        errors.append("SAMPLING_INTERVAL_DAYS must be positive.")
    if SCALE_SPACING_DJ <= 0:
        errors.append("SCALE_SPACING_DJ must be positive.")
    if MORLET_OMEGA0 <= 0:
        errors.append("MORLET_OMEGA0 must be positive.")
    if not 0 < SIGNIFICANCE_LEVEL < 1:
        errors.append("SIGNIFICANCE_LEVEL must be between zero and one.")
    if SCALOGRAM_PERIOD_UNIT not in ("days", "years"):
        errors.append("Invalid SCALOGRAM_PERIOD_UNIT.")
    if GLOBAL_SPECTRUM_PERIOD_UNIT not in ("days", "years"):
        errors.append("Invalid GLOBAL_SPECTRUM_PERIOD_UNIT.")
    if BAND_OVERLAY_NORMALIZATION not in (
        "absolute", "mean", "standard_score", "maximum"
    ):
        errors.append("Invalid BAND_OVERLAY_NORMALIZATION.")

    if errors:
        raise ValueError(
            "Wavelet configuration validation failed:\n"
            + "\n".join(f"  - {error}" for error in errors)
        )

    messages = [
        f"Enabled wavelet bands: {len(bands)}",
        (
            "Expected figures per variable and output format: "
            f"{expected_figure_count_per_variable()}"
        ),
        "Global spectrum: Fourier period on X; power on Y.",
        (
            "Global spectrum X scale: "
            f"{'logarithmic base 2' if GLOBAL_SPECTRUM_LOG2_X_AXIS else 'linear'}"
        ),
    ]
    messages.extend(f"WARNING: {warning}" for warning in warnings)
    return tuple(messages)


def wavelet_configuration_summary_lines() -> tuple[str, ...]:
    """Return a human-readable wavelet configuration summary."""
    lines = [
        f"Mother wavelet              : {MOTHER_WAVELET_NAME}",
        f"Morlet omega0               : {MORLET_OMEGA0}",
        f"Sampling interval           : {SAMPLING_INTERVAL_DAYS:.6f} day",
        f"Scale spacing dj            : {SCALE_SPACING_DJ:.12f}",
        f"Maximum Fourier period      : {MAXIMUM_PERIOD_YEARS:.3f} years",
        f"Significance level          : {SIGNIFICANCE_LEVEL:.3f}",
        f"Enabled spectral bands      : {len(enabled_wavelet_bands())}",
        (
            "Expected figures/variable   : "
            f"{expected_figure_count_per_variable()}"
        ),
        (
            "Global spectrum X axis      : "
            f"Fourier period ({GLOBAL_SPECTRUM_PERIOD_UNIT})"
        ),
        "Global spectrum Y axis      : Global wavelet power",
    ]
    for index, band in enumerate(enabled_wavelet_bands(), start=1):
        lines.append(
            f"Band {index:<22d}: {band.label} [{band.period_range_text}]"
        )
    return tuple(lines)


def main() -> None:
    """Validate and print the wavelet configuration."""
    messages = validate_wavelet_configuration()
    print("=" * 78)
    print("PACIFIC WARM POOL WAVELET CONFIGURATION")
    print("=" * 78)
    for line in wavelet_configuration_summary_lines():
        print(line)
    print()
    print("VALIDATION")
    print("-" * 78)
    for message in messages:
        print(message)
    print()
    print("=" * 78)
    print("WAVELET CONFIGURATION VALIDATION COMPLETED SUCCESSFULLY.")
    print("=" * 78)


if __name__ == "__main__":
    main()