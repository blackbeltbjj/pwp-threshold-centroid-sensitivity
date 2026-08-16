# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Ocean Spectral Analysis Framework (OSAF)
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

FILE
    config/cross_wavelet_nino_soi_config.py

VERSION
    1.0.0

PURPOSE
    Dedicated scientific configuration for the external climate-coupling branch
    of Program 19:

        PWP properties × weekly Niño-region indices
        PWP properties × monthly Southern Oscillation Index (SOI)

ARCHITECTURAL PRINCIPLE
    The validated internal Program 19 is not modified. This module is a sibling
    analysis branch with its own output namespace:

        cross_wavelet_nino_soi

DATA CONTRACT
    Niño input is produced by Program 21:
        data/processed/climate_indices/nino/weekly_nino_indices.csv

    SOI input is produced by Program 21:
        data/processed/climate_indices/soi/monthly_soi.csv

TEMPORAL POLICY
    PWP daily data are aggregated downward to the native temporal resolution of
    each climate index:

        daily PWP -> weekly PWP  <-> weekly Niño
        daily PWP -> monthly PWP <-> monthly SOI

    Climate indices are never artificially upsampled.

PRIMARY ENSO FIELD
    Niño SSTA is the default scientific field for ENSO coupling. Absolute SST
    is retained as an optional diagnostic mode because it contains substantial
    seasonal-cycle information.

PHASE CONVENTION
    X = PWP variable
    Y = climate index

    XWT = W_X * conjugate(W_Y)
    relative phase = phase(X) - phase(Y)

    Positive equivalent lag means:
        PWP variable leads climate index.

    Negative equivalent lag means:
        climate index leads PWP variable.

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

from typing import Final, Literal


NinoField = Literal["SSTA", "SST"]
AggregationStatistic = Literal["mean", "median"]


# =============================================================================
# DATA PRODUCTS FROM PROGRAM 21
# =============================================================================

NINO_RELATIVE_PATH: Final[tuple[str, ...]] = (
    "data",
    "processed",
    "climate_indices",
    "nino",
    "weekly_nino_indices.csv",
)

SOI_RELATIVE_PATH: Final[tuple[str, ...]] = (
    "data",
    "processed",
    "climate_indices",
    "soi",
    "monthly_soi.csv",
)

OUTPUT_NAMESPACE: Final[str] = "cross_wavelet_nino_soi"


# =============================================================================
# SCIENTIFIC ANALYSIS SELECTION
# =============================================================================

ENABLE_WEEKLY_NINO_ANALYSIS: Final[bool] = True
ENABLE_MONTHLY_SOI_ANALYSIS: Final[bool] = True

NINO_ANALYSIS_FIELD: Final[NinoField] = "SSTA"

PWP_VARIABLE_KEYS: Final[tuple[str, ...]] = (
    "longitude",
    "latitude",
    "area",
    "mean_sst",
)

NINO_REGION_KEYS: Final[tuple[str, ...]] = (
    "nino12",
    "nino3",
    "nino34",
    "nino4",
)

PWP_LABELS: Final[dict[str, str]] = {
    "longitude": "PWP centroid longitude",
    "latitude": "PWP centroid latitude",
    "area": "PWP total area",
    "mean_sst": "Area-weighted mean PWP SST",
}

NINO_LABELS: Final[dict[str, str]] = {
    "nino12": "Niño 1+2",
    "nino3": "Niño 3",
    "nino34": "Niño 3.4",
    "nino4": "Niño 4",
}

SOI_LABEL: Final[str] = "Southern Oscillation Index"


# =============================================================================
# TEMPORAL ALIGNMENT
# =============================================================================

# The Program-21 weekly Niño dates are Wednesdays.
WEEKLY_AGGREGATION_FREQUENCY: Final[str] = "W-WED"
WEEKLY_AGGREGATION_STATISTIC: Final[AggregationStatistic] = "mean"
MINIMUM_VALID_DAYS_PER_WEEK: Final[int] = 5

# Program 21 stores SOI on the first day of each month.
MONTHLY_AGGREGATION_FREQUENCY: Final[str] = "MS"
MONTHLY_AGGREGATION_STATISTIC: Final[AggregationStatistic] = "mean"
MINIMUM_VALID_DAYS_PER_MONTH: Final[int] = 20

# External climate-index products are already QC'd and continuous. Do not fill
# missing index values or fabricate observations.
WEEKLY_PAIR_INTERPOLATION_LIMIT: Final[int] = 0
MONTHLY_PAIR_INTERPOLATION_LIMIT: Final[int] = 0

MINIMUM_COMMON_WEEKLY_RECORDS: Final[int] = 128
MINIMUM_COMMON_MONTHLY_RECORDS: Final[int] = 96


# =============================================================================
# PROGRAM-21 CONTRACT VALIDATION
# =============================================================================

REQUIRE_PROGRAM21_QC_PASS: Final[bool] = True

EXPECTED_NINO_COLUMNS: Final[tuple[str, ...]] = (
    "date",
    "nino12_sst",
    "nino12_ssta",
    "nino3_sst",
    "nino3_ssta",
    "nino34_sst",
    "nino34_ssta",
    "nino4_sst",
    "nino4_ssta",
    "qc_status",
)

EXPECTED_SOI_COLUMNS: Final[tuple[str, ...]] = (
    "date",
    "year",
    "month",
    "month_name",
    "soi",
    "source_year_status",
    "qc_status",
)


# =============================================================================
# CROSS-THRESHOLD PRIORITY FIGURES
# =============================================================================

GENERATE_CROSS_THRESHOLD_PRIORITY_FIGURES: Final[bool] = True

PRIORITY_PAIR_KEYS: Final[tuple[str, ...]] = (
    "longitude_x_nino34_ssta",
    "area_x_nino34_ssta",
    "mean_sst_x_nino34_ssta",
    "longitude_x_soi",
    "area_x_soi",
)


# =============================================================================
# SCIENTIFIC REPORTING
# =============================================================================

GENERATE_COUPLING_EVIDENCE_TABLE: Final[bool] = True

COUPLING_EVIDENCE_COLUMNS: Final[tuple[str, ...]] = (
    "threshold_c",
    "pair_key",
    "pair_label",
    "group",
    "band_key",
    "band_label",
    "mean_cross_power_inside_coi",
    "mean_coherence_inside_coi",
    "fraction_coherence_above_threshold",
    "circular_mean_phase_radians",
    "phase_resultant_length",
    "equivalent_lag_days_positive_x_leads_y",
    "lead_lag_valid_coefficients",
)


def nino_column(region_key: str) -> str:
    """Return the Program-21 column name selected for a Niño region."""

    suffix = NINO_ANALYSIS_FIELD.lower()

    if region_key not in NINO_REGION_KEYS:
        raise KeyError(
            f"Unsupported Niño region key: {region_key}"
        )

    return f"{region_key}_{suffix}"


def validate_cross_wavelet_nino_soi_configuration() -> tuple[str, ...]:
    """Validate the external climate-coupling configuration."""

    errors: list[str] = []

    if NINO_ANALYSIS_FIELD not in ("SSTA", "SST"):
        errors.append(
            "NINO_ANALYSIS_FIELD must be either 'SSTA' or 'SST'."
        )

    if not PWP_VARIABLE_KEYS:
        errors.append(
            "At least one PWP variable must be enabled."
        )

    if not NINO_REGION_KEYS:
        errors.append(
            "At least one Niño region must be enabled."
        )

    if MINIMUM_VALID_DAYS_PER_WEEK < 1:
        errors.append(
            "MINIMUM_VALID_DAYS_PER_WEEK must be positive."
        )

    if MINIMUM_VALID_DAYS_PER_MONTH < 1:
        errors.append(
            "MINIMUM_VALID_DAYS_PER_MONTH must be positive."
        )

    if WEEKLY_PAIR_INTERPOLATION_LIMIT != 0:
        errors.append(
            "Weekly external-index interpolation must remain zero."
        )

    if MONTHLY_PAIR_INTERPOLATION_LIMIT != 0:
        errors.append(
            "Monthly external-index interpolation must remain zero."
        )

    if errors:
        raise ValueError(
            "Cross-wavelet Niño/SOI configuration validation failed:\n"
            + "\n".join(
                f"  - {error}"
                for error in errors
            )
        )

    return (
        f"Niño analysis field          : {NINO_ANALYSIS_FIELD}",
        f"Weekly Niño enabled          : {ENABLE_WEEKLY_NINO_ANALYSIS}",
        f"Monthly SOI enabled          : {ENABLE_MONTHLY_SOI_ANALYSIS}",
        (
            "Weekly aggregation          : "
            f"{WEEKLY_AGGREGATION_FREQUENCY} / "
            f"{WEEKLY_AGGREGATION_STATISTIC}"
        ),
        (
            "Monthly aggregation         : "
            f"{MONTHLY_AGGREGATION_FREQUENCY} / "
            f"{MONTHLY_AGGREGATION_STATISTIC}"
        ),
        "Climate-index upsampling       : prohibited",
        f"Output namespace             : {OUTPUT_NAMESPACE}",
    )


if __name__ == "__main__":
    print("=" * 78)
    print("OSAF CROSS-WAVELET NIÑO/SOI CONFIGURATION")
    print("=" * 78)

    for line in validate_cross_wavelet_nino_soi_configuration():
        print(line)

    print("=" * 78)
    print("CONFIGURATION VALIDATION COMPLETED SUCCESSFULLY.")
    print("=" * 78)