#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
===============================================================================
PROJECT
    Ocean Spectral Analysis Framework (OSAF)
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

PROGRAM
    33_pwp_enso_wavelet_band_modulation.py

TITLE
    ENSO Modulation of Time-Varying PWP Longitude Wavelet-Band Variance

VERSION
    1.0.0

AUTHOR
    Fabio Vieira Machado

PYTHON
    Python 3.10+

===============================================================================
PURPOSE
===============================================================================

    Program 33 is the FINAL planned analytical program for AUDIT 06 before
    scientific freeze and manuscript production.

    It tests whether the time-varying wavelet variance of Pacific Warm Pool
    (PWP) centroid longitude is statistically associated with ENSO state and/or
    ENSO magnitude across the FINAL wavelet bands defined by Program 32.

    Program 33 DOES NOT recompute the Continuous Wavelet Transform.

    It reads the authoritative Program-32 daily band-power product:

        outputs/tables/threshold_comparison/
            pwp_longitude_wavelet_band_variance/
                pwp_longitude_wavelet_band_power_daily.csv

    and the Program-32 band definitions:

        pwp_longitude_wavelet_band_definitions.csv

    It then combines these data, at MONTHLY resolution, with the canonical
    Program-21 ENSO indices:

        Niño 1+2 SSTA
        Niño 3 SSTA
        Niño 3.4 SSTA
        Niño 4 SSTA
        Southern Oscillation Index (SOI)

    The final PWP threshold hierarchy is:

        28.0 °C  = PRIMARY PWP DEFINITION
        28.5 °C  = threshold sensitivity
        29.0 °C  = warm-core / threshold sensitivity

===============================================================================
WHY PROGRAM 33 EXISTS
===============================================================================

    Earlier Program 28 tested ENSO modulation of a legacy annual band:

        330–400 days

    and treated 29 °C as the primary threshold.

    Subsequent Audits 01–05 changed the final scientific architecture:

        1. 28 °C was selected as the primary PWP definition after long-term
           spatial-connectivity and centroid-sensitivity tests.

        2. Program 32 introduced the final, explicitly defined hierarchy of
           annual and interannual wavelet bands.

        3. Program 32 imposed a conservative strict cone-of-influence rule:
           the COMPLETE band must lie inside the reliable COI.

    Program 33 therefore repeats and extends the modulation test using the
    FINAL scientific definitions rather than the legacy Program-28 settings.

===============================================================================
SCIENTIFIC QUESTIONS
===============================================================================

    Q1. ZERO-LAG MODULATION

        Is scale-averaged PWP longitude wavelet power associated
        contemporaneously with ENSO?

        For band B and ENSO predictor E:

            r_B,E(0) = Corr[P_B(t), E(t)]

    Q2. ENSO STATE VERSUS ENSO MAGNITUDE

        SIGNED ENSO STATE:

            E(t)

        asks whether band power systematically increases or decreases toward
        one ENSO phase.

        ENSO MAGNITUDE:

            |E(t)|

        asks whether stronger ENSO departures, irrespective of warm/cold sign,
        are associated with stronger PWP band variance.

        These are different hypotheses and are NEVER merged into one score.

    Q3. LAGGED MODULATION

        Is band power associated with an earlier/later ENSO state?

            r_B,E(τ) = Corr[P_B(t + τ), E(t)]

        Positive lag:

            τ > 0

        means:

            ENSO predictor LEADS PWP wavelet-band power.

        Example:

            lag = +8 months

        means ENSO(t) is correlated with PWP band power at t + 8 months.

    Q4. BAND DEPENDENCE

        Does association differ between:

            annual
            quasi-biennial
            lower interannual
            ~4-year interannual
            integrated interannual
            broad ENSO scales?

    Q5. THRESHOLD ROBUSTNESS

        Do effect sign, magnitude and lag structure persist across
        28.0, 28.5 and 29.0 °C PWP definitions?

===============================================================================
FINAL PROGRAM-32 WAVELET BANDS
===============================================================================

    Program 33 validates these values against the Program-32 definition table
    before analysis.

    1. ANNUAL BROAD
        280–448 days
        40–64 weeks

    2. ANNUAL STRICT
        336–378 days
        48–54 weeks

    3. QUASI-BIENNIAL
        550–1096 days
        approximately 1.5–3.0 years

    4. INTERANNUAL 1
        805–980 days
        approximately 2.2–2.7 years

    5. INTERANNUAL 2
        1280–1642 days
        approximately 3.5–4.5 years

    6. INTERANNUAL 3
        805–1642 days
        approximately 2.2–4.5 years

    7. ENSO BROAD
        1095–2922 days
        approximately 3–8 years

    IMPORTANT:
        These bands overlap intentionally. They are diagnostic views of the
        time-frequency field, not mutually exclusive variance partitions.

===============================================================================
PROGRAM-32 BAND-POWER EQUATION
===============================================================================

    Program 32 reconstructed scale-averaged physical wavelet power as:

                               j2
                               ----
                               \     |W_t(s_j)|²
        P_B(t) = Δj Δt / Cδ     >    -----------
                               /         s_j
                               ----
                               j=j1

    where:

        P_B(t)
            time-varying wavelet variance within band B;

        |W_t(s_j)|²
            physical Morlet wavelet power;

        s_j
            wavelet scale;

        Δj
            logarithmic scale spacing;

        Δt
            sampling interval;

        Cδ
            Morlet reconstruction factor.

    Program 33 reads P_B(t) from Program 32 and never changes its wavelet
    definition.

===============================================================================
CONE OF INFLUENCE POLICY
===============================================================================

    Program 32 provides:

        band_power_strict_coi
        strict_coi_valid

    Program 33 uses ONLY:

        strict_coi_valid == True
        AND finite band_power_strict_coi

    for monthly aggregation.

    Thus edge-contaminated daily band power does not enter the main analysis.

    This is especially important for the broad ENSO band, where long periods
    necessarily reduce the usable record near both endpoints.

===============================================================================
COMMON TEMPORAL RESOLUTION
===============================================================================

    EVERYTHING IS ANALYZED MONTHLY.

    PWP wavelet band power:
        daily strict-COI-valid observations
            -> monthly mean / median / count

    Niño 1+2 / Niño 3 / Niño 3.4 / Niño 4:
        Program-21 weekly SSTA
            -> monthly mean

    SOI:
        Program-21 monthly
            -> unchanged

    NO climate-index upsampling is permitted.

    A Program-32 band-power month is retained only when at least:

        15 strict-COI-valid daily observations

    contribute to that month.

===============================================================================
STATISTICAL METHOD
===============================================================================

    For every:

        threshold × band × ENSO index × test type

    Program 33 calculates at zero lag:

        Pearson correlation
        Spearman rank correlation
        AR(1)-adjusted effective sample size
        autocorrelation-adjusted p-value

    AR(1) effective sample size:

                       1 - ρx ρy
        N_eff = N -------------------
                       1 + ρx ρy

    where ρx and ρy are lag-1 autocorrelations of the paired series.

    The correlation t statistic is:

                     r sqrt(N_eff - 2)
        t = --------------------------------
                sqrt(1 - r²)

    with approximately:

        df = N_eff - 2

===============================================================================
ZERO-LAG MULTIPLE-TEST CONTROL
===============================================================================

    Because Program 33 analyzes more bands than Program 28, it provides TWO
    explicit Benjamini-Hochberg FDR controls.

    A. WITHIN-BAND FDR

        For each:
            threshold × band × test type

        FDR is applied across the five ENSO indices.

        This reproduces the logic of Program 28 for each independently defined
        band.

    B. GLOBAL-BAND FDR

        For each:
            threshold × test type

        FDR is applied across all:

            7 bands × 5 ENSO indices = 35 relationships.

        This is the more conservative exploratory-family control.

    Program 33 reports BOTH values.

    The primary hypothesis summary uses GLOBAL-BAND FDR for zero-lag support.

===============================================================================
LAGGED ANALYSIS
===============================================================================

    Lag window:

        -24 to +24 calendar months

    IMPORTANT:
        Lag alignment is performed using CALENDAR MONTHS, not by blindly
        shifting array positions.

        This prevents a missing month from being interpreted as a one-month
        lag.

    For each individual:

        threshold × band × predictor × test type

    Benjamini-Hochberg FDR is applied across the 49 tested lags.

    Lagged tests are treated as exploratory temporal diagnostics because the
    best lag is selected after examining the lag family.

    A significant lag is NOT interpreted as proof of causal direction.

===============================================================================
PRIMARY VERSUS SECONDARY HYPOTHESES
===============================================================================

    PRIMARY annual-cycle hypotheses:

        Annual Broad
        Annual Strict

    They test whether ENSO is associated with the intensity of the annual
    zonal cycle under two independently specified annual-band definitions.

    SECONDARY scale-specific hypotheses:

        Quasi-biennial
        Interannual 1
        Interannual 2
        Interannual 3
        ENSO Broad

    They test how ENSO association changes across interannual scales.

===============================================================================
SCIENTIFIC INTERPRETATION POLICY
===============================================================================

    Allowed language:

        association
        covariation
        temporal correspondence
        lagged association
        modulation supported statistically
        modulation not supported
        threshold robustness

    Prohibited inference from Program 33 alone:

        ENSO causes ...
        ENSO forces ...
        PWP causes ENSO ...
        lag proves propagation ...
        coherence/correlation proves mechanism ...

    A non-significant result is retained as a scientific result.

===============================================================================
AUTHORITATIVE INPUTS
===============================================================================

    1. Program 32 daily strict-COI band power

       outputs/tables/threshold_comparison/
           pwp_longitude_wavelet_band_variance/
               pwp_longitude_wavelet_band_power_daily.csv

    2. Program 32 band definitions

       outputs/tables/threshold_comparison/
           pwp_longitude_wavelet_band_variance/
               pwp_longitude_wavelet_band_definitions.csv

    3. Program 21 weekly Niño indices

       Resolved through:
           config.cross_wavelet_nino_soi_config.NINO_RELATIVE_PATH

    4. Program 21 monthly SOI

       Resolved through:
           config.cross_wavelet_nino_soi_config.SOI_RELATIVE_PATH

===============================================================================
OUTPUT DIRECTORIES
===============================================================================

    TABLES

        outputs/tables/threshold_comparison/
            pwp_enso_wavelet_band_modulation/

    FIGURES

        outputs/figures/threshold_comparison/
            pwp_enso_wavelet_band_modulation/

    REPORTS

        outputs/reports/threshold_comparison/
            pwp_enso_wavelet_band_modulation/

===============================================================================
OUTPUT TABLES
===============================================================================

    1. pwp_enso_wavelet_band_monthly_common_data.csv

       Monthly long-format common dataset.

    2. pwp_enso_wavelet_band_zero_lag_association.csv

       Zero-lag Pearson/Spearman results with:
            within-band FDR
            global-band FDR

    3. pwp_enso_wavelet_band_lagged_correlation_curves.csv

       Complete -24..+24 month Pearson lag curves.

    4. pwp_enso_wavelet_band_lagged_correlation_summary.csv

       Strongest absolute lag and strongest FDR-significant lag.

    5. pwp_enso_wavelet_band_threshold_robustness.csv

       Cross-threshold robustness for each band/predictor/test type.

    6. pwp_enso_wavelet_band_hypothesis_summary.csv

       Primary 28 °C band-level hypothesis decisions.

===============================================================================
OUTPUT FIGURES
===============================================================================

    1. pwp_enso_wavelet_band_zero_lag_heatmaps_28C.png/pdf

       A) signed ENSO state
       B) ENSO magnitude

       Rows:
           seven Program-32 bands

       Columns:
           Niño 1+2, Niño 3, Niño 3.4, Niño 4, SOI

       Cell values:
           Pearson r

       A marker identifies relationships surviving global-band FDR.

    2. pwp_enso_wavelet_key_band_nino34_lag_curves_28C.png/pdf

       Four key bands × two tests:

           Annual Broad
           Annual Strict
           Interannual 3
           ENSO Broad

       Left:
           signed Niño 3.4 state

       Right:
           |Niño 3.4|

       FDR-significant lags are explicitly marked.

    3. pwp_enso_wavelet_annual_band_time_series_28C.png/pdf

       A) Annual Broad power
       B) Annual Strict power
       C) Niño 3.4 SSTA
       D) |Niño 3.4 SSTA|

===============================================================================
REPORTS
===============================================================================

    PROGRAM33_PWP_ENSO_WAVELET_BAND_MODULATION.txt
    PROGRAM33_PWP_ENSO_WAVELET_BAND_MODULATION.json

    The report explicitly distinguishes:

        displacement coupling established by earlier programs
        versus
        wavelet-band-power modulation tested here.

===============================================================================
DEPENDENCIES
===============================================================================

    Standard library:
        json
        math
        sys
        datetime
        pathlib
        typing

    Third-party:
        numpy
        pandas
        scipy
        matplotlib

    Project:
        config.config
        config.cross_wavelet_nino_soi_config

===============================================================================
EXECUTION
===============================================================================

    Recommended from project root:

        python src/threshold_analysis/33_pwp_enso_wavelet_band_modulation.py

===============================================================================
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# =============================================================================
# PROJECT DISCOVERY
# =============================================================================

SCRIPT_FILE = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_FILE.parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# CANONICAL PROJECT CONFIGURATION
# =============================================================================

try:
    from config.config import (  # noqa: E402
        FIGURE_DPI,
        GRID_ALPHA,
        GRID_LINESTYLE,
        GRID_LINEWIDTH,
        PROJECT_DIR,
        SAVE_BBOX,
        SAVE_PAD_INCHES,
        SAVE_TRANSPARENT,
        validate_project_configuration,
    )

    from config.cross_wavelet_nino_soi_config import (  # noqa: E402
        EXPECTED_NINO_COLUMNS,
        EXPECTED_SOI_COLUMNS,
        NINO_ANALYSIS_FIELD,
        NINO_LABELS,
        NINO_REGION_KEYS,
        NINO_RELATIVE_PATH,
        SOI_LABEL,
        SOI_RELATIVE_PATH,
        nino_column,
        validate_cross_wavelet_nino_soi_configuration,
    )
except Exception as error:
    raise RuntimeError(
        "Unable to import canonical OSAF configuration for Program 33.\n"
        f"Expected project root: {PROJECT_ROOT}\n"
        f"Original error: {type(error).__name__}: {error}"
    ) from error


# =============================================================================
# PROGRAM CONSTANTS
# =============================================================================

PROGRAM_NAME = "OSAF PROGRAM 33 — PWP–ENSO WAVELET-BAND MODULATION"
PROGRAM_VERSION = "1.0.0"

THRESHOLDS_C = (28.0, 28.5, 29.0)
PRIMARY_THRESHOLD_C = 28.0

MONTHLY_MAX_LAG = 24
MINIMUM_MONTHLY_RECORDS = 120
MINIMUM_VALID_DAYS_PER_MONTH = 15
FDR_ALPHA = 0.05

YEAR_TICK_ANCHOR = 1982
YEAR_TICK_INTERVAL = 4

PROGRAM32_TABLE_DIR = (
    Path(PROJECT_DIR)
    / "outputs"
    / "tables"
    / "threshold_comparison"
    / "pwp_longitude_wavelet_band_variance"
)

PROGRAM32_DAILY_FILE = (
    PROGRAM32_TABLE_DIR
    / "pwp_longitude_wavelet_band_power_daily.csv"
)

PROGRAM32_BAND_FILE = (
    PROGRAM32_TABLE_DIR
    / "pwp_longitude_wavelet_band_definitions.csv"
)

WEEKLY_NINO_FILE = (
    Path(PROJECT_DIR)
    .joinpath(*NINO_RELATIVE_PATH)
)

MONTHLY_SOI_FILE = (
    Path(PROJECT_DIR)
    .joinpath(*SOI_RELATIVE_PATH)
)

TABLE_DIR = (
    Path(PROJECT_DIR)
    / "outputs"
    / "tables"
    / "threshold_comparison"
    / "pwp_enso_wavelet_band_modulation"
)

FIGURE_DIR = (
    Path(PROJECT_DIR)
    / "outputs"
    / "figures"
    / "threshold_comparison"
    / "pwp_enso_wavelet_band_modulation"
)

REPORT_DIR = (
    Path(PROJECT_DIR)
    / "outputs"
    / "reports"
    / "threshold_comparison"
    / "pwp_enso_wavelet_band_modulation"
)

MONTHLY_DATA_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_monthly_common_data.csv"
)

ZERO_LAG_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_zero_lag_association.csv"
)

LAG_CURVES_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_lagged_correlation_curves.csv"
)

LAG_SUMMARY_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_lagged_correlation_summary.csv"
)

ROBUSTNESS_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_threshold_robustness.csv"
)

HYPOTHESIS_FILE = (
    TABLE_DIR
    / "pwp_enso_wavelet_band_hypothesis_summary.csv"
)

REPORT_TXT = (
    REPORT_DIR
    / "PROGRAM33_PWP_ENSO_WAVELET_BAND_MODULATION.txt"
)

REPORT_JSON = (
    REPORT_DIR
    / "PROGRAM33_PWP_ENSO_WAVELET_BAND_MODULATION.json"
)

EXPECTED_BANDS = (
    {
        "band_key": "annual_broad",
        "band_label": "Annual Broad",
        "minimum_period_days": 280.0,
        "maximum_period_days": 448.0,
        "hypothesis_family": "primary_annual",
    },
    {
        "band_key": "annual_strict",
        "band_label": "Annual Strict",
        "minimum_period_days": 336.0,
        "maximum_period_days": 378.0,
        "hypothesis_family": "primary_annual",
    },
    {
        "band_key": "quasi_biennial",
        "band_label": "Quasi-biennial",
        "minimum_period_days": 550.0,
        "maximum_period_days": 1096.0,
        "hypothesis_family": "secondary_interannual",
    },
    {
        "band_key": "interannual_1",
        "band_label": "Interannual 1",
        "minimum_period_days": 805.0,
        "maximum_period_days": 980.0,
        "hypothesis_family": "secondary_interannual",
    },
    {
        "band_key": "interannual_2",
        "band_label": "Interannual 2",
        "minimum_period_days": 1280.0,
        "maximum_period_days": 1642.0,
        "hypothesis_family": "secondary_interannual",
    },
    {
        "band_key": "interannual_3",
        "band_label": "Interannual 3",
        "minimum_period_days": 805.0,
        "maximum_period_days": 1642.0,
        "hypothesis_family": "secondary_interannual",
    },
    {
        "band_key": "enso_broad",
        "band_label": "ENSO Broad",
        "minimum_period_days": 1095.0,
        "maximum_period_days": 2922.0,
        "hypothesis_family": "secondary_interannual",
    },
)

BAND_ORDER = {
    row["band_key"]: i
    for i, row in enumerate(EXPECTED_BANDS)
}

INDEX_ORDER = {
    "nino12": 0,
    "nino3": 1,
    "nino34": 2,
    "nino4": 3,
    "soi": 4,
}


# =============================================================================
# TERMINAL / FIGURE UTILITIES
# =============================================================================

def rule(
    character: str = "=",
    width: int = 78,
) -> None:
    print(character * width)


def section(
    title: str,
) -> None:
    print()
    print(title)
    rule("-")


def item(
    label: str,
    value: object,
) -> None:
    print(
        f"{label:<43s}: {value}"
    )


def threshold_equal(
    values: pd.Series | np.ndarray,
    threshold_c: float,
) -> np.ndarray:
    return np.isclose(
        np.asarray(
            values,
            dtype=float,
        ),
        float(threshold_c),
    )


def add_panel_label(
    axis: plt.Axes,
    label: str,
) -> None:
    axis.text(
        0.01,
        0.95,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        fontweight="bold",
    )


def configure_year_axis(
    axis: plt.Axes,
    minimum_date: pd.Timestamp,
    maximum_date: pd.Timestamp,
    show_labels: bool,
) -> None:
    first_year = max(
        YEAR_TICK_ANCHOR,
        int(
            minimum_date.year
        ),
    )

    offset = (
        first_year
        - YEAR_TICK_ANCHOR
    ) % YEAR_TICK_INTERVAL

    if offset:
        first_year += (
            YEAR_TICK_INTERVAL
            - offset
        )

    years = np.arange(
        first_year,
        int(
            maximum_date.year
        ) + 1,
        YEAR_TICK_INTERVAL,
        dtype=int,
    )

    ticks = pd.DatetimeIndex(
        pd.to_datetime(
            [
                f"{year}-01-01"
                for year in years
            ]
        )
    )

    inside = (
        (ticks >= minimum_date)
        & (ticks <= maximum_date)
    )

    ticks = ticks[inside]
    years = years[inside]

    axis.set_xlim(
        minimum_date,
        maximum_date,
    )
    axis.margins(
        x=0.0
    )
    axis.set_xticks(
        ticks
    )

    if show_labels:
        axis.set_xticklabels(
            [
                str(year)
                for year in years
            ]
        )
        axis.set_xlabel(
            "Year"
        )
    else:
        axis.tick_params(
            axis="x",
            labelbottom=False,
        )

    axis.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )


def save_figure(
    figure: plt.Figure,
    stem: str,
) -> tuple[Path, Path]:
    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    png = (
        FIGURE_DIR
        / f"{stem}.png"
    )
    pdf = (
        FIGURE_DIR
        / f"{stem}.pdf"
    )

    figure.savefig(
        png,
        dpi=FIGURE_DPI,
        bbox_inches=SAVE_BBOX,
        pad_inches=SAVE_PAD_INCHES,
        transparent=SAVE_TRANSPARENT,
    )

    figure.savefig(
        pdf,
        bbox_inches=SAVE_BBOX,
        pad_inches=SAVE_PAD_INCHES,
        transparent=SAVE_TRANSPARENT,
    )

    plt.close(
        figure
    )

    return (
        png,
        pdf,
    )


# =============================================================================
# STATISTICAL UTILITIES
# =============================================================================

def finite_pairs(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    x = np.asarray(
        x,
        dtype=float,
    )
    y = np.asarray(
        y,
        dtype=float,
    )

    valid = (
        np.isfinite(x)
        & np.isfinite(y)
    )

    return (
        x[valid],
        y[valid],
    )


def lag1_autocorrelation(
    values: np.ndarray,
) -> float:
    values = np.asarray(
        values,
        dtype=float,
    )

    if values.size < 3:
        return np.nan

    x, y = finite_pairs(
        values[:-1],
        values[1:],
    )

    if x.size < 3:
        return np.nan

    if (
        np.nanstd(x) <= 0.0
        or np.nanstd(y) <= 0.0
    ):
        return np.nan

    return float(
        np.corrcoef(
            x,
            y,
        )[0, 1]
    )


def effective_sample_size(
    x: np.ndarray,
    y: np.ndarray,
) -> float:
    x, y = finite_pairs(
        x,
        y,
    )

    n = x.size

    if n < 3:
        return np.nan

    r1_x = lag1_autocorrelation(
        x
    )
    r1_y = lag1_autocorrelation(
        y
    )

    if (
        not np.isfinite(r1_x)
        or not np.isfinite(r1_y)
    ):
        return float(n)

    product = float(
        r1_x
        * r1_y
    )

    denominator = (
        1.0
        + product
    )

    if denominator <= 0.0:
        return float(n)

    n_eff = (
        n
        * (
            1.0
            - product
        )
        / denominator
    )

    return float(
        np.clip(
            n_eff,
            3.0,
            float(n),
        )
    )


def correlation_p_from_neff(
    r: float,
    n_eff: float,
) -> float:
    if (
        not np.isfinite(r)
        or not np.isfinite(n_eff)
        or n_eff <= 2.0
    ):
        return np.nan

    r_clipped = float(
        np.clip(
            r,
            -0.999999999999,
            0.999999999999,
        )
    )

    t_value = (
        r_clipped
        * math.sqrt(
            (
                n_eff
                - 2.0
            )
            / max(
                1.0e-15,
                1.0
                - r_clipped**2,
            )
        )
    )

    return float(
        2.0
        * stats.t.sf(
            abs(t_value),
            df=(
                n_eff
                - 2.0
            ),
        )
    )


def benjamini_hochberg(
    p_values: np.ndarray,
) -> np.ndarray:
    p_values = np.asarray(
        p_values,
        dtype=float,
    )

    q_values = np.full(
        p_values.shape,
        np.nan,
        dtype=float,
    )

    valid = np.isfinite(
        p_values
    )

    if not np.any(valid):
        return q_values

    p = p_values[valid]

    order = np.argsort(
        p
    )
    ranked = p[order]

    m = len(
        ranked
    )

    adjusted = (
        ranked
        * float(m)
        / np.arange(
            1,
            m + 1,
            dtype=float,
        )
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.clip(
        adjusted,
        0.0,
        1.0,
    )

    restored = np.empty_like(
        adjusted
    )
    restored[order] = adjusted

    q_values[valid] = restored

    return q_values


def correlation_metrics(
    x: np.ndarray,
    y: np.ndarray,
) -> dict[
    str,
    float,
]:
    x, y = finite_pairs(
        x,
        y,
    )

    if x.size < 3:
        return {
            "n": float(
                x.size
            ),
            "n_eff": np.nan,
            "pearson_r": np.nan,
            "pearson_p_eff": np.nan,
            "spearman_rho": np.nan,
            "spearman_p_eff": np.nan,
        }

    if (
        np.nanstd(x) <= 0.0
        or np.nanstd(y) <= 0.0
    ):
        return {
            "n": float(
                x.size
            ),
            "n_eff": np.nan,
            "pearson_r": np.nan,
            "pearson_p_eff": np.nan,
            "spearman_rho": np.nan,
            "spearman_p_eff": np.nan,
        }

    pearson_r = float(
        stats.pearsonr(
            x,
            y,
        ).statistic
    )

    spearman_rho = float(
        stats.spearmanr(
            x,
            y,
        ).statistic
    )

    n_eff = effective_sample_size(
        x,
        y,
    )

    return {
        "n": float(
            x.size
        ),
        "n_eff": n_eff,
        "pearson_r": pearson_r,
        "pearson_p_eff": correlation_p_from_neff(
            pearson_r,
            n_eff,
        ),
        "spearman_rho": spearman_rho,
        "spearman_p_eff": correlation_p_from_neff(
            spearman_rho,
            n_eff,
        ),
    }


# =============================================================================
# CONFIGURATION / INPUT VALIDATION
# =============================================================================

def validate_configuration() -> None:
    if (
        Path(
            PROJECT_DIR
        ).resolve()
        != PROJECT_ROOT.resolve()
    ):
        raise ValueError(
            "Project-root mismatch:\n"
            f"  config.PROJECT_DIR : {PROJECT_DIR}\n"
            f"  script root        : {PROJECT_ROOT}"
        )

    validate_project_configuration()
    validate_cross_wavelet_nino_soi_configuration()

    for directory in (
        TABLE_DIR,
        FIGURE_DIR,
        REPORT_DIR,
    ):
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    for source in (
        PROGRAM32_DAILY_FILE,
        PROGRAM32_BAND_FILE,
        WEEKLY_NINO_FILE,
        MONTHLY_SOI_FILE,
    ):
        if not source.is_file():
            raise FileNotFoundError(
                "Required Program-33 input is missing:\n"
                f"{source}"
            )


def load_and_validate_band_definitions() -> pd.DataFrame:
    data = pd.read_csv(
        PROGRAM32_BAND_FILE,
        low_memory=False,
    )

    required = {
        "band_key",
        "band_label",
        "minimum_period_days",
        "maximum_period_days",
    }

    missing = sorted(
        required.difference(
            data.columns
        )
    )

    if missing:
        raise KeyError(
            "Program-32 band-definition table is missing columns: "
            f"{missing}"
        )

    output = data.copy()

    output[
        "minimum_period_days"
    ] = pd.to_numeric(
        output[
            "minimum_period_days"
        ],
        errors="raise",
    )

    output[
        "maximum_period_days"
    ] = pd.to_numeric(
        output[
            "maximum_period_days"
        ],
        errors="raise",
    )

    if output[
        "band_key"
    ].duplicated().any():
        raise ValueError(
            "Duplicate Program-32 band keys."
        )

    expected_keys = {
        row["band_key"]
        for row in EXPECTED_BANDS
    }

    observed_keys = set(
        output[
            "band_key"
        ].astype(str)
    )

    if observed_keys != expected_keys:
        raise ValueError(
            "Program-32 band keys do not match the frozen Program-33 set.\n"
            f"Expected: {sorted(expected_keys)}\n"
            f"Observed: {sorted(observed_keys)}"
        )

    for expected in EXPECTED_BANDS:
        row = output.loc[
            output[
                "band_key"
            ].eq(
                expected[
                    "band_key"
                ]
            )
        ].iloc[0]

        if not np.isclose(
            float(
                row[
                    "minimum_period_days"
                ]
            ),
            float(
                expected[
                    "minimum_period_days"
                ]
            ),
        ):
            raise ValueError(
                "Program-32 minimum band period changed for "
                f"{expected['band_key']}."
            )

        if not np.isclose(
            float(
                row[
                    "maximum_period_days"
                ]
            ),
            float(
                expected[
                    "maximum_period_days"
                ]
            ),
        ):
            raise ValueError(
                "Program-32 maximum band period changed for "
                f"{expected['band_key']}."
            )

    output[
        "band_order"
    ] = output[
        "band_key"
    ].map(
        BAND_ORDER
    )

    output[
        "hypothesis_family"
    ] = output[
        "band_key"
    ].map(
        {
            row["band_key"]:
            row["hypothesis_family"]
            for row in EXPECTED_BANDS
        }
    )

    return (
        output.sort_values(
            "band_order"
        )
        .reset_index(
            drop=True
        )
    )


def parse_bool_series(
    series: pd.Series,
) -> pd.Series:
    if pd.api.types.is_bool_dtype(
        series
    ):
        return series.astype(
            bool
        )

    return (
        series
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(
            (
                "true",
                "1",
                "yes",
            )
        )
    )


def load_program32_daily_power() -> pd.DataFrame:
    data = pd.read_csv(
        PROGRAM32_DAILY_FILE,
        low_memory=False,
    )

    required = {
        "date",
        "threshold_c",
        "band_key",
        "band_label",
        "minimum_period_days",
        "maximum_period_days",
        "band_power_strict_coi",
        "strict_coi_valid",
    }

    missing = sorted(
        required.difference(
            data.columns
        )
    )

    if missing:
        raise KeyError(
            "Program-32 daily band-power table is missing columns: "
            f"{missing}"
        )

    output = data.copy()

    output[
        "date"
    ] = pd.to_datetime(
        output[
            "date"
        ],
        errors="raise",
    )

    output[
        "threshold_c"
    ] = pd.to_numeric(
        output[
            "threshold_c"
        ],
        errors="raise",
    )

    output[
        "band_power_strict_coi"
    ] = pd.to_numeric(
        output[
            "band_power_strict_coi"
        ],
        errors="coerce",
    )

    output[
        "strict_coi_valid"
    ] = parse_bool_series(
        output[
            "strict_coi_valid"
        ]
    )

    allowed_thresholds = np.array(
        THRESHOLDS_C,
        dtype=float,
    )

    threshold_ok = np.any(
        np.isclose(
            output[
                "threshold_c"
            ].to_numpy(
                dtype=float
            )[:, np.newaxis],
            allowed_thresholds[
                np.newaxis,
                :
            ],
        ),
        axis=1,
    )

    if not np.all(
        threshold_ok
    ):
        invalid = sorted(
            set(
                output.loc[
                    ~threshold_ok,
                    "threshold_c",
                ].astype(
                    float
                )
            )
        )

        raise ValueError(
            "Unexpected thresholds in Program-32 daily table: "
            f"{invalid}"
        )

    expected_keys = set(
        BAND_ORDER
    )

    observed_keys = set(
        output[
            "band_key"
        ].astype(str)
    )

    if observed_keys != expected_keys:
        raise ValueError(
            "Program-32 daily band keys do not match the frozen set.\n"
            f"Expected: {sorted(expected_keys)}\n"
            f"Observed: {sorted(observed_keys)}"
        )

    duplicates = output.duplicated(
        subset=[
            "date",
            "threshold_c",
            "band_key",
        ]
    )

    if duplicates.any():
        raise ValueError(
            "Duplicate date × threshold × band rows in Program-32 daily table."
        )

    return (
        output.sort_values(
            [
                "threshold_c",
                "band_key",
                "date",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def validate_index_table(
    data: pd.DataFrame,
    expected_columns: tuple[str, ...],
    source: Path,
) -> pd.DataFrame:
    missing = [
        column
        for column in expected_columns
        if column not in data.columns
    ]

    if missing:
        raise ValueError(
            f"Program-21 contract violation in {source}: "
            f"missing columns {missing}"
        )

    output = data.copy()

    output[
        "date"
    ] = pd.to_datetime(
        output[
            "date"
        ],
        errors="raise",
    )

    if output[
        "date"
    ].duplicated().any():
        raise ValueError(
            f"Duplicate dates in Program-21 product: {source}"
        )

    return (
        output.sort_values(
            "date"
        )
        .reset_index(
            drop=True
        )
    )


def load_weekly_nino() -> pd.DataFrame:
    data = pd.read_csv(
        WEEKLY_NINO_FILE,
        low_memory=False,
    )

    data = validate_index_table(
        data,
        EXPECTED_NINO_COLUMNS,
        WEEKLY_NINO_FILE,
    )

    output = pd.DataFrame(
        {
            "date":
            data[
                "date"
            ]
        }
    )

    for region in NINO_REGION_KEYS:
        column = nino_column(
            region
        )

        output[
            column
        ] = pd.to_numeric(
            data[
                column
            ],
            errors="coerce",
        )

    return output


def load_monthly_soi() -> pd.DataFrame:
    data = pd.read_csv(
        MONTHLY_SOI_FILE,
        low_memory=False,
    )

    data = validate_index_table(
        data,
        EXPECTED_SOI_COLUMNS,
        MONTHLY_SOI_FILE,
    )

    return pd.DataFrame(
        {
            "date":
            data[
                "date"
            ],
            "soi":
            pd.to_numeric(
                data[
                    "soi"
                ],
                errors="coerce",
            ),
        }
    )


# =============================================================================
# MONTHLY COMMON DATA
# =============================================================================

def monthly_program32_power(
    daily: pd.DataFrame,
) -> pd.DataFrame:
    valid = daily.loc[
        daily[
            "strict_coi_valid"
        ]
        & np.isfinite(
            daily[
                "band_power_strict_coi"
            ].to_numpy(
                dtype=float
            )
        )
    ].copy()

    valid[
        "date"
    ] = (
        valid[
            "date"
        ]
        .dt.to_period(
            "M"
        )
        .dt.to_timestamp()
    )

    grouped = (
        valid.groupby(
            [
                "threshold_c",
                "band_key",
                "band_label",
                "minimum_period_days",
                "maximum_period_days",
                "date",
            ],
            as_index=False,
            observed=True,
        )
        .agg(
            band_power_mean=(
                "band_power_strict_coi",
                "mean",
            ),
            band_power_median=(
                "band_power_strict_coi",
                "median",
            ),
            band_power_std=(
                "band_power_strict_coi",
                "std",
            ),
            band_power_days=(
                "band_power_strict_coi",
                "count",
            ),
        )
    )

    grouped = grouped.loc[
        grouped[
            "band_power_days"
        ]
        >= MINIMUM_VALID_DAYS_PER_MONTH
    ].copy()

    grouped[
        "band_order"
    ] = grouped[
        "band_key"
    ].map(
        BAND_ORDER
    )

    grouped[
        "hypothesis_family"
    ] = grouped[
        "band_key"
    ].map(
        {
            row["band_key"]:
            row["hypothesis_family"]
            for row in EXPECTED_BANDS
        }
    )

    return (
        grouped.sort_values(
            [
                "threshold_c",
                "band_order",
                "date",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def monthly_nino(
    weekly: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        nino_column(
            region
        )
        for region in NINO_REGION_KEYS
    ]

    return (
        weekly.set_index(
            "date"
        )[columns]
        .resample(
            "MS"
        )
        .mean()
        .reset_index()
    )


def normalize_month_start(
    data: pd.DataFrame,
) -> pd.DataFrame:
    output = data.copy()

    output[
        "date"
    ] = (
        output[
            "date"
        ]
        .dt.to_period(
            "M"
        )
        .dt.to_timestamp()
    )

    return output


def build_monthly_common(
    monthly_power: pd.DataFrame,
    weekly_nino: pd.DataFrame,
    monthly_soi: pd.DataFrame,
) -> pd.DataFrame:
    nino = monthly_nino(
        weekly_nino
    )

    soi = normalize_month_start(
        monthly_soi
    )

    climate = (
        nino.merge(
            soi,
            on="date",
            how="inner",
            validate="one_to_one",
        )
        .sort_values(
            "date"
        )
        .reset_index(
            drop=True
        )
    )

    for region in NINO_REGION_KEYS:
        column = nino_column(
            region
        )

        climate[
            f"{column}_magnitude"
        ] = np.abs(
            climate[
                column
            ].to_numpy(
                dtype=float
            )
        )

    climate[
        "soi_magnitude"
    ] = np.abs(
        climate[
            "soi"
        ].to_numpy(
            dtype=float
        )
    )

    merged = monthly_power.merge(
        climate,
        on="date",
        how="inner",
        validate="many_to_one",
    )

    # Validate that every threshold x band retains a scientifically useful
    # monthly sample.
    counts = (
        merged.groupby(
            [
                "threshold_c",
                "band_key",
            ],
            observed=True,
        )[
            "date"
        ]
        .count()
    )

    insufficient = counts.loc[
        counts
        < MINIMUM_MONTHLY_RECORDS
    ]

    if not insufficient.empty:
        raise ValueError(
            "One or more threshold × band combinations have insufficient "
            "monthly common data:\n"
            f"{insufficient}"
        )

    return (
        merged.sort_values(
            [
                "threshold_c",
                "band_order",
                "date",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def predictor_definitions() -> list[
    dict[
        str,
        str,
    ]
]:
    definitions: list[
        dict[
            str,
            str,
        ]
    ] = []

    for region in NINO_REGION_KEYS:
        column = nino_column(
            region
        )

        label = (
            f"{NINO_LABELS[region]} "
            f"{NINO_ANALYSIS_FIELD}"
        )

        definitions.append(
            {
                "index_key": region,
                "index_label": label,
                "test_type": "signed_state",
                "predictor_column": column,
                "predictor_label": label,
            }
        )

        definitions.append(
            {
                "index_key": region,
                "index_label": label,
                "test_type": "enso_magnitude",
                "predictor_column":
                f"{column}_magnitude",
                "predictor_label":
                f"|{label}|",
            }
        )

    definitions.append(
        {
            "index_key": "soi",
            "index_label": SOI_LABEL,
            "test_type": "signed_state",
            "predictor_column": "soi",
            "predictor_label": SOI_LABEL,
        }
    )

    definitions.append(
        {
            "index_key": "soi",
            "index_label": SOI_LABEL,
            "test_type": "enso_magnitude",
            "predictor_column": "soi_magnitude",
            "predictor_label": "|SOI|",
        }
    )

    return definitions


# =============================================================================
# ZERO-LAG ANALYSIS
# =============================================================================

def zero_lag_analysis(
    monthly_common: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[
        dict[
            str,
            Any,
        ]
    ] = []

    grouping = (
        monthly_common[
            [
                "threshold_c",
                "band_key",
                "band_label",
                "minimum_period_days",
                "maximum_period_days",
                "band_order",
                "hypothesis_family",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "threshold_c",
                "band_order",
            ]
        )
    )

    for band_row in grouping.itertuples(
        index=False
    ):
        subset = monthly_common.loc[
            threshold_equal(
                monthly_common[
                    "threshold_c"
                ],
                float(
                    band_row.threshold_c
                ),
            )
            & monthly_common[
                "band_key"
            ].eq(
                band_row.band_key
            )
        ].copy()

        for definition in predictor_definitions():
            metrics = correlation_metrics(
                subset[
                    definition[
                        "predictor_column"
                    ]
                ].to_numpy(
                    dtype=float
                ),
                subset[
                    "band_power_mean"
                ].to_numpy(
                    dtype=float
                ),
            )

            rows.append(
                {
                    "threshold_c":
                    float(
                        band_row.threshold_c
                    ),
                    "band_key":
                    band_row.band_key,
                    "band_label":
                    band_row.band_label,
                    "band_order":
                    int(
                        band_row.band_order
                    ),
                    "hypothesis_family":
                    band_row.hypothesis_family,
                    "minimum_period_days":
                    float(
                        band_row.minimum_period_days
                    ),
                    "maximum_period_days":
                    float(
                        band_row.maximum_period_days
                    ),
                    **definition,
                    **metrics,
                }
            )

    output = pd.DataFrame(
        rows
    )

    output[
        "pearson_q_within_band"
    ] = np.nan
    output[
        "spearman_q_within_band"
    ] = np.nan
    output[
        "pearson_q_global_bands"
    ] = np.nan
    output[
        "spearman_q_global_bands"
    ] = np.nan

    # Program-28-like FDR: five indices within each threshold x band x test type.
    within_groups = output.groupby(
        [
            "threshold_c",
            "band_key",
            "test_type",
        ],
        sort=False,
    ).groups

    for _, indices in within_groups.items():
        idx = np.asarray(
            list(
                indices
            ),
            dtype=int,
        )

        output.loc[
            idx,
            "pearson_q_within_band",
        ] = benjamini_hochberg(
            output.loc[
                idx,
                "pearson_p_eff",
            ].to_numpy(
                dtype=float
            )
        )

        output.loc[
            idx,
            "spearman_q_within_band",
        ] = benjamini_hochberg(
            output.loc[
                idx,
                "spearman_p_eff",
            ].to_numpy(
                dtype=float
            )
        )

    # More conservative Program-33 family:
    # 7 bands x 5 ENSO indices for each threshold x test type.
    global_groups = output.groupby(
        [
            "threshold_c",
            "test_type",
        ],
        sort=False,
    ).groups

    for _, indices in global_groups.items():
        idx = np.asarray(
            list(
                indices
            ),
            dtype=int,
        )

        output.loc[
            idx,
            "pearson_q_global_bands",
        ] = benjamini_hochberg(
            output.loc[
                idx,
                "pearson_p_eff",
            ].to_numpy(
                dtype=float
            )
        )

        output.loc[
            idx,
            "spearman_q_global_bands",
        ] = benjamini_hochberg(
            output.loc[
                idx,
                "spearman_p_eff",
            ].to_numpy(
                dtype=float
            )
        )

    output[
        "pearson_significant_within_band"
    ] = (
        output[
            "pearson_q_within_band"
        ]
        <= FDR_ALPHA
    )

    output[
        "spearman_significant_within_band"
    ] = (
        output[
            "spearman_q_within_band"
        ]
        <= FDR_ALPHA
    )

    output[
        "pearson_significant_global_bands"
    ] = (
        output[
            "pearson_q_global_bands"
        ]
        <= FDR_ALPHA
    )

    output[
        "spearman_significant_global_bands"
    ] = (
        output[
            "spearman_q_global_bands"
        ]
        <= FDR_ALPHA
    )

    output[
        "index_order"
    ] = output[
        "index_key"
    ].map(
        INDEX_ORDER
    )

    return (
        output.sort_values(
            [
                "threshold_c",
                "band_order",
                "test_type",
                "index_order",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# CALENDAR-MONTH LAG ANALYSIS
# =============================================================================

def calendar_lag_pairs(
    monthly: pd.DataFrame,
    predictor_column: str,
    lag_months: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    pd.DatetimeIndex,
]:
    """
    Positive lag:
        ENSO predictor at time t leads PWP band power at t + lag.

    Calendar-month matching is used, so missing months cannot silently change
    the meaning of the lag.
    """

    predictor = monthly[
        [
            "date",
            predictor_column,
        ]
    ].copy()

    predictor[
        "response_date"
    ] = (
        predictor[
            "date"
        ]
        + pd.DateOffset(
            months=int(
                lag_months
            )
        )
    )

    response = (
        monthly[
            [
                "date",
                "band_power_mean",
            ]
        ]
        .rename(
            columns={
                "date":
                "response_date",
            }
        )
    )

    merged = predictor.merge(
        response,
        on="response_date",
        how="inner",
        validate="one_to_one",
    )

    x = pd.to_numeric(
        merged[
            predictor_column
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    y = pd.to_numeric(
        merged[
            "band_power_mean"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    valid = (
        np.isfinite(x)
        & np.isfinite(y)
    )

    return (
        x[valid],
        y[valid],
        pd.DatetimeIndex(
            merged.loc[
                valid,
                "response_date",
            ]
        ),
    )


def lag_curve(
    monthly: pd.DataFrame,
    definition: dict[
        str,
        str,
    ],
) -> pd.DataFrame:
    threshold_c = float(
        monthly[
            "threshold_c"
        ].iloc[0]
    )

    band_key = str(
        monthly[
            "band_key"
        ].iloc[0]
    )

    band_label = str(
        monthly[
            "band_label"
        ].iloc[0]
    )

    band_order = int(
        monthly[
            "band_order"
        ].iloc[0]
    )

    hypothesis_family = str(
        monthly[
            "hypothesis_family"
        ].iloc[0]
    )

    minimum_period_days = float(
        monthly[
            "minimum_period_days"
        ].iloc[0]
    )

    maximum_period_days = float(
        monthly[
            "maximum_period_days"
        ].iloc[0]
    )

    rows: list[
        dict[
            str,
            Any,
        ]
    ] = []

    for lag in range(
        -MONTHLY_MAX_LAG,
        MONTHLY_MAX_LAG + 1,
    ):
        x, y, paired_dates = calendar_lag_pairs(
            monthly=monthly,
            predictor_column=definition[
                "predictor_column"
            ],
            lag_months=lag,
        )

        if len(x) < MINIMUM_MONTHLY_RECORDS:
            continue

        metrics = correlation_metrics(
            x,
            y,
        )

        rows.append(
            {
                "threshold_c":
                threshold_c,
                "band_key":
                band_key,
                "band_label":
                band_label,
                "band_order":
                band_order,
                "hypothesis_family":
                hypothesis_family,
                "minimum_period_days":
                minimum_period_days,
                "maximum_period_days":
                maximum_period_days,
                "index_key":
                definition[
                    "index_key"
                ],
                "index_label":
                definition[
                    "index_label"
                ],
                "test_type":
                definition[
                    "test_type"
                ],
                "predictor_column":
                definition[
                    "predictor_column"
                ],
                "predictor_label":
                definition[
                    "predictor_label"
                ],
                "lag_months":
                int(
                    lag
                ),
                "n":
                int(
                    metrics[
                        "n"
                    ]
                ),
                "n_eff":
                metrics[
                    "n_eff"
                ],
                "pearson_r":
                metrics[
                    "pearson_r"
                ],
                "p_eff":
                metrics[
                    "pearson_p_eff"
                ],
                "paired_start_date":
                (
                    paired_dates.min()
                    if len(
                        paired_dates
                    )
                    else pd.NaT
                ),
                "paired_end_date":
                (
                    paired_dates.max()
                    if len(
                        paired_dates
                    )
                    else pd.NaT
                ),
            }
        )

    output = pd.DataFrame(
        rows
    )

    if output.empty:
        return output

    output[
        "q_fdr_within_lag_search"
    ] = benjamini_hochberg(
        output[
            "p_eff"
        ].to_numpy(
            dtype=float
        )
    )

    output[
        "fdr_significant"
    ] = (
        output[
            "q_fdr_within_lag_search"
        ]
        <= FDR_ALPHA
    )

    return output


def all_lag_curves(
    monthly_common: pd.DataFrame,
) -> pd.DataFrame:
    frames: list[
        pd.DataFrame
    ] = []

    groups = (
        monthly_common[
            [
                "threshold_c",
                "band_key",
                "band_order",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "threshold_c",
                "band_order",
            ]
        )
    )

    for group in groups.itertuples(
        index=False
    ):
        monthly = monthly_common.loc[
            threshold_equal(
                monthly_common[
                    "threshold_c"
                ],
                float(
                    group.threshold_c
                ),
            )
            & monthly_common[
                "band_key"
            ].eq(
                group.band_key
            )
        ].copy()

        for definition in predictor_definitions():
            curve = lag_curve(
                monthly,
                definition,
            )

            if not curve.empty:
                frames.append(
                    curve
                )

    if not frames:
        raise ValueError(
            "Program 33 generated no valid lag curves."
        )

    output = pd.concat(
        frames,
        ignore_index=True,
    )

    output[
        "index_order"
    ] = output[
        "index_key"
    ].map(
        INDEX_ORDER
    )

    return (
        output.sort_values(
            [
                "threshold_c",
                "band_order",
                "test_type",
                "index_order",
                "lag_months",
            ]
        )
        .reset_index(
            drop=True
        )
    )


def summarize_lag_curves(
    curves: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[
        dict[
            str,
            Any,
        ]
    ] = []

    grouping = [
        "threshold_c",
        "band_key",
        "band_label",
        "band_order",
        "hypothesis_family",
        "minimum_period_days",
        "maximum_period_days",
        "index_key",
        "index_label",
        "test_type",
        "predictor_label",
    ]

    for keys, subset in curves.groupby(
        grouping,
        sort=False,
    ):
        (
            threshold_c,
            band_key,
            band_label,
            band_order,
            hypothesis_family,
            minimum_period_days,
            maximum_period_days,
            index_key,
            index_label,
            test_type,
            predictor_label,
        ) = keys

        valid = subset.loc[
            np.isfinite(
                subset[
                    "pearson_r"
                ].to_numpy(
                    dtype=float
                )
            )
        ].copy()

        if valid.empty:
            continue

        strongest_absolute = valid.iloc[
            np.nanargmax(
                np.abs(
                    valid[
                        "pearson_r"
                    ].to_numpy(
                        dtype=float
                    )
                )
            )
        ]

        significant = valid.loc[
            valid[
                "fdr_significant"
            ]
        ].copy()

        if significant.empty:
            best_lag = np.nan
            best_r = np.nan
            best_q = np.nan
            best_p = np.nan
            best_n_eff = np.nan
        else:
            best = significant.iloc[
                np.nanargmax(
                    np.abs(
                        significant[
                            "pearson_r"
                        ].to_numpy(
                            dtype=float
                        )
                    )
                )
            ]

            best_lag = float(
                best[
                    "lag_months"
                ]
            )
            best_r = float(
                best[
                    "pearson_r"
                ]
            )
            best_q = float(
                best[
                    "q_fdr_within_lag_search"
                ]
            )
            best_p = float(
                best[
                    "p_eff"
                ]
            )
            best_n_eff = float(
                best[
                    "n_eff"
                ]
            )

        rows.append(
            {
                "threshold_c":
                float(
                    threshold_c
                ),
                "band_key":
                band_key,
                "band_label":
                band_label,
                "band_order":
                int(
                    band_order
                ),
                "hypothesis_family":
                hypothesis_family,
                "minimum_period_days":
                float(
                    minimum_period_days
                ),
                "maximum_period_days":
                float(
                    maximum_period_days
                ),
                "index_key":
                index_key,
                "index_label":
                index_label,
                "index_order":
                INDEX_ORDER[
                    index_key
                ],
                "test_type":
                test_type,
                "predictor_label":
                predictor_label,
                "maximum_absolute_lag_months":
                int(
                    strongest_absolute[
                        "lag_months"
                    ]
                ),
                "maximum_absolute_r":
                float(
                    strongest_absolute[
                        "pearson_r"
                    ]
                ),
                "maximum_absolute_p_eff":
                float(
                    strongest_absolute[
                        "p_eff"
                    ]
                ),
                "maximum_absolute_q_fdr":
                float(
                    strongest_absolute[
                        "q_fdr_within_lag_search"
                    ]
                ),
                "best_fdr_significant_lag_months":
                best_lag,
                "best_fdr_significant_r":
                best_r,
                "best_fdr_significant_p_eff":
                best_p,
                "best_fdr_significant_q":
                best_q,
                "best_fdr_significant_n_eff":
                best_n_eff,
                "significant_lag_count":
                int(
                    significant.shape[0]
                ),
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "threshold_c",
                "band_order",
                "test_type",
                "index_order",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# THRESHOLD ROBUSTNESS
# =============================================================================

def threshold_robustness(
    zero_lag: pd.DataFrame,
    lag_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[
        dict[
            str,
            Any,
        ]
    ] = []

    grouping = [
        "band_key",
        "band_label",
        "band_order",
        "hypothesis_family",
        "index_key",
        "index_label",
        "index_order",
        "test_type",
        "predictor_label",
    ]

    for keys, zero_subset in zero_lag.groupby(
        grouping,
        sort=False,
    ):
        (
            band_key,
            band_label,
            band_order,
            hypothesis_family,
            index_key,
            index_label,
            index_order,
            test_type,
            predictor_label,
        ) = keys

        zero_subset = zero_subset.sort_values(
            "threshold_c"
        )

        r_values = zero_subset[
            "pearson_r"
        ].to_numpy(
            dtype=float
        )

        finite_r = r_values[
            np.isfinite(
                r_values
            )
        ]

        same_zero_sign = bool(
            finite_r.size > 0
            and np.all(
                np.sign(
                    finite_r
                )
                == np.sign(
                    finite_r[0]
                )
            )
        )

        lag_subset = lag_summary.loc[
            lag_summary[
                "band_key"
            ].eq(
                band_key
            )
            & lag_summary[
                "index_key"
            ].eq(
                index_key
            )
            & lag_summary[
                "test_type"
            ].eq(
                test_type
            )
        ].copy()

        lag_valid = lag_subset.loc[
            np.isfinite(
                lag_subset[
                    "best_fdr_significant_r"
                ].to_numpy(
                    dtype=float
                )
            )
        ].copy()

        if lag_valid.empty:
            lag_count = 0
            mean_lag_r = np.nan
            mean_lag = np.nan
            lag_std = np.nan
            same_lag_r_sign = False
            lag_sign_consistent = False
        else:
            lag_r_values = lag_valid[
                "best_fdr_significant_r"
            ].to_numpy(
                dtype=float
            )

            lag_values = lag_valid[
                "best_fdr_significant_lag_months"
            ].to_numpy(
                dtype=float
            )

            nonzero_signs = np.sign(
                lag_values[
                    ~np.isclose(
                        lag_values,
                        0.0,
                    )
                ]
            )

            lag_count = int(
                len(
                    lag_valid
                )
            )

            mean_lag_r = float(
                np.mean(
                    lag_r_values
                )
            )

            mean_lag = float(
                np.mean(
                    lag_values
                )
            )

            lag_std = (
                float(
                    np.std(
                        lag_values,
                        ddof=1,
                    )
                )
                if len(
                    lag_values
                ) >= 2
                else 0.0
            )

            same_lag_r_sign = bool(
                np.all(
                    np.sign(
                        lag_r_values
                    )
                    == np.sign(
                        lag_r_values[0]
                    )
                )
            )

            lag_sign_consistent = bool(
                nonzero_signs.size <= 1
                or np.all(
                    nonzero_signs
                    == nonzero_signs[0]
                )
            )

        rows.append(
            {
                "band_key":
                band_key,
                "band_label":
                band_label,
                "band_order":
                int(
                    band_order
                ),
                "hypothesis_family":
                hypothesis_family,
                "index_key":
                index_key,
                "index_label":
                index_label,
                "index_order":
                int(
                    index_order
                ),
                "test_type":
                test_type,
                "predictor_label":
                predictor_label,
                "thresholds_available":
                int(
                    zero_subset[
                        "threshold_c"
                    ].nunique()
                ),
                "mean_zero_lag_pearson_r":
                float(
                    np.nanmean(
                        r_values
                    )
                ),
                "minimum_zero_lag_pearson_r":
                float(
                    np.nanmin(
                        r_values
                    )
                ),
                "maximum_zero_lag_pearson_r":
                float(
                    np.nanmax(
                        r_values
                    )
                ),
                "std_zero_lag_pearson_r":
                (
                    float(
                        np.nanstd(
                            r_values,
                            ddof=1,
                        )
                    )
                    if np.count_nonzero(
                        np.isfinite(
                            r_values
                        )
                    ) >= 2
                    else 0.0
                ),
                "same_zero_lag_correlation_sign":
                same_zero_sign,
                "within_band_significant_thresholds":
                int(
                    np.count_nonzero(
                        zero_subset[
                            "pearson_significant_within_band"
                        ].to_numpy(
                            dtype=bool
                        )
                    )
                ),
                "global_band_significant_thresholds":
                int(
                    np.count_nonzero(
                        zero_subset[
                            "pearson_significant_global_bands"
                        ].to_numpy(
                            dtype=bool
                        )
                    )
                ),
                "lag_significant_thresholds":
                lag_count,
                "mean_best_lag_r":
                mean_lag_r,
                "mean_best_lag_months":
                mean_lag,
                "lag_std_months":
                lag_std,
                "same_best_lag_correlation_sign":
                same_lag_r_sign,
                "lag_sign_consistent":
                lag_sign_consistent,
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "band_order",
                "test_type",
                "index_order",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# PRIMARY HYPOTHESIS SUMMARY
# =============================================================================

def build_hypothesis_summary(
    zero_lag: pd.DataFrame,
    lag_summary: pd.DataFrame,
    robustness: pd.DataFrame,
) -> pd.DataFrame:
    primary = zero_lag.loc[
        threshold_equal(
            zero_lag[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
    ].copy()

    rows: list[
        dict[
            str,
            Any,
        ]
    ] = []

    for band in EXPECTED_BANDS:
        for test_type in (
            "signed_state",
            "enso_magnitude",
        ):
            subset = primary.loc[
                primary[
                    "band_key"
                ].eq(
                    band[
                        "band_key"
                    ]
                )
                & primary[
                    "test_type"
                ].eq(
                    test_type
                )
            ].copy()

            if subset.empty:
                continue

            subset[
                "absolute_r"
            ] = np.abs(
                subset[
                    "pearson_r"
                ].to_numpy(
                    dtype=float
                )
            )

            strongest = (
                subset.sort_values(
                    "absolute_r",
                    ascending=False,
                )
                .iloc[0]
            )

            any_global = bool(
                subset[
                    "pearson_significant_global_bands"
                ].any()
            )

            any_within = bool(
                subset[
                    "pearson_significant_within_band"
                ].any()
            )

            if any_global:
                zero_status = (
                    "SUPPORTED_ZERO_LAG_GLOBAL_FDR"
                )
            elif any_within:
                zero_status = (
                    "SUPPORTED_WITHIN_BAND_ONLY"
                )
            else:
                zero_status = (
                    "NOT_SUPPORTED_AT_ZERO_LAG"
                )

            lag_subset = lag_summary.loc[
                threshold_equal(
                    lag_summary[
                        "threshold_c"
                    ],
                    PRIMARY_THRESHOLD_C,
                )
                & lag_summary[
                    "band_key"
                ].eq(
                    band[
                        "band_key"
                    ]
                )
                & lag_summary[
                    "test_type"
                ].eq(
                    test_type
                )
                & np.isfinite(
                    lag_summary[
                        "best_fdr_significant_r"
                    ].to_numpy(
                        dtype=float
                    )
                )
            ].copy()

            if lag_subset.empty:
                lag_status = (
                    "NO_FDR_SIGNIFICANT_LAG"
                )
                strongest_lag_predictor = ""
                strongest_lag = np.nan
                strongest_lag_r = np.nan
                strongest_lag_q = np.nan
            else:
                lag_subset[
                    "absolute_best_lag_r"
                ] = np.abs(
                    lag_subset[
                        "best_fdr_significant_r"
                    ].to_numpy(
                        dtype=float
                    )
                )

                lag_best = (
                    lag_subset.sort_values(
                        "absolute_best_lag_r",
                        ascending=False,
                    )
                    .iloc[0]
                )

                lag_status = (
                    "EXPLORATORY_LAGGED_ASSOCIATION"
                )
                strongest_lag_predictor = str(
                    lag_best[
                        "predictor_label"
                    ]
                )
                strongest_lag = float(
                    lag_best[
                        "best_fdr_significant_lag_months"
                    ]
                )
                strongest_lag_r = float(
                    lag_best[
                        "best_fdr_significant_r"
                    ]
                )
                strongest_lag_q = float(
                    lag_best[
                        "best_fdr_significant_q"
                    ]
                )

            robust = robustness.loc[
                robustness[
                    "band_key"
                ].eq(
                    band[
                        "band_key"
                    ]
                )
                & robustness[
                    "index_key"
                ].eq(
                    strongest[
                        "index_key"
                    ]
                )
                & robustness[
                    "test_type"
                ].eq(
                    test_type
                )
            ]

            robust_row = (
                robust.iloc[0]
                if not robust.empty
                else None
            )

            rows.append(
                {
                    "band_key":
                    band[
                        "band_key"
                    ],
                    "band_label":
                    band[
                        "band_label"
                    ],
                    "band_order":
                    BAND_ORDER[
                        band[
                            "band_key"
                        ]
                    ],
                    "hypothesis_family":
                    band[
                        "hypothesis_family"
                    ],
                    "test_type":
                    test_type,
                    "primary_threshold_c":
                    PRIMARY_THRESHOLD_C,
                    "strongest_zero_lag_index_key":
                    strongest[
                        "index_key"
                    ],
                    "strongest_zero_lag_predictor":
                    strongest[
                        "predictor_label"
                    ],
                    "zero_lag_pearson_r":
                    float(
                        strongest[
                            "pearson_r"
                        ]
                    ),
                    "zero_lag_pearson_q_within_band":
                    float(
                        strongest[
                            "pearson_q_within_band"
                        ]
                    ),
                    "zero_lag_pearson_q_global_bands":
                    float(
                        strongest[
                            "pearson_q_global_bands"
                        ]
                    ),
                    "zero_lag_spearman_rho":
                    float(
                        strongest[
                            "spearman_rho"
                        ]
                    ),
                    "zero_lag_spearman_q_global_bands":
                    float(
                        strongest[
                            "spearman_q_global_bands"
                        ]
                    ),
                    "zero_lag_scientific_status":
                    zero_status,
                    "lag_scientific_status":
                    lag_status,
                    "strongest_lag_predictor":
                    strongest_lag_predictor,
                    "strongest_fdr_lag_months":
                    strongest_lag,
                    "strongest_fdr_lag_r":
                    strongest_lag_r,
                    "strongest_fdr_lag_q":
                    strongest_lag_q,
                    "strongest_zero_lag_global_significant_thresholds":
                    (
                        int(
                            robust_row[
                                "global_band_significant_thresholds"
                            ]
                        )
                        if robust_row is not None
                        else 0
                    ),
                    "strongest_zero_lag_same_sign_across_thresholds":
                    (
                        bool(
                            robust_row[
                                "same_zero_lag_correlation_sign"
                            ]
                        )
                        if robust_row is not None
                        else False
                    ),
                }
            )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "band_order",
                "test_type",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# OUTPUT TABLES
# =============================================================================

def export_tables(
    monthly_common: pd.DataFrame,
    zero_lag: pd.DataFrame,
    lag_curves: pd.DataFrame,
    lag_summary: pd.DataFrame,
    robustness: pd.DataFrame,
    hypothesis: pd.DataFrame,
) -> tuple[
    Path,
    ...,
]:
    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    monthly_out = monthly_common.copy()

    monthly_out[
        "date"
    ] = pd.to_datetime(
        monthly_out[
            "date"
        ]
    ).dt.strftime(
        "%Y-%m-%d"
    )

    lag_curves_out = lag_curves.copy()

    for column in (
        "paired_start_date",
        "paired_end_date",
    ):
        lag_curves_out[
            column
        ] = pd.to_datetime(
            lag_curves_out[
                column
            ]
        ).dt.strftime(
            "%Y-%m-%d"
        )

    monthly_out.to_csv(
        MONTHLY_DATA_FILE,
        index=False,
        float_format="%.10g",
    )

    zero_lag.to_csv(
        ZERO_LAG_FILE,
        index=False,
        float_format="%.10g",
    )

    lag_curves_out.to_csv(
        LAG_CURVES_FILE,
        index=False,
        float_format="%.10g",
    )

    lag_summary.to_csv(
        LAG_SUMMARY_FILE,
        index=False,
        float_format="%.10g",
    )

    robustness.to_csv(
        ROBUSTNESS_FILE,
        index=False,
        float_format="%.10g",
    )

    hypothesis.to_csv(
        HYPOTHESIS_FILE,
        index=False,
        float_format="%.10g",
    )

    return (
        MONTHLY_DATA_FILE,
        ZERO_LAG_FILE,
        LAG_CURVES_FILE,
        LAG_SUMMARY_FILE,
        ROBUSTNESS_FILE,
        HYPOTHESIS_FILE,
    )


# =============================================================================
# FIGURE 1 — ZERO-LAG HEATMAPS
# =============================================================================

def display_index_label(
    index_key: str,
) -> str:
    if index_key == "soi":
        return "SOI"

    return (
        NINO_LABELS[
            index_key
        ]
        .replace(
            "Niño ",
            "Niño "
        )
    )


def plot_zero_lag_heatmaps(
    zero_lag: pd.DataFrame,
) -> tuple[
    Path,
    Path,
]:
    primary = zero_lag.loc[
        threshold_equal(
            zero_lag[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
    ].copy()

    band_keys = [
        row[
            "band_key"
        ]
        for row in EXPECTED_BANDS
    ]

    band_labels = [
        row[
            "band_label"
        ]
        for row in EXPECTED_BANDS
    ]

    index_keys = [
        "nino12",
        "nino3",
        "nino34",
        "nino4",
        "soi",
    ]

    index_labels = [
        display_index_label(
            key
        )
        for key in index_keys
    ]

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(
            14.5,
            7.8,
        ),
        constrained_layout=True,
    )

    image = None

    for panel_index, (
        axis,
        test_type,
        title,
    ) in enumerate(
        (
            (
                axes[0],
                "signed_state",
                "Signed ENSO state",
            ),
            (
                axes[1],
                "enso_magnitude",
                "ENSO magnitude",
            ),
        )
    ):
        matrix = np.full(
            (
                len(
                    band_keys
                ),
                len(
                    index_keys
                ),
            ),
            np.nan,
            dtype=float,
        )

        significant = np.zeros(
            matrix.shape,
            dtype=bool,
        )

        for i, band_key in enumerate(
            band_keys
        ):
            for j, index_key in enumerate(
                index_keys
            ):
                row = primary.loc[
                    primary[
                        "band_key"
                    ].eq(
                        band_key
                    )
                    & primary[
                        "index_key"
                    ].eq(
                        index_key
                    )
                    & primary[
                        "test_type"
                    ].eq(
                        test_type
                    )
                ]

                if row.empty:
                    continue

                record = row.iloc[0]

                matrix[
                    i,
                    j,
                ] = float(
                    record[
                        "pearson_r"
                    ]
                )

                significant[
                    i,
                    j,
                ] = bool(
                    record[
                        "pearson_significant_global_bands"
                    ]
                )

        image = axis.imshow(
            matrix,
            aspect="auto",
            vmin=-1.0,
            vmax=1.0,
        )

        axis.set_xticks(
            np.arange(
                len(
                    index_labels
                )
            )
        )

        axis.set_xticklabels(
            index_labels,
            rotation=35,
            ha="right",
        )

        axis.set_yticks(
            np.arange(
                len(
                    band_labels
                )
            )
        )

        if panel_index == 0:
            axis.set_yticklabels(
                band_labels
            )
        else:
            axis.set_yticklabels(
                []
            )

        axis.set_title(
            title
        )

        for i in range(
            matrix.shape[0]
        ):
            for j in range(
                matrix.shape[1]
            ):
                if not np.isfinite(
                    matrix[
                        i,
                        j,
                    ]
                ):
                    continue

                text = (
                    f"{matrix[i, j]:+.2f}"
                    + (
                        "\n*"
                        if significant[
                            i,
                            j
                        ]
                        else ""
                    )
                )

                axis.text(
                    j,
                    i,
                    text,
                    ha="center",
                    va="center",
                    fontsize=8,
                )

        add_panel_label(
            axis,
            (
                "A)"
                if panel_index == 0
                else "B)"
            ),
        )

    if image is not None:
        cbar = figure.colorbar(
            image,
            ax=axes,
            orientation="vertical",
            fraction=0.03,
            pad=0.025,
        )

        cbar.set_label(
            "Pearson correlation, r"
        )

    figure.suptitle(
        (
            "ENSO association with time-varying PWP longitude wavelet variance "
            "— 28 °C"
        ),
        fontsize=14,
    )

    figure.text(
        0.5,
        0.01,
        "* survives Program-33 global-band FDR (7 bands × 5 indices)",
        ha="center",
        va="bottom",
        fontsize=9,
    )

    return save_figure(
        figure,
        "pwp_enso_wavelet_band_zero_lag_heatmaps_28C",
    )


# =============================================================================
# FIGURE 2 — KEY BAND × NIÑO 3.4 LAG CURVES
# =============================================================================

def plot_key_band_nino34_lag_curves(
    lag_curves: pd.DataFrame,
) -> tuple[
    Path,
    Path,
]:
    key_bands = (
        "annual_broad",
        "annual_strict",
        "interannual_3",
        "enso_broad",
    )

    primary = lag_curves.loc[
        threshold_equal(
            lag_curves[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
        & lag_curves[
            "index_key"
        ].eq(
            "nino34"
        )
        & lag_curves[
            "band_key"
        ].isin(
            key_bands
        )
    ].copy()

    figure, axes = plt.subplots(
        nrows=len(
            key_bands
        ),
        ncols=2,
        figsize=(
            13.5,
            13.0,
        ),
        sharex=True,
        sharey=False,
        constrained_layout=True,
    )

    test_specs = (
        (
            "signed_state",
            "Signed Niño 3.4 state",
        ),
        (
            "enso_magnitude",
            "|Niño 3.4|",
        ),
    )

    for row_index, band_key in enumerate(
        key_bands
    ):
        band_label = next(
            row[
                "band_label"
            ]
            for row in EXPECTED_BANDS
            if row[
                "band_key"
            ] == band_key
        )

        for column_index, (
            test_type,
            test_label,
        ) in enumerate(
            test_specs
        ):
            axis = axes[
                row_index,
                column_index,
            ]

            subset = primary.loc[
                primary[
                    "band_key"
                ].eq(
                    band_key
                )
                & primary[
                    "test_type"
                ].eq(
                    test_type
                )
            ].sort_values(
                "lag_months"
            )

            axis.plot(
                subset[
                    "lag_months"
                ],
                subset[
                    "pearson_r"
                ],
                linewidth=1.1,
                label=test_label,
            )

            significant = subset.loc[
                subset[
                    "fdr_significant"
                ]
            ]

            axis.scatter(
                significant[
                    "lag_months"
                ],
                significant[
                    "pearson_r"
                ],
                s=16,
                label="FDR-significant lag",
            )

            axis.axhline(
                0.0,
                linewidth=0.7,
                linestyle="--",
            )

            axis.axvline(
                0.0,
                linewidth=0.7,
                linestyle=":",
            )

            axis.grid(
                linestyle=GRID_LINESTYLE,
                linewidth=GRID_LINEWIDTH,
                alpha=GRID_ALPHA,
            )

            if row_index == 0:
                axis.set_title(
                    test_label
                )

            if column_index == 0:
                axis.set_ylabel(
                    f"{band_label}\nPearson r"
                )

            if row_index == (
                len(
                    key_bands
                )
                - 1
            ):
                axis.set_xlabel(
                    "Lag (months; positive = ENSO leads PWP band power)"
                )

            panel_number = (
                row_index
                * 2
                + column_index
            )

            add_panel_label(
                axis,
                f"{chr(65 + panel_number)})",
            )

            if (
                row_index == 0
                and column_index == 1
            ):
                axis.legend(
                    frameon=False,
                    loc="best",
                    fontsize=8,
                )

    figure.suptitle(
        (
            "Niño 3.4 versus PWP longitude wavelet-band power — "
            "calendar-month lag structure — 28 °C"
        ),
        fontsize=14,
    )

    return save_figure(
        figure,
        "pwp_enso_wavelet_key_band_nino34_lag_curves_28C",
    )


# =============================================================================
# FIGURE 3 — ANNUAL BANDS AND NIÑO 3.4
# =============================================================================

def plot_annual_band_time_series(
    monthly_common: pd.DataFrame,
) -> tuple[
    Path,
    Path,
]:
    primary = monthly_common.loc[
        threshold_equal(
            monthly_common[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
        & monthly_common[
            "band_key"
        ].isin(
            (
                "annual_broad",
                "annual_strict",
            )
        )
    ].copy()

    broad = primary.loc[
        primary[
            "band_key"
        ].eq(
            "annual_broad"
        )
    ].sort_values(
        "date"
    )

    strict = primary.loc[
        primary[
            "band_key"
        ].eq(
            "annual_strict"
        )
    ].sort_values(
        "date"
    )

    nino34_column = nino_column(
        "nino34"
    )

    # Use the broad-band rows as the climate-index time base.
    climate = broad[
        [
            "date",
            nino34_column,
            f"{nino34_column}_magnitude",
        ]
    ].copy()

    figure, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(
            14.5,
            10.5,
        ),
        sharex=True,
        constrained_layout=True,
    )

    axes[0].plot(
        broad[
            "date"
        ],
        broad[
            "band_power_mean"
        ],
        linewidth=0.9,
        label="Annual Broad (280–448 d)",
    )

    axes[0].set_ylabel(
        "Annual Broad\npower"
    )
    axes[0].legend(
        frameon=False,
        loc="upper right",
    )
    add_panel_label(
        axes[0],
        "A)",
    )

    axes[1].plot(
        strict[
            "date"
        ],
        strict[
            "band_power_mean"
        ],
        linewidth=0.9,
        label="Annual Strict (336–378 d)",
    )

    axes[1].set_ylabel(
        "Annual Strict\npower"
    )
    axes[1].legend(
        frameon=False,
        loc="upper right",
    )
    add_panel_label(
        axes[1],
        "B)",
    )

    axes[2].plot(
        climate[
            "date"
        ],
        climate[
            nino34_column
        ],
        linewidth=0.9,
        label="Niño 3.4 SSTA",
    )

    axes[2].axhline(
        0.0,
        linewidth=0.7,
        linestyle="--",
    )

    axes[2].set_ylabel(
        "Niño 3.4\nSSTA (°C)"
    )
    axes[2].legend(
        frameon=False,
        loc="upper right",
    )
    add_panel_label(
        axes[2],
        "C)",
    )

    axes[3].plot(
        climate[
            "date"
        ],
        climate[
            f"{nino34_column}_magnitude"
        ],
        linewidth=0.9,
        label="|Niño 3.4 SSTA|",
    )

    axes[3].set_ylabel(
        "|Niño 3.4|\n(°C)"
    )
    axes[3].legend(
        frameon=False,
        loc="upper right",
    )
    add_panel_label(
        axes[3],
        "D)",
    )

    minimum_date = max(
        broad[
            "date"
        ].min(),
        strict[
            "date"
        ].min(),
        climate[
            "date"
        ].min(),
    )

    maximum_date = min(
        broad[
            "date"
        ].max(),
        strict[
            "date"
        ].max(),
        climate[
            "date"
        ].max(),
    )

    for index, axis in enumerate(
        axes
    ):
        configure_year_axis(
            axis,
            minimum_date,
            maximum_date,
            show_labels=(
                index
                == len(
                    axes
                )
                - 1
            ),
        )

    figure.suptitle(
        (
            "PWP annual-cycle wavelet variance and Niño 3.4 — "
            "28 °C primary definition"
        ),
        fontsize=14,
    )

    return save_figure(
        figure,
        "pwp_enso_wavelet_annual_band_time_series_28C",
    )


# =============================================================================
# REPORTING
# =============================================================================

def serializable_value(
    value: Any,
) -> Any:
    if isinstance(
        value,
        (
            np.integer,
        ),
    ):
        return int(
            value
        )

    if isinstance(
        value,
        (
            np.floating,
        ),
    ):
        if np.isfinite(
            value
        ):
            return float(
                value
            )
        return None

    if isinstance(
        value,
        (
            np.bool_,
        ),
    ):
        return bool(
            value
        )

    if isinstance(
        value,
        pd.Timestamp,
    ):
        return value.isoformat()

    if pd.isna(
        value
    ):
        return None

    return value


def dataframe_records_json(
    data: pd.DataFrame,
) -> list[
    dict[
        str,
        Any,
    ]
]:
    records: list[
        dict[
            str,
            Any,
        ]
    ] = []

    for record in data.to_dict(
        orient="records"
    ):
        records.append(
            {
                key:
                serializable_value(
                    value
                )
                for key, value in record.items()
            }
        )

    return records


def write_reports(
    band_definitions: pd.DataFrame,
    monthly_common: pd.DataFrame,
    zero_lag: pd.DataFrame,
    lag_summary: pd.DataFrame,
    robustness: pd.DataFrame,
    hypothesis: pd.DataFrame,
    created_files: list[
        Path
    ],
) -> tuple[
    Path,
    Path,
]:
    generated = datetime.now(
        timezone.utc
    ).isoformat()

    lines: list[
        str
    ] = [
        PROGRAM_NAME,
        "=" * 78,
        "",
        "1. PROGRAM",
        "-" * 78,
        f"Version                      : {PROGRAM_VERSION}",
        f"Generated UTC                : {generated}",
        f"Project root                 : {PROJECT_DIR}",
        "",
        "2. SCIENTIFIC PURPOSE",
        "-" * 78,
        (
            "Test ENSO association with FINAL Program-32 time-varying "
            "wavelet-band variance of PWP centroid longitude."
        ),
        (
            "Primary PWP definition = 28.0 °C; "
            "28.5 and 29.0 °C = sensitivity tests."
        ),
        "",
        "3. TEMPORAL / STATISTICAL POLICY",
        "-" * 78,
        "PWP band power              : daily strict-COI -> monthly mean",
        "Niño SSTA                   : weekly -> monthly mean",
        "SOI                         : monthly unchanged",
        "Climate-index upsampling    : NO",
        (
            f"Minimum valid PWP days/month: "
            f"{MINIMUM_VALID_DAYS_PER_MONTH}"
        ),
        (
            f"Lag window                  : "
            f"±{MONTHLY_MAX_LAG} calendar months"
        ),
        "Positive lag                : ENSO leads PWP band power",
        "Serial correlation          : AR(1) effective sample size",
        "Zero-lag FDR                : within-band + global-band families",
        "Lag FDR                     : within each 49-lag relationship",
        "",
        "4. FINAL PROGRAM-32 BANDS",
        "-" * 78,
    ]

    for row in band_definitions.sort_values(
        "band_order"
    ).itertuples(
        index=False
    ):
        lines.append(
            (
                f"{row.band_label:<20s}: "
                f"{row.minimum_period_days:.0f}–"
                f"{row.maximum_period_days:.0f} days"
            )
        )

    lines.extend(
        [
            "",
            "5. PRIMARY 28 °C ZERO-LAG HYPOTHESIS SUMMARY",
            "-" * 78,
        ]
    )

    for row in hypothesis.itertuples(
        index=False
    ):
        lines.append(
            (
                f"{row.band_label} [{row.test_type}] | "
                f"strongest={row.strongest_zero_lag_predictor} | "
                f"r={row.zero_lag_pearson_r:+.3f} | "
                f"q_within={row.zero_lag_pearson_q_within_band:.3g} | "
                f"q_global={row.zero_lag_pearson_q_global_bands:.3g} | "
                f"{row.zero_lag_scientific_status}"
            )
        )

    lines.extend(
        [
            "",
            "6. PRIMARY 28 °C EXPLORATORY BEST FDR-SIGNIFICANT LAGS",
            "-" * 78,
        ]
    )

    primary_lags = lag_summary.loc[
        threshold_equal(
            lag_summary[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
    ].copy()

    for row in primary_lags.itertuples(
        index=False
    ):
        if np.isfinite(
            row.best_fdr_significant_lag_months
        ):
            lines.append(
                (
                    f"{row.band_label} | "
                    f"{row.predictor_label} [{row.test_type}] | "
                    f"lag={row.best_fdr_significant_lag_months:+.0f} months | "
                    f"r={row.best_fdr_significant_r:+.3f} | "
                    f"q={row.best_fdr_significant_q:.3g}"
                )
            )

    lines.extend(
        [
            "",
            "7. THRESHOLD ROBUSTNESS",
            "-" * 78,
        ]
    )

    for row in robustness.itertuples(
        index=False
    ):
        lines.append(
            (
                f"{row.band_label} | "
                f"{row.predictor_label} [{row.test_type}] | "
                f"mean zero-lag r={row.mean_zero_lag_pearson_r:+.3f} | "
                f"same sign={row.same_zero_lag_correlation_sign} | "
                f"global-significant thresholds="
                f"{row.global_band_significant_thresholds}/"
                f"{row.thresholds_available} | "
                f"lag-significant thresholds="
                f"{row.lag_significant_thresholds}/"
                f"{row.thresholds_available}"
            )
        )

    lines.extend(
        [
            "",
            "8. INTERPRETIVE SAFEGUARDS",
            "-" * 78,
            (
                "1. Program 33 tests wavelet-band-power modulation, not the "
                "already established direct PWP-longitude/ENSO displacement "
                "association."
            ),
            (
                "2. Zero-lag signed-state and ENSO-magnitude hypotheses are "
                "separate."
            ),
            (
                "3. Global-band FDR is the conservative primary zero-lag "
                "multiple-testing control."
            ),
            (
                "4. Lagged findings are exploratory even after FDR across the "
                "49-lag search because the strongest lag is selected from the "
                "lag family."
            ),
            (
                "5. Positive lag means ENSO leads PWP band power; lag does not "
                "prove causal direction."
            ),
            (
                "6. Program-32 strict COI filtering is retained. Edge-affected "
                "band power does not enter the monthly response."
            ),
            (
                "7. Overlapping bands are not independent partitions of total "
                "variance."
            ),
            (
                "8. Non-significant hypotheses remain valid scientific results."
            ),
            "",
            "9. DATA INVENTORY",
            "-" * 78,
            (
                f"Monthly common rows total  : "
                f"{len(monthly_common):,}"
            ),
            (
                f"Thresholds                  : "
                f"{monthly_common['threshold_c'].nunique()}"
            ),
            (
                f"Bands                       : "
                f"{monthly_common['band_key'].nunique()}"
            ),
            "",
            "10. FILES CREATED",
            "-" * 78,
        ]
    )

    for path in created_files:
        lines.append(
            str(
                path
            )
        )

    lines.extend(
        [
            "",
            "11. FINAL STATUS",
            "-" * 78,
            "Machine status               : PASS",
            "Scientific status            : READY FOR AUDIT-06 INTERPRETATION",
            "",
            "=" * 78,
        ]
    )

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_TXT.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )

    payload = {
        "program":
        PROGRAM_NAME,
        "version":
        PROGRAM_VERSION,
        "generated_utc":
        generated,
        "project_root":
        str(
            PROJECT_DIR
        ),
        "primary_threshold_c":
        PRIMARY_THRESHOLD_C,
        "thresholds_c":
        list(
            THRESHOLDS_C
        ),
        "monthly_max_lag":
        MONTHLY_MAX_LAG,
        "minimum_valid_days_per_month":
        MINIMUM_VALID_DAYS_PER_MONTH,
        "fdr_alpha":
        FDR_ALPHA,
        "band_definitions":
        dataframe_records_json(
            band_definitions
        ),
        "hypothesis_summary":
        dataframe_records_json(
            hypothesis
        ),
        "threshold_robustness":
        dataframe_records_json(
            robustness
        ),
        "interpretive_safeguards":
        [
            (
                "Association and lag do not establish causal direction."
            ),
            (
                "Signed ENSO state and ENSO magnitude are separate tests."
            ),
            (
                "Global-band FDR is the primary zero-lag family control."
            ),
            (
                "Lagged relationships are exploratory after within-relation "
                "FDR."
            ),
            (
                "Program-32 strict-COI filtering is retained."
            ),
            (
                "Overlapping wavelet bands are diagnostic and non-independent."
            ),
        ],
        "files_created":
        [
            str(
                path
            )
            for path in created_files
        ],
    }

    REPORT_JSON.write_text(
        json.dumps(
            payload,
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    return (
        REPORT_TXT,
        REPORT_JSON,
    )


# =============================================================================
# TERMINAL SUMMARY
# =============================================================================

def print_primary_summary(
    hypothesis: pd.DataFrame,
    lag_summary: pd.DataFrame,
) -> None:
    section(
        "## PRIMARY 28 °C — ZERO-LAG HYPOTHESIS SUMMARY"
    )

    for row in hypothesis.itertuples(
        index=False
    ):
        print(
            f"{row.band_label:<20s} | "
            f"{row.test_type:<14s} | "
            f"{row.strongest_zero_lag_predictor:<25s} | "
            f"r={row.zero_lag_pearson_r:+.3f} | "
            f"q_global={row.zero_lag_pearson_q_global_bands:.3g} | "
            f"{row.zero_lag_scientific_status}"
        )

    section(
        "## PRIMARY 28 °C — STRONGEST EXPLORATORY LAGGED ASSOCIATIONS"
    )

    primary_lags = lag_summary.loc[
        threshold_equal(
            lag_summary[
                "threshold_c"
            ],
            PRIMARY_THRESHOLD_C,
        )
        & np.isfinite(
            lag_summary[
                "best_fdr_significant_r"
            ].to_numpy(
                dtype=float
            )
        )
    ].copy()

    if primary_lags.empty:
        print(
            "No FDR-significant lagged relationships."
        )
        return

    for band_key in BAND_ORDER:
        subset = primary_lags.loc[
            primary_lags[
                "band_key"
            ].eq(
                band_key
            )
        ].copy()

        if subset.empty:
            continue

        subset[
            "abs_r"
        ] = np.abs(
            subset[
                "best_fdr_significant_r"
            ].to_numpy(
                dtype=float
            )
        )

        best = (
            subset.sort_values(
                "abs_r",
                ascending=False,
            )
            .iloc[0]
        )

        print(
            f"{best['band_label']:<20s} | "
            f"{best['predictor_label']:<25s} | "
            f"{best['test_type']:<14s} | "
            f"lag={best['best_fdr_significant_lag_months']:+.0f} months | "
            f"r={best['best_fdr_significant_r']:+.3f} | "
            f"q={best['best_fdr_significant_q']:.3g}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    rule()
    print(
        PROGRAM_NAME
    )
    rule()

    section(
        "## CONFIGURATION"
    )

    item(
        "Program version",
        PROGRAM_VERSION,
    )
    item(
        "Project root",
        PROJECT_ROOT,
    )
    item(
        "Primary PWP threshold",
        f"{PRIMARY_THRESHOLD_C:.1f} °C",
    )
    item(
        "Thresholds",
        ", ".join(
            f"{value:.1f}"
            for value in THRESHOLDS_C
        )
        + " °C",
    )
    item(
        "Common resolution",
        "monthly",
    )
    item(
        "Climate-index upsampling",
        "PROHIBITED",
    )
    item(
        "Program-32 CWT recomputation",
        "NO",
    )
    item(
        "Program-32 strict COI",
        "RETAINED",
    )
    item(
        "Minimum valid PWP days/month",
        MINIMUM_VALID_DAYS_PER_MONTH,
    )
    item(
        "Lag window",
        f"±{MONTHLY_MAX_LAG} months",
    )
    item(
        "Positive lag",
        "ENSO leads PWP wavelet-band power",
    )
    item(
        "FDR alpha",
        FDR_ALPHA,
    )
    item(
        "Zero-lag FDR families",
        "within-band + global-band",
    )
    item(
        "Lag FDR family",
        "49 lags within each relationship",
    )

    validate_configuration()

    section(
        "## LOADING AUTHORITATIVE INPUTS"
    )

    band_definitions = (
        load_and_validate_band_definitions()
    )

    daily_power = (
        load_program32_daily_power()
    )

    weekly_nino = (
        load_weekly_nino()
    )

    monthly_soi = (
        load_monthly_soi()
    )

    item(
        "Program-32 daily file",
        PROGRAM32_DAILY_FILE,
    )
    item(
        "Program-32 band file",
        PROGRAM32_BAND_FILE,
    )
    item(
        "Program-21 weekly Niño file",
        WEEKLY_NINO_FILE,
    )
    item(
        "Program-21 monthly SOI file",
        MONTHLY_SOI_FILE,
    )
    item(
        "Program-32 daily rows",
        f"{len(daily_power):,}",
    )
    item(
        "Program-32 bands",
        len(
            band_definitions
        ),
    )

    section(
        "## FINAL PROGRAM-32 BANDS"
    )

    for row in band_definitions.sort_values(
        "band_order"
    ).itertuples(
        index=False
    ):
        print(
            f"{row.band_label:<20s} : "
            f"{row.minimum_period_days:.0f}–"
            f"{row.maximum_period_days:.0f} days"
        )

    section(
        "## BUILDING MONTHLY COMMON DATA"
    )

    monthly_power = (
        monthly_program32_power(
            daily_power
        )
    )

    monthly_common = (
        build_monthly_common(
            monthly_power=monthly_power,
            weekly_nino=weekly_nino,
            monthly_soi=monthly_soi,
        )
    )

    for threshold_c in THRESHOLDS_C:
        print()
        print(
            f"{threshold_c:.1f} °C"
        )

        subset = monthly_common.loc[
            threshold_equal(
                monthly_common[
                    "threshold_c"
                ],
                threshold_c,
            )
        ]

        for band in EXPECTED_BANDS:
            band_data = subset.loc[
                subset[
                    "band_key"
                ].eq(
                    band[
                        "band_key"
                    ]
                )
            ]

            print(
                f"  {band['band_label']:<20s} | "
                f"N={len(band_data):4d} | "
                f"{band_data['date'].min():%Y-%m} to "
                f"{band_data['date'].max():%Y-%m}"
            )

    section(
        "## ZERO-LAG ASSOCIATION"
    )

    zero_lag = zero_lag_analysis(
        monthly_common
    )

    item(
        "Zero-lag relationships",
        f"{len(zero_lag):,}",
    )

    section(
        "## CALENDAR-MONTH LAG ANALYSIS"
    )

    lag_curves = all_lag_curves(
        monthly_common
    )

    lag_summary = summarize_lag_curves(
        lag_curves
    )

    item(
        "Lag-curve rows",
        f"{len(lag_curves):,}",
    )
    item(
        "Lag relationships summarized",
        f"{len(lag_summary):,}",
    )

    section(
        "## THRESHOLD ROBUSTNESS"
    )

    robustness = threshold_robustness(
        zero_lag=zero_lag,
        lag_summary=lag_summary,
    )

    item(
        "Robustness relationships",
        f"{len(robustness):,}",
    )

    section(
        "## PRIMARY HYPOTHESIS DECISIONS"
    )

    hypothesis = build_hypothesis_summary(
        zero_lag=zero_lag,
        lag_summary=lag_summary,
        robustness=robustness,
    )

    item(
        "Band-level hypothesis rows",
        len(
            hypothesis
        ),
    )

    section(
        "## EXPORTING TABLES"
    )

    created_files: list[
        Path
    ] = list(
        export_tables(
            monthly_common=monthly_common,
            zero_lag=zero_lag,
            lag_curves=lag_curves,
            lag_summary=lag_summary,
            robustness=robustness,
            hypothesis=hypothesis,
        )
    )

    for path in created_files:
        print(
            path
        )

    section(
        "## GENERATING FIGURES"
    )

    figure_groups = (
        plot_zero_lag_heatmaps(
            zero_lag
        ),
        plot_key_band_nino34_lag_curves(
            lag_curves
        ),
        plot_annual_band_time_series(
            monthly_common
        ),
    )

    for group in figure_groups:
        created_files.extend(
            group
        )

        for path in group:
            print(
                path
            )

    print_primary_summary(
        hypothesis=hypothesis,
        lag_summary=lag_summary,
    )

    section(
        "## WRITING SCIENTIFIC REPORTS"
    )

    report_txt, report_json = write_reports(
        band_definitions=band_definitions,
        monthly_common=monthly_common,
        zero_lag=zero_lag,
        lag_summary=lag_summary,
        robustness=robustness,
        hypothesis=hypothesis,
        created_files=created_files,
    )

    created_files.extend(
        [
            report_txt,
            report_json,
        ]
    )

    print(
        report_txt
    )
    print(
        report_json
    )

    print()
    rule()
    print(
        "PROGRAM 33 COMPLETED SUCCESSFULLY."
    )
    print()
    print(
        "Final Program-32 annual/interannual wavelet-band variance was tested "
        "against monthly ENSO state and magnitude for 28.0, 28.5 and 29.0 °C."
    )
    print(
        "Zero-lag inference includes both within-band and conservative "
        "global-band FDR control."
    )
    print(
        "Lagged relationships use calendar-month alignment and are treated as "
        "exploratory temporal associations, not causal evidence."
    )
    print()
    print(
        "AUDIT 06 can now be closed after scientific interpretation of the "
        "Program-33 tables/report."
    )
    rule()


if __name__ == "__main__":
    main()
