#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
===============================================================================
PROJECT
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

PROGRAM
    32_pwp_longitude_wavelet_band_variance_analysis.py

TITLE
    PWP Centroid Longitude — Time-Varying Wavelet Variance Across
    Annual and Interannual Bands

VERSION
    1.0.0

AUTHOR
    Fabio Vieira Machado

PYTHON
    Python 3.10+

===============================================================================
PURPOSE
===============================================================================

    Program 32 quantifies how the variance of Pacific Warm Pool (PWP) centroid
    longitude evolves through time across a hierarchy of annual and interannual
    wavelet bands.

    The program is designed specifically for AUDIT 05 of the PWP project.

    It DOES NOT recompute the Continuous Wavelet Transform (CWT).

    Instead, it reads the authoritative complex-Morlet CWT archives produced
    by Program 11:

        src_a/11_wavelet_analysis.py

    and reconstructs scale-averaged wavelet power directly from the archived
    physical wavelet-power matrix.

    The analysis is performed independently for:

        28.0 °C   — PRIMARY PWP DEFINITION
        28.5 °C   — THRESHOLD SENSITIVITY
        29.0 °C   — WARM-CORE / THRESHOLD SENSITIVITY

===============================================================================
SCIENTIFIC MOTIVATION
===============================================================================

    Fourier/Welch analysis describes how variance is distributed globally
    across frequency over the complete record, but it does not identify when a
    particular frequency band becomes stronger or weaker.

    The Continuous Wavelet Transform preserves both time and scale.

    Program 32 therefore addresses the question:

        How does the variance of PWP zonal displacement evolve through time
        across annual, quasi-biennial, and interannual/ENSO-related bands?

    A second question concerns threshold robustness:

        Does the time-frequency structure persist when the operational PWP SST
        threshold is changed from 28.0 to 28.5 or 29.0 °C?

    This is scientifically important because Programs 30 and 31 showed that the
    SST threshold materially affects PWP geometry and spatial connectivity.
    Program 32 tests whether those geometric differences artificially generate,
    or merely modify, the time-frequency structure of centroid longitude.

===============================================================================
PRIMARY VARIABLE
===============================================================================

    PWP centroid longitude:

        lon_360(t)

    expressed in degrees east on [0, 360).

    Program 32 analyzes only centroid longitude.

    Latitude, area, and area-weighted mean PWP SST remain available elsewhere
    in the pipeline but are not recomputed here.

===============================================================================
OFFICIAL WAVELET BANDS USED BY PROGRAM 32
===============================================================================

    1. ANNUAL BROAD

        40–64 weeks
        = 280–448 days

        Purpose:
            broad representation of variance surrounding the annual cycle.

    2. ANNUAL STRICT

        48–54 weeks
        = 336–378 days

        Purpose:
            narrow representation concentrated close to the 52-week annual
            period.

    3. QUASI-BIENNIAL

        550–1096 days
        approximately 1.51–3.00 years

        Purpose:
            broad low-interannual / quasi-biennial variability.

    4. INTERANNUAL 1

        805–980 days
        approximately 2.20–2.68 years

        Purpose:
            narrower lower-interannual band.

    5. INTERANNUAL 2

        1280–1642 days
        approximately 3.50–4.50 years

        Purpose:
            narrower band centred on the approximately four-year variability
            identified by the global Fourier/Welch analysis.

    6. INTERANNUAL 3

        805–1642 days
        approximately 2.20–4.50 years

        Purpose:
            integrated interannual band spanning Interannual 1 and
            Interannual 2.

    7. ENSO BROAD

        1095–2922 days
        approximately 3.00–8.00 years

        Purpose:
            broad ENSO-scale variance.

    IMPORTANT:
        These bands intentionally overlap.

        The overlapping definitions are diagnostic, not mutually exclusive
        partitions of total variance. Each band answers a different scientific
        question about the concentration and temporal modulation of variance.

===============================================================================
AUTHORITATIVE INPUT
===============================================================================

    Program 11 transform archives:

        data/processed/28/wavelet/transform/
            pwp_wavelet_transform_lon_360.npz

        data/processed/28.5/wavelet/transform/
            pwp_wavelet_transform_lon_360.npz

        data/processed/29/wavelet/transform/
            pwp_wavelet_transform_lon_360.npz

    Program 32 resolves these paths through:

        config.config.get_threshold_paths()

    rather than hard-coding threshold-specific project paths.

===============================================================================
REQUIRED PROGRAM-11 ARCHIVE ARRAYS
===============================================================================

    dates
        Daily observation dates.

    scales_days
        Morlet wavelet scales in days.

    fourier_period_days
        Fourier-equivalent period associated with each Morlet scale.

    coi_days
        Cone-of-influence period at every time step.

    physical_power
        Wavelet power in physical variance units.

    local_significance_ratio
        Local power divided by the Program-11 red-noise significance level.

    Program 32 does not alter these arrays.

===============================================================================
WAVELET METHOD
===============================================================================

    Program 11 uses a complex Morlet mother wavelet.

    The Continuous Wavelet Transform may be represented as:

                  N-1
        W_n(s) =  Σ x_n' ψ* [ (n' - n) Δt / s ]
                  n'=0

    where:

        W_n(s) = wavelet coefficient at time n and scale s
        x_n'   = analyzed time series
        ψ*     = complex conjugate of the scaled/translated mother wavelet
        Δt     = sampling interval
        s      = wavelet scale

    Local wavelet power is:

        |W_n(s)|²

    Program 11 stores physical power after restoring the variance scale of the
    original series.

===============================================================================
SCALE-AVERAGED WAVELET POWER
===============================================================================

    For each diagnostic band B = [s_j1, s_j2], Program 32 reconstructs the
    Torrence–Compo scale-averaged wavelet power using the SAME weighting already
    used by Programs 11 and 20:

                           j2
                           ----
                           \    |W_n(s_j)|²
        P_B(n) = Δj Δt / Cδ  >   -----------
                           /        s_j
                           ----
                           j=j1

    where:

        P_B(n)
            scale-averaged wavelet power for band B at time n;

        Δj
            logarithmic spacing between adjacent wavelet scales;

        Δt
            sampling interval in days;

        Cδ
            reconstruction factor of the Morlet wavelet;

        |W_n(s_j)|²
            PHYSICAL wavelet power archived by Program 11;

        s_j
            wavelet scale in days.

    Because Program 11's physical_power already contains the variance factor,
    Program 32 MUST NOT multiply by the time-series variance again.

    The resulting P_B(t) is a time series describing the evolution of variance
    contained within the selected band.

===============================================================================
CONE OF INFLUENCE
===============================================================================

    Wavelet estimates near the beginning and end of the record are affected by
    edge effects.

    Program 11 archives:

        coi_days(t)

    For each band, Program 32 reports:

        coi_valid_fraction(t)

            fraction of the CWT scales within that diagnostic band whose
            Fourier-equivalent period satisfies:

                period <= coi_days(t)

        strict_coi_valid(t)

            True only when the ENTIRE selected band lies inside the reliable
            cone of influence:

                max(selected band periods) <= coi_days(t)

    Scientific summaries and primary time-series figures use only
    strict_coi_valid observations.

    Raw scale-averaged power is also retained in the output CSV for audit and
    reproducibility, but edge-contaminated values are not used as primary
    evidence.

===============================================================================
LOCAL SIGNIFICANCE
===============================================================================

    Program 11 archives local_significance_ratio:

        R_n(s) = wavelet_power / local_red_noise_significance

    Program 32 calculates the mean local-significance ratio across the selected
    scales for each band and time:

        R_B(t) = mean_s [ R_t(s) ]

    It also reports the fraction of selected scales whose local significance
    ratio is >= 1.

    IMPORTANT:
        Scale-averaged band power and local pointwise significance are not the
        same statistical object.

        Program 32 therefore exports these quantities as diagnostics and does
        not claim that mean local-significance ratio >= 1 is equivalent to a
        formal Torrence–Compo scale-averaged significance test.

===============================================================================
THRESHOLD ROBUSTNESS
===============================================================================

    Threshold robustness is evaluated by comparing the strict-COI-valid
    band-power time series among:

        28.0 °C
        28.5 °C
        29.0 °C

    For every band and threshold pair, Program 32 reports:

        Pearson correlation
        Spearman correlation
        number of common valid dates

    High temporal agreement across thresholds supports the interpretation that
    the identified band modulation is not artificially generated by the
    operational PWP threshold.

    Correlation is evidence of robustness/association, not causality.

===============================================================================
INPUT PROGRAMS / DATA PROVENANCE
===============================================================================

    PWP centroid and area:
        src_a/05_calculate_pwp_centroid.py

    Quality control:
        src_a/06_quality_control_and_diagnostics.py

    Authoritative CWT:
        src_a/11_wavelet_analysis.py

    Scientific spectral synthesis:
        src_a/20_scientific_interpretation_and_spectral_synthesis.py

    Program 20 independently reconstructs diagnostic band time series from the
    same Program-11 archive using:

        Δj * Δt / Cδ * sum(physical_power / scale)

    Program 32 intentionally follows the same formulation.

===============================================================================
OUTPUT DIRECTORIES
===============================================================================

    TABLES

        outputs/tables/threshold_comparison/
            pwp_longitude_wavelet_band_variance/

    FIGURES

        outputs/figures/threshold_comparison/
            pwp_longitude_wavelet_band_variance/

    REPORTS

        outputs/reports/threshold_comparison/
            pwp_longitude_wavelet_band_variance/

===============================================================================
OUTPUT TABLES
===============================================================================

    1. pwp_longitude_wavelet_band_power_daily.csv

       Long-format daily table:

            date
            threshold_c
            band_key
            band_label
            minimum_period_days
            maximum_period_days
            scale_count
            band_power_raw
            band_power_strict_coi
            coi_valid_fraction
            strict_coi_valid
            mean_local_significance_ratio
            significant_scale_fraction

    2. pwp_longitude_wavelet_band_summary.csv

       Per threshold × band:

            number of records
            strict-COI-valid records
            valid fraction
            mean
            median
            standard deviation
            P05
            P25
            P75
            P95
            P99
            maximum
            date of maximum

    3. pwp_longitude_wavelet_band_threshold_robustness.csv

       Pairwise threshold comparison for every band.

    4. pwp_longitude_wavelet_band_definitions.csv

       Machine-readable record of all official Program-32 bands.

===============================================================================
OUTPUT FIGURES
===============================================================================

    1. pwp_longitude_wavelet_bands_28C.png/pdf

       Seven-panel primary figure for SST >= 28 °C.

       A) Annual Broad
       B) Annual Strict
       C) Quasi-biennial
       D) Interannual 1
       E) Interannual 2
       F) Interannual 3
       G) ENSO Broad

    2. pwp_longitude_wavelet_band_threshold_comparison.png/pdf

       Seven-panel threshold-robustness figure.
       Each panel overlays 28.0, 28.5, and 29.0 °C.

    3. pwp_longitude_wavelet_annual_vs_interannual_28C.png/pdf

       Compact synthesis showing:

            Annual Broad
            Annual Strict
            Interannual 3
            ENSO Broad

       for the primary 28 °C definition.

===============================================================================
OUTPUT REPORT
===============================================================================

    PROGRAM32_PWP_LONGITUDE_WAVELET_BAND_VARIANCE.txt

    Includes:

        configuration
        authoritative input files
        band definitions
        scale counts
        strict COI coverage
        summary statistics
        dates and values of maximum valid band power
        threshold robustness
        interpretive safeguards

===============================================================================
FIGURE POLICY
===============================================================================

    - Publication-oriented figures.
    - Transparent legends.
    - Four-year x-axis ticks.
    - X-axis date labels shown only on the bottom panel of multi-panel figures.
    - Panel labels A), B), C), ...
    - No smoothing of band-power time series unless explicitly documented.
    - Edge-contaminated observations are excluded from primary plotted lines.
    - Threshold-comparison panels retain raw physical power units.
    - No normalization across thresholds in the primary figures.

===============================================================================
SCIENTIFIC SAFEGUARDS
===============================================================================

    1. Program 32 does not recompute the CWT.

    2. Band limits are defined a priori in this program and are not moved to
       maximize agreement with ENSO events.

    3. Wavelet power identifies time-varying variance, not causal forcing.

    4. Temporal coincidence with El Niño or La Niña does not prove ENSO
       causality.

    5. The broad ENSO band intentionally overlaps other interannual bands.

    6. Fourier/Welch and global wavelet spectra answer global-frequency
       questions; Program 32 answers when variance in predefined bands changes.

    7. All main quantitative summaries exclude times where the complete band is
       outside the reliable cone of influence.

===============================================================================
DEPENDENCIES
===============================================================================

    Standard library:
        dataclasses
        pathlib
        sys
        typing

    Third-party:
        numpy
        pandas
        scipy
        matplotlib
        pycwt

    Project:
        config.config
        config.wavelet_config

===============================================================================
EXECUTION
===============================================================================

    Recommended from project root:

        python src/threshold_analysis/32_pwp_longitude_wavelet_band_variance_analysis.py

===============================================================================
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

try:
    import pycwt as wavelet
except ImportError as error:
    raise ImportError(
        "Program 32 requires pycwt because the Morlet reconstruction factor "
        "C_delta must match the Program-11 wavelet definition."
    ) from error


# =============================================================================
# PROJECT ROOT
# =============================================================================

SCRIPT_FILE = Path(__file__).resolve()
ROOT = SCRIPT_FILE.parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# =============================================================================
# PROJECT CONFIGURATION
# =============================================================================

from config.config import (  # noqa: E402
    FIGURE_DPI,
    GRID_ALPHA,
    GRID_LINESTYLE,
    GRID_LINEWIDTH,
    PANEL_LABEL_FONT_SIZE,
    PANEL_LABEL_FONT_WEIGHT,
    PANEL_LABEL_X,
    PANEL_LABEL_Y,
    PROJECT_DIR,
    SAVE_BBOX,
    SAVE_PAD_INCHES,
    SAVE_TRANSPARENT,
    YEAR_TICK_INTERVAL,
    get_threshold_paths,
    validate_project_configuration,
)

from config.wavelet_config import (  # noqa: E402
    MORLET_OMEGA0,
    SAMPLING_INTERVAL_DAYS,
    SCALE_SPACING_DJ,
    validate_wavelet_configuration,
)


# =============================================================================
# PROGRAM CONFIGURATION
# =============================================================================

PROGRAM_NAME = (
    "OSAF PROGRAM 32 — PWP LONGITUDE WAVELET BAND-VARIANCE ANALYSIS"
)
PROGRAM_VERSION = "1.0.0"

VARIABLE = "lon_360"
VARIABLE_LABEL = "PWP centroid longitude"

THRESHOLDS_C = (28.0, 28.5, 29.0)
PRIMARY_THRESHOLD_C = 28.0

WEEKS_TO_DAYS = 7.0

TABLE_DIR = (
    ROOT
    / "outputs"
    / "tables"
    / "threshold_comparison"
    / "pwp_longitude_wavelet_band_variance"
)

FIGURE_DIR = (
    ROOT
    / "outputs"
    / "figures"
    / "threshold_comparison"
    / "pwp_longitude_wavelet_band_variance"
)

REPORT_DIR = (
    ROOT
    / "outputs"
    / "reports"
    / "threshold_comparison"
    / "pwp_longitude_wavelet_band_variance"
)

DAILY_TABLE = TABLE_DIR / "pwp_longitude_wavelet_band_power_daily.csv"
SUMMARY_TABLE = TABLE_DIR / "pwp_longitude_wavelet_band_summary.csv"
ROBUSTNESS_TABLE = (
    TABLE_DIR
    / "pwp_longitude_wavelet_band_threshold_robustness.csv"
)
BAND_DEFINITION_TABLE = (
    TABLE_DIR
    / "pwp_longitude_wavelet_band_definitions.csv"
)

REPORT_FILE = (
    REPORT_DIR
    / "PROGRAM32_PWP_LONGITUDE_WAVELET_BAND_VARIANCE.txt"
)


# =============================================================================
# BAND DEFINITIONS
# =============================================================================

@dataclass(frozen=True)
class Band:
    key: str
    label: str
    minimum_period_days: float
    maximum_period_days: float
    scientific_role: str

    @property
    def period_range_text(self) -> str:
        return (
            f"{self.minimum_period_days:.0f}–"
            f"{self.maximum_period_days:.0f} days"
        )


BANDS: tuple[Band, ...] = (
    Band(
        key="annual_broad",
        label="Annual Broad",
        minimum_period_days=40.0 * WEEKS_TO_DAYS,
        maximum_period_days=64.0 * WEEKS_TO_DAYS,
        scientific_role=(
            "Broad variance surrounding the annual cycle "
            "(40–64 weeks)."
        ),
    ),
    Band(
        key="annual_strict",
        label="Annual Strict",
        minimum_period_days=48.0 * WEEKS_TO_DAYS,
        maximum_period_days=54.0 * WEEKS_TO_DAYS,
        scientific_role=(
            "Narrow annual-cycle variance near 52 weeks "
            "(48–54 weeks)."
        ),
    ),
    Band(
        key="quasi_biennial",
        label="Quasi-biennial",
        minimum_period_days=550.0,
        maximum_period_days=1096.0,
        scientific_role=(
            "Broad quasi-biennial / lower-interannual variability."
        ),
    ),
    Band(
        key="interannual_1",
        label="Interannual 1",
        minimum_period_days=805.0,
        maximum_period_days=980.0,
        scientific_role=(
            "Narrow lower-interannual band (~2.2–2.7 years)."
        ),
    ),
    Band(
        key="interannual_2",
        label="Interannual 2",
        minimum_period_days=1280.0,
        maximum_period_days=1642.0,
        scientific_role=(
            "Narrow ~3.5–4.5-year band surrounding the "
            "approximately four-year Fourier/Welch feature."
        ),
    ),
    Band(
        key="interannual_3",
        label="Interannual 3",
        minimum_period_days=805.0,
        maximum_period_days=1642.0,
        scientific_role=(
            "Integrated ~2.2–4.5-year interannual band."
        ),
    ),
    Band(
        key="enso_broad",
        label="ENSO Broad",
        minimum_period_days=1095.0,
        maximum_period_days=2922.0,
        scientific_role=(
            "Broad ~3–8-year ENSO-scale variance."
        ),
    ),
)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class TransformArchive:
    threshold_c: float
    path: Path
    dates: pd.DatetimeIndex
    scales_days: np.ndarray
    periods_days: np.ndarray
    coi_days: np.ndarray
    physical_power: np.ndarray
    local_significance_ratio: np.ndarray


# =============================================================================
# HELPERS
# =============================================================================

def threshold_slug(threshold_c: float) -> str:
    return f"{threshold_c:.1f}".replace(".", "p")


def add_panel_label(axis: plt.Axes, text: str) -> None:
    axis.text(
        PANEL_LABEL_X,
        PANEL_LABEL_Y,
        text,
        transform=axis.transAxes,
        fontsize=PANEL_LABEL_FONT_SIZE,
        fontweight=PANEL_LABEL_FONT_WEIGHT,
        va="top",
        ha="left",
    )


def configure_date_axis(
    axis: plt.Axes,
    show_labels: bool,
) -> None:
    axis.xaxis.set_major_locator(
        mdates.YearLocator(base=YEAR_TICK_INTERVAL)
    )
    axis.xaxis.set_major_formatter(
        mdates.DateFormatter("%Y")
    )

    if not show_labels:
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

    png = FIGURE_DIR / f"{stem}.png"
    pdf = FIGURE_DIR / f"{stem}.pdf"

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

    plt.close(figure)

    return png, pdf


def program11_transform_path(threshold_c: float) -> Path:
    paths = get_threshold_paths(
        float(threshold_c)
    )

    return (
        paths.wavelet.processed_dir
        / "transform"
        / f"pwp_wavelet_transform_{VARIABLE}.npz"
    )


# =============================================================================
# VALIDATION
# =============================================================================

def validate_band_definitions() -> None:
    keys = [
        band.key
        for band in BANDS
    ]

    if len(keys) != len(set(keys)):
        raise ValueError(
            "Program-32 wavelet band keys must be unique."
        )

    for band in BANDS:
        if (
            not np.isfinite(
                band.minimum_period_days
            )
            or not np.isfinite(
                band.maximum_period_days
            )
            or band.minimum_period_days <= 0.0
            or band.maximum_period_days
            <= band.minimum_period_days
        ):
            raise ValueError(
                f"Invalid band definition: {band}"
            )


def validate_project() -> None:
    if Path(PROJECT_DIR).resolve() != ROOT.resolve():
        raise ValueError(
            "Project-root mismatch:\n"
            f"  Config : {PROJECT_DIR}\n"
            f"  Script : {ROOT}"
        )

    validate_project_configuration()
    validate_wavelet_configuration()
    validate_band_definitions()

    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


# =============================================================================
# PROGRAM-11 ARCHIVE LOADING
# =============================================================================

def load_transform_archive(
    threshold_c: float,
) -> TransformArchive:
    path = program11_transform_path(
        threshold_c
    )

    if not path.is_file():
        raise FileNotFoundError(
            "Authoritative Program-11 transform archive was not found:\n"
            f"{path}"
        )

    with np.load(
        path,
        allow_pickle=False,
    ) as archive:
        required = (
            "dates",
            "scales_days",
            "fourier_period_days",
            "coi_days",
            "physical_power",
            "local_significance_ratio",
        )

        missing = [
            key
            for key in required
            if key not in archive.files
        ]

        if missing:
            raise KeyError(
                f"{path} is missing required arrays: {missing}\n"
                f"Available arrays: {archive.files}"
            )

        dates = pd.DatetimeIndex(
            pd.to_datetime(
                archive["dates"]
            )
        )

        scales = np.asarray(
            archive["scales_days"],
            dtype=float,
        )

        periods = np.asarray(
            archive["fourier_period_days"],
            dtype=float,
        )

        coi = np.asarray(
            archive["coi_days"],
            dtype=float,
        )

        physical_power = np.asarray(
            archive["physical_power"],
            dtype=float,
        )

        local_ratio = np.asarray(
            archive["local_significance_ratio"],
            dtype=float,
        )

    if dates.has_duplicates:
        raise ValueError(
            f"Duplicate dates in Program-11 archive: {path}"
        )

    if not dates.is_monotonic_increasing:
        raise ValueError(
            f"Dates are not monotonic in Program-11 archive: {path}"
        )

    if scales.ndim != 1 or periods.ndim != 1 or coi.ndim != 1:
        raise ValueError(
            f"Unexpected Program-11 coordinate dimensions: {path}"
        )

    expected_shape = (
        scales.size,
        dates.size,
    )

    if physical_power.shape != expected_shape:
        raise ValueError(
            "physical_power has unexpected shape.\n"
            f"Expected: {expected_shape}\n"
            f"Found   : {physical_power.shape}"
        )

    if local_ratio.shape != expected_shape:
        raise ValueError(
            "local_significance_ratio has unexpected shape.\n"
            f"Expected: {expected_shape}\n"
            f"Found   : {local_ratio.shape}"
        )

    if periods.size != scales.size:
        raise ValueError(
            "Program-11 scales and periods have inconsistent lengths."
        )

    if coi.size != dates.size:
        raise ValueError(
            "Program-11 COI and dates have inconsistent lengths."
        )

    if np.any(~np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError(
            f"Invalid wavelet scales in {path}"
        )

    if np.any(~np.isfinite(periods)) or np.any(periods <= 0.0):
        raise ValueError(
            f"Invalid Fourier-equivalent periods in {path}"
        )

    if np.any(~np.isfinite(coi)) or np.any(coi <= 0.0):
        raise ValueError(
            f"Invalid cone-of-influence periods in {path}"
        )

    return TransformArchive(
        threshold_c=float(threshold_c),
        path=path,
        dates=dates,
        scales_days=scales,
        periods_days=periods,
        coi_days=coi,
        physical_power=physical_power,
        local_significance_ratio=local_ratio,
    )


# =============================================================================
# SCALE-AVERAGED BAND POWER
# =============================================================================

def band_indices(
    periods_days: np.ndarray,
    band: Band,
) -> np.ndarray:
    indices = np.flatnonzero(
        (
            periods_days
            >= band.minimum_period_days
        )
        & (
            periods_days
            <= band.maximum_period_days
        )
    )

    if indices.size < 2:
        nearest = periods_days[
            np.argsort(
                np.abs(
                    periods_days
                    - (
                        band.minimum_period_days
                        + band.maximum_period_days
                    )
                    / 2.0
                )
            )[:6]
        ]

        raise ValueError(
            f"Band {band.label} ({band.period_range_text}) "
            "contains fewer than two Program-11 CWT scales.\n"
            f"Nearest available periods: {np.sort(nearest)}"
        )

    return indices


def calculate_band_timeseries(
    archive: TransformArchive,
    band: Band,
    cdelta: float,
) -> pd.DataFrame:
    indices = band_indices(
        archive.periods_days,
        band,
    )

    selected_power = (
        archive.physical_power[
            indices,
            :
        ]
    )

    selected_scales = (
        archive.scales_days[
            indices
        ]
    )

    selected_periods = (
        archive.periods_days[
            indices
        ]
    )

    selected_local_ratio = (
        archive.local_significance_ratio[
            indices,
            :
        ]
    )

    # Torrence-Compo scale-averaged physical wavelet power.
    band_power = (
        SCALE_SPACING_DJ
        * SAMPLING_INTERVAL_DAYS
        / cdelta
        * np.sum(
            selected_power
            / selected_scales[
                :,
                np.newaxis,
            ],
            axis=0,
        )
    )

    # A scale is within the reliable COI when its period is <= COI(t).
    scale_is_coi_valid = (
        selected_periods[
            :,
            np.newaxis
        ]
        <= archive.coi_days[
            np.newaxis,
            :
        ]
    )

    coi_valid_fraction = np.mean(
        scale_is_coi_valid,
        axis=0,
    )

    strict_coi_valid = np.all(
        scale_is_coi_valid,
        axis=0,
    )

    mean_local_ratio = np.nanmean(
        selected_local_ratio,
        axis=0,
    )

    significant_scale_fraction = np.mean(
        selected_local_ratio >= 1.0,
        axis=0,
    )

    strict_power = np.where(
        strict_coi_valid,
        band_power,
        np.nan,
    )

    return pd.DataFrame(
        {
            "date": archive.dates,
            "threshold_c": (
                archive.threshold_c
            ),
            "band_key": band.key,
            "band_label": band.label,
            "minimum_period_days": (
                band.minimum_period_days
            ),
            "maximum_period_days": (
                band.maximum_period_days
            ),
            "scale_count": int(
                indices.size
            ),
            "minimum_selected_period_days": float(
                np.min(
                    selected_periods
                )
            ),
            "maximum_selected_period_days": float(
                np.max(
                    selected_periods
                )
            ),
            "band_power_raw": band_power,
            "band_power_strict_coi": (
                strict_power
            ),
            "coi_valid_fraction": (
                coi_valid_fraction
            ),
            "strict_coi_valid": (
                strict_coi_valid
            ),
            "mean_local_significance_ratio": (
                mean_local_ratio
            ),
            "significant_scale_fraction": (
                significant_scale_fraction
            ),
        }
    )


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def build_summary(
    daily: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []

    grouped = daily.groupby(
        [
            "threshold_c",
            "band_key",
            "band_label",
            "minimum_period_days",
            "maximum_period_days",
        ],
        sort=True,
        observed=True,
    )

    for keys, group in grouped:
        (
            threshold_c,
            band_key,
            band_label,
            minimum_period_days,
            maximum_period_days,
        ) = keys

        valid = group.loc[
            group["strict_coi_valid"].astype(bool)
            & group["band_power_strict_coi"].notna()
        ].copy()

        if valid.empty:
            raise ValueError(
                f"No strict-COI-valid observations for "
                f"{threshold_c:.1f} °C / {band_label}."
            )

        values = valid[
            "band_power_strict_coi"
        ].astype(float)

        max_index = values.idxmax()
        max_row = valid.loc[
            max_index
        ]

        rows.append(
            {
                "threshold_c": float(
                    threshold_c
                ),
                "band_key": str(
                    band_key
                ),
                "band_label": str(
                    band_label
                ),
                "minimum_period_days": float(
                    minimum_period_days
                ),
                "maximum_period_days": float(
                    maximum_period_days
                ),
                "scale_count": int(
                    group[
                        "scale_count"
                    ].iloc[0]
                ),
                "n_total": int(
                    len(group)
                ),
                "n_strict_coi_valid": int(
                    len(valid)
                ),
                "strict_coi_valid_pct": float(
                    100.0
                    * len(valid)
                    / len(group)
                ),
                "mean_power": float(
                    values.mean()
                ),
                "median_power": float(
                    values.median()
                ),
                "std_power": float(
                    values.std(ddof=1)
                ),
                "p05_power": float(
                    values.quantile(0.05)
                ),
                "p25_power": float(
                    values.quantile(0.25)
                ),
                "p75_power": float(
                    values.quantile(0.75)
                ),
                "p95_power": float(
                    values.quantile(0.95)
                ),
                "p99_power": float(
                    values.quantile(0.99)
                ),
                "maximum_power": float(
                    max_row[
                        "band_power_strict_coi"
                    ]
                ),
                "maximum_power_date": (
                    pd.Timestamp(
                        max_row["date"]
                    ).strftime(
                        "%Y-%m-%d"
                    )
                ),
                "mean_local_significance_ratio": float(
                    valid[
                        "mean_local_significance_ratio"
                    ].mean()
                ),
                "mean_significant_scale_fraction": float(
                    valid[
                        "significant_scale_fraction"
                    ].mean()
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
                "minimum_period_days",
                "maximum_period_days",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# THRESHOLD ROBUSTNESS
# =============================================================================

def safe_correlation(
    x: pd.Series,
    y: pd.Series,
) -> tuple[
    int,
    float,
    float,
    float,
    float,
]:
    frame = pd.DataFrame(
        {
            "x": pd.to_numeric(
                x,
                errors="coerce",
            ),
            "y": pd.to_numeric(
                y,
                errors="coerce",
            ),
        }
    ).dropna()

    n = int(
        len(frame)
    )

    if n < 3:
        return (
            n,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    pearson_r, pearson_p = pearsonr(
        frame["x"],
        frame["y"],
    )

    spearman_rho, spearman_p = spearmanr(
        frame["x"],
        frame["y"],
    )

    return (
        n,
        float(pearson_r),
        float(pearson_p),
        float(spearman_rho),
        float(spearman_p),
    )


def build_threshold_robustness(
    daily: pd.DataFrame,
) -> pd.DataFrame:
    pairs = (
        (28.0, 28.5),
        (28.0, 29.0),
        (28.5, 29.0),
    )

    rows: list[dict] = []

    for band in BANDS:
        subset = daily.loc[
            daily["band_key"]
            == band.key
        ].copy()

        for threshold_a, threshold_b in pairs:
            a = subset.loc[
                subset["threshold_c"]
                == threshold_a,
                [
                    "date",
                    "band_power_strict_coi",
                ],
            ].rename(
                columns={
                    "band_power_strict_coi":
                    "power_a",
                }
            )

            b = subset.loc[
                subset["threshold_c"]
                == threshold_b,
                [
                    "date",
                    "band_power_strict_coi",
                ],
            ].rename(
                columns={
                    "band_power_strict_coi":
                    "power_b",
                }
            )

            merged = a.merge(
                b,
                on="date",
                how="inner",
                validate="one_to_one",
            ).dropna(
                subset=[
                    "power_a",
                    "power_b",
                ]
            )

            (
                n,
                pearson_r,
                pearson_p,
                spearman_rho,
                spearman_p,
            ) = safe_correlation(
                merged["power_a"],
                merged["power_b"],
            )

            rows.append(
                {
                    "band_key": band.key,
                    "band_label": band.label,
                    "minimum_period_days": (
                        band.minimum_period_days
                    ),
                    "maximum_period_days": (
                        band.maximum_period_days
                    ),
                    "threshold_a_c": (
                        threshold_a
                    ),
                    "threshold_b_c": (
                        threshold_b
                    ),
                    "n_common_strict_coi": n,
                    "pearson_r": (
                        pearson_r
                    ),
                    "pearson_p": (
                        pearson_p
                    ),
                    "spearman_rho": (
                        spearman_rho
                    ),
                    "spearman_p": (
                        spearman_p
                    ),
                }
            )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# OUTPUT TABLES
# =============================================================================

def band_definition_frame() -> pd.DataFrame:
    rows: list[dict] = []

    for order, band in enumerate(
        BANDS,
        start=1,
    ):
        rows.append(
            {
                "order": order,
                "band_key": band.key,
                "band_label": band.label,
                "minimum_period_days": (
                    band.minimum_period_days
                ),
                "maximum_period_days": (
                    band.maximum_period_days
                ),
                "minimum_period_weeks": (
                    band.minimum_period_days
                    / 7.0
                ),
                "maximum_period_weeks": (
                    band.maximum_period_days
                    / 7.0
                ),
                "minimum_period_years": (
                    band.minimum_period_days
                    / 365.2425
                ),
                "maximum_period_years": (
                    band.maximum_period_days
                    / 365.2425
                ),
                "scientific_role": (
                    band.scientific_role
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


def export_tables(
    daily: pd.DataFrame,
    summary: pd.DataFrame,
    robustness: pd.DataFrame,
) -> tuple[Path, ...]:
    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    daily_out = daily.copy()

    daily_out["date"] = pd.to_datetime(
        daily_out["date"]
    ).dt.strftime(
        "%Y-%m-%d"
    )

    daily_out.to_csv(
        DAILY_TABLE,
        index=False,
        float_format="%.10g",
    )

    summary.to_csv(
        SUMMARY_TABLE,
        index=False,
        float_format="%.10g",
    )

    robustness.to_csv(
        ROBUSTNESS_TABLE,
        index=False,
        float_format="%.10g",
    )

    band_definition_frame().to_csv(
        BAND_DEFINITION_TABLE,
        index=False,
        float_format="%.10g",
    )

    return (
        DAILY_TABLE,
        SUMMARY_TABLE,
        ROBUSTNESS_TABLE,
        BAND_DEFINITION_TABLE,
    )


# =============================================================================
# FIGURES
# =============================================================================

def plot_primary_bands(
    daily: pd.DataFrame,
) -> tuple[Path, Path]:
    primary = daily.loc[
        daily["threshold_c"]
        == PRIMARY_THRESHOLD_C
    ].copy()

    figure, axes = plt.subplots(
        len(BANDS),
        1,
        figsize=(
            15,
            2.65 * len(BANDS),
        ),
        sharex=True,
        constrained_layout=True,
    )

    for index, (
        axis,
        band,
    ) in enumerate(
        zip(
            axes,
            BANDS,
        )
    ):
        group = (
            primary.loc[
                primary[
                    "band_key"
                ]
                == band.key
            ]
            .sort_values(
                "date"
            )
        )

        axis.plot(
            group["date"],
            group[
                "band_power_strict_coi"
            ],
            linewidth=0.8,
            label=(
                f"{band.label} "
                f"({band.period_range_text})"
            ),
        )

        axis.set_ylabel(
            "Scale-averaged\npower"
        )

        axis.set_title(
            (
                f"{band.label}: "
                f"{band.period_range_text}"
            ),
            loc="left",
            fontsize=10,
        )

        axis.legend(
            frameon=False,
            loc="upper right",
            fontsize=8,
        )

        configure_date_axis(
            axis,
            show_labels=(
                index
                == len(BANDS) - 1
            ),
        )

        add_panel_label(
            axis,
            f"{chr(65 + index)})",
        )

    axes[-1].set_xlabel(
        "Year"
    )

    figure.suptitle(
        (
            "Time-varying wavelet variance of PWP centroid longitude "
            "— SST ≥ 28 °C"
        ),
        fontsize=14,
    )

    return save_figure(
        figure,
        "pwp_longitude_wavelet_bands_28C",
    )


def plot_threshold_comparison(
    daily: pd.DataFrame,
) -> tuple[Path, Path]:
    figure, axes = plt.subplots(
        len(BANDS),
        1,
        figsize=(
            15,
            2.75 * len(BANDS),
        ),
        sharex=True,
        constrained_layout=True,
    )

    for index, (
        axis,
        band,
    ) in enumerate(
        zip(
            axes,
            BANDS,
        )
    ):
        subset = daily.loc[
            daily["band_key"]
            == band.key
        ]

        for threshold in THRESHOLDS_C:
            group = (
                subset.loc[
                    subset["threshold_c"]
                    == threshold
                ]
                .sort_values(
                    "date"
                )
            )

            axis.plot(
                group["date"],
                group[
                    "band_power_strict_coi"
                ],
                linewidth=0.75,
                label=(
                    f"{threshold:.1f} °C"
                ),
            )

        axis.set_ylabel(
            "Scale-averaged\npower"
        )

        axis.set_title(
            (
                f"{band.label}: "
                f"{band.period_range_text}"
            ),
            loc="left",
            fontsize=10,
        )

        configure_date_axis(
            axis,
            show_labels=(
                index
                == len(BANDS) - 1
            ),
        )

        add_panel_label(
            axis,
            f"{chr(65 + index)})",
        )

        if index == 0:
            axis.legend(
                frameon=False,
                ncol=3,
                loc="upper right",
            )

    axes[-1].set_xlabel(
        "Year"
    )

    figure.suptitle(
        (
            "Threshold robustness of time-varying PWP longitude "
            "wavelet variance"
        ),
        fontsize=14,
    )

    return save_figure(
        figure,
        "pwp_longitude_wavelet_band_threshold_comparison",
    )


def plot_compact_primary_synthesis(
    daily: pd.DataFrame,
) -> tuple[Path, Path]:
    keys = (
        "annual_broad",
        "annual_strict",
        "interannual_3",
        "enso_broad",
    )

    selected_bands = [
        next(
            band
            for band in BANDS
            if band.key == key
        )
        for key in keys
    ]

    primary = daily.loc[
        daily["threshold_c"]
        == PRIMARY_THRESHOLD_C
    ].copy()

    figure, axes = plt.subplots(
        len(selected_bands),
        1,
        figsize=(15, 11),
        sharex=True,
        constrained_layout=True,
    )

    for index, (
        axis,
        band,
    ) in enumerate(
        zip(
            axes,
            selected_bands,
        )
    ):
        group = (
            primary.loc[
                primary[
                    "band_key"
                ]
                == band.key
            ]
            .sort_values(
                "date"
            )
        )

        axis.plot(
            group["date"],
            group[
                "band_power_strict_coi"
            ],
            linewidth=0.9,
            label=band.label,
        )

        axis.set_ylabel(
            "Scale-averaged\npower"
        )

        axis.set_title(
            (
                f"{band.label}: "
                f"{band.period_range_text}"
            ),
            loc="left",
            fontsize=10,
        )

        axis.legend(
            frameon=False,
            loc="upper right",
        )

        configure_date_axis(
            axis,
            show_labels=(
                index
                == len(
                    selected_bands
                ) - 1
            ),
        )

        add_panel_label(
            axis,
            f"{chr(65 + index)})",
        )

    axes[-1].set_xlabel(
        "Year"
    )

    figure.suptitle(
        (
            "Annual and interannual wavelet-variance modulation of "
            "PWP centroid longitude — SST ≥ 28 °C"
        ),
        fontsize=14,
    )

    return save_figure(
        figure,
        "pwp_longitude_wavelet_annual_vs_interannual_28C",
    )


# =============================================================================
# REPORT
# =============================================================================

def write_report(
    archives: dict[
        float,
        TransformArchive,
    ],
    summary: pd.DataFrame,
    robustness: pd.DataFrame,
    created_files: Iterable[
        Path
    ],
) -> Path:
    lines: list[str] = []

    lines.extend(
        [
            PROGRAM_NAME,
            "=" * 78,
            "",
            "1. CONFIGURATION",
            "-" * 78,
            f"Program version              : {PROGRAM_VERSION}",
            f"Project root                 : {ROOT}",
            f"Variable                     : {VARIABLE_LABEL}",
            "Thresholds                   : 28.0, 28.5, 29.0 °C",
            "Primary threshold            : 28.0 °C",
            f"Morlet omega0                : {MORLET_OMEGA0}",
            f"Sampling interval            : {SAMPLING_INTERVAL_DAYS} days",
            f"Scale spacing dj             : {SCALE_SPACING_DJ}",
            "",
            "2. BAND DEFINITIONS",
            "-" * 78,
        ]
    )

    for band in BANDS:
        lines.append(
            f"{band.label:<22s}: "
            f"{band.minimum_period_days:.0f}–"
            f"{band.maximum_period_days:.0f} days"
        )

    lines.extend(
        [
            "",
            "3. AUTHORITATIVE PROGRAM-11 INPUTS",
            "-" * 78,
        ]
    )

    for threshold in THRESHOLDS_C:
        archive = archives[
            threshold
        ]

        lines.append(
            f"{threshold:4.1f} °C : {archive.path}"
        )

        lines.append(
            "           "
            f"{len(archive.dates):,} dates | "
            f"{archive.dates.min():%Y-%m-%d} to "
            f"{archive.dates.max():%Y-%m-%d} | "
            f"{len(archive.scales_days)} CWT scales"
        )

    lines.extend(
        [
            "",
            "4. STRICT-COI SUMMARY",
            "-" * 78,
        ]
    )

    for threshold in THRESHOLDS_C:
        lines.append(
            f"\nSST >= {threshold:.1f} °C"
        )

        subset = summary.loc[
            summary["threshold_c"]
            == threshold
        ]

        for band in BANDS:
            row = subset.loc[
                subset["band_key"]
                == band.key
            ].iloc[0]

            lines.append(
                (
                    f"  {band.label:<20s} | "
                    f"valid={row['strict_coi_valid_pct']:6.2f}% | "
                    f"mean={row['mean_power']:12.5g} | "
                    f"median={row['median_power']:12.5g} | "
                    f"P95={row['p95_power']:12.5g} | "
                    f"max={row['maximum_power']:12.5g} | "
                    f"date={row['maximum_power_date']}"
                )
            )

    lines.extend(
        [
            "",
            "5. THRESHOLD ROBUSTNESS",
            "-" * 78,
        ]
    )

    for band in BANDS:
        lines.append(
            f"\n{band.label}"
        )

        subset = robustness.loc[
            robustness["band_key"]
            == band.key
        ]

        for _, row in subset.iterrows():
            lines.append(
                (
                    f"  {row['threshold_a_c']:.1f} vs "
                    f"{row['threshold_b_c']:.1f} °C | "
                    f"N={int(row['n_common_strict_coi']):,} | "
                    f"Pearson r={row['pearson_r']:.4f} | "
                    f"Spearman rho={row['spearman_rho']:.4f}"
                )
            )

    lines.extend(
        [
            "",
            "6. METHOD",
            "-" * 78,
            (
                "Scale-averaged physical wavelet power was reconstructed "
                "directly from the authoritative Program-11 CWT archive using:"
            ),
            "",
            "    P_B(t) = dj * dt / C_delta * sum_j[ physical_power_j(t) / scale_j ]",
            "",
            (
                "Program-11 physical_power already contains the physical "
                "variance factor and was not multiplied by variance again."
            ),
            (
                "Primary summaries use only dates for which the complete "
                "diagnostic band is inside the Program-11 cone of influence."
            ),
            "",
            "7. INTERPRETIVE SAFEGUARDS",
            "-" * 78,
            (
                "1. Wavelet band power describes time-varying variance; it "
                "does not establish causal forcing."
            ),
            (
                "2. Band limits were defined a priori and were not optimized "
                "against ENSO-event dates."
            ),
            (
                "3. The diagnostic bands intentionally overlap and must not "
                "be interpreted as mutually exclusive variance partitions."
            ),
            (
                "4. Threshold correlation tests temporal robustness of band "
                "modulation, not physical equivalence of PWP geometries."
            ),
            (
                "5. Edge-contaminated observations are retained only for "
                "audit in the raw-power column and excluded from primary "
                "strict-COI summaries."
            ),
            "",
            "8. FILES CREATED",
            "-" * 78,
        ]
    )

    for path in created_files:
        lines.append(
            str(path)
        )

    lines.extend(
        [
            "",
            "=" * 78,
            "PROGRAM 32 COMPLETED SUCCESSFULLY.",
            "=" * 78,
        ]
    )

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_FILE.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )

    return REPORT_FILE


# =============================================================================
# TERMINAL OUTPUT
# =============================================================================

def print_summary(
    summary: pd.DataFrame,
    robustness: pd.DataFrame,
) -> None:
    print()
    print("## PRIMARY 28 °C — STRICT-COI BAND SUMMARY")
    print()

    primary = summary.loc[
        summary["threshold_c"]
        == PRIMARY_THRESHOLD_C
    ]

    for band in BANDS:
        row = primary.loc[
            primary["band_key"]
            == band.key
        ].iloc[0]

        print(
            f"{band.label:<20s} | "
            f"valid={row['strict_coi_valid_pct']:6.2f}% | "
            f"mean={row['mean_power']:10.4f} | "
            f"median={row['median_power']:10.4f} | "
            f"P95={row['p95_power']:10.4f} | "
            f"max={row['maximum_power']:10.4f} | "
            f"{row['maximum_power_date']}"
        )

    print()
    print("## THRESHOLD ROBUSTNESS — 28 °C VS OTHER DEFINITIONS")
    print()

    focus = robustness.loc[
        robustness["threshold_a_c"]
        == 28.0
    ]

    for _, row in focus.iterrows():
        print(
            f"{row['band_label']:<20s} | "
            f"28.0 vs {row['threshold_b_c']:.1f} °C | "
            f"N={int(row['n_common_strict_coi']):,} | "
            f"Pearson r={row['pearson_r']:.4f} | "
            f"Spearman rho={row['spearman_rho']:.4f}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 78)
    print(PROGRAM_NAME)
    print("=" * 78)
    print()

    print("## CONFIGURATION")
    print()
    print(
        f"Program version                         : {PROGRAM_VERSION}"
    )
    print(
        f"Project root                            : {ROOT}"
    )
    print(
        f"Variable                                : {VARIABLE_LABEL}"
    )
    print(
        "Thresholds                              : 28.0, 28.5, 29.0 °C"
    )
    print(
        f"Primary threshold                       : {PRIMARY_THRESHOLD_C:.1f} °C"
    )
    print(
        f"Morlet omega0                           : {MORLET_OMEGA0}"
    )
    print(
        f"Sampling interval                       : {SAMPLING_INTERVAL_DAYS} days"
    )
    print(
        f"Scale spacing dj                        : {SCALE_SPACING_DJ}"
    )
    print(
        "Program-11 CWT recomputation            : NO"
    )
    print(
        "COI policy                              : complete band required inside COI"
    )

    validate_project()

    mother = wavelet.Morlet(
        MORLET_OMEGA0
    )

    cdelta = float(
        mother.cdelta
    )

    if (
        not np.isfinite(
            cdelta
        )
        or cdelta <= 0.0
    ):
        raise ValueError(
            f"Invalid Morlet reconstruction factor C_delta: {cdelta}"
        )

    print(
        f"Morlet reconstruction factor C_delta    : {cdelta}"
    )

    print()
    print("## OFFICIAL DIAGNOSTIC BANDS")
    print()

    for band in BANDS:
        print(
            f"{band.label:<20s} : "
            f"{band.minimum_period_days:.0f}–"
            f"{band.maximum_period_days:.0f} days"
        )

    print()
    print("## LOADING AUTHORITATIVE PROGRAM-11 TRANSFORMS")
    print()

    archives: dict[
        float,
        TransformArchive,
    ] = {}

    for threshold in THRESHOLDS_C:
        archive = load_transform_archive(
            threshold
        )

        archives[
            threshold
        ] = archive

        print(
            f"{threshold:4.1f} °C : "
            f"{len(archive.dates):,} dates | "
            f"{len(archive.scales_days)} scales | "
            f"{archive.path}"
        )

    # Require identical temporal records across thresholds so direct threshold
    # comparison is scientifically transparent.
    primary_dates = archives[
        PRIMARY_THRESHOLD_C
    ].dates

    for threshold in THRESHOLDS_C:
        if not archives[
            threshold
        ].dates.equals(
            primary_dates
        ):
            raise ValueError(
                "Program-11 transform dates differ among thresholds. "
                "Program 32 requires identical daily records for direct "
                "threshold comparison."
            )

    print()
    print("## CALCULATING SCALE-AVERAGED BAND POWER")
    print()

    frames: list[
        pd.DataFrame
    ] = []

    for threshold in THRESHOLDS_C:
        archive = archives[
            threshold
        ]

        for band in BANDS:
            frame = calculate_band_timeseries(
                archive=archive,
                band=band,
                cdelta=cdelta,
            )

            frames.append(
                frame
            )

            print(
                f"{threshold:4.1f} °C | "
                f"{band.label:<20s} | "
                f"scales={int(frame['scale_count'].iloc[0]):3d} | "
                f"strict COI="
                f"{100.0 * frame['strict_coi_valid'].mean():6.2f}%"
            )

    daily = pd.concat(
        frames,
        ignore_index=True,
    )

    summary = build_summary(
        daily
    )

    robustness = build_threshold_robustness(
        daily
    )

    print()
    print("## EXPORTING TABLES")
    print()

    created_files: list[
        Path
    ] = list(
        export_tables(
            daily=daily,
            summary=summary,
            robustness=robustness,
        )
    )

    for path in created_files:
        print(path)

    print()
    print("## GENERATING FIGURES")
    print()

    for files in (
        plot_primary_bands(
            daily
        ),
        plot_threshold_comparison(
            daily
        ),
        plot_compact_primary_synthesis(
            daily
        ),
    ):
        created_files.extend(
            files
        )

        for path in files:
            print(path)

    print_summary(
        summary=summary,
        robustness=robustness,
    )

    report = write_report(
        archives=archives,
        summary=summary,
        robustness=robustness,
        created_files=created_files,
    )

    created_files.append(
        report
    )

    print()
    print("## REPORT")
    print()
    print(report)

    print()
    print("=" * 78)
    print("PROGRAM 32 COMPLETED SUCCESSFULLY.")
    print()
    print(
        "Seven predefined annual/interannual wavelet bands were reconstructed "
        "from Program-11 physical CWT power for the 28.0, 28.5, and 29.0 °C "
        "PWP definitions."
    )
    print(
        "Primary summaries exclude dates for which the complete band lies "
        "outside the reliable cone of influence."
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

