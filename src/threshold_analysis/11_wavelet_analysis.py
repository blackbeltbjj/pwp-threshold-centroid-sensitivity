# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

PROGRAM
    11_wavelet_analysis.py

VERSION
    4.0.0

PURPOSE
    Perform a complete, configurable, publication-oriented Continuous Wavelet
    Transform analysis of quality-controlled daily Pacific Warm Pool series.

SCIENTIFIC BASIS
    Torrence and Compo (1998); Domingues et al. (2005); Gu and Philander
    (1995); Torrence and Webster (1999).

CONFIGURATION
    Scientific bands and figure families are controlled by:

        config/wavelet_config.py

    Between one and six period bands may be enabled.

FIGURE FAMILIES
    - traditional complete;
    - band multipanel;
    - one figure per band;
    - all bands overlaid;
    - series + scalogram + bands + aligned global spectrum;
    - series + scalogram;
    - isolated global spectrum with period on X;
    - isolated real component.

ALIGNMENT POLICY
    Temporal axes share identical date limits. Adjacent global spectra share
    the exact period axis of their scalogram. Colour bars occupy dedicated
    GridSpec columns and never resize scientific panels.

GLOBAL SPECTRUM POLICY
    Isolated figure:
        X = Fourier period;
        Y = global wavelet power;
        X scale = logarithmic base 2;
        period unit = days or years, selected in wavelet_config.py.

INPUT
    data/processed/<threshold>/quality_control/pwp_centroid_series_qc.csv

OUTPUT
    data/processed/<threshold>/wavelet/
    outputs/tables/<threshold>/wavelet/
    outputs/figures/<threshold>/wavelet/
    outputs/reports/<threshold>/wavelet/

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal

try:
    import pycwt as wavelet
except ImportError as error:
    raise ImportError(
        "Program 11 requires pycwt in the active Python environment."
    ) from error


SCRIPT_FILE = Path(__file__).resolve()
ROOT = SCRIPT_FILE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.metadata import sha256, slug  # noqa: E402
from core.reporting import item, rule, section  # noqa: E402
from core.scientific_variables import (  # noqa: E402
    PWP_VARIABLE_KEYS as VARIABLES,
    VARIABLE_LABELS,
    VARIABLE_UNITS,
)

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
    RUN_ALL_THRESHOLDS,
    SAVE_BBOX,
    SAVE_PAD_INCHES,
    SAVE_TRANSPARENT,
    YEAR_TICK_INTERVAL,
    get_threshold_paths,
    thresholds_to_run,
    validate_project_configuration,
)

from config.wavelet_config import (  # noqa: E402
    ALLOW_SHORT_GAP_INTERPOLATION,
    AR1_MAXIMUM,
    AR1_MINIMUM,
    BAND_MULTIPANEL_SHARE_Y,
    BAND_OVERLAY_NORMALIZATION,
    DAYS_PER_YEAR,
    EXPORT_COMPLEX_COEFFICIENTS,
    EXPORT_IMAGINARY_COEFFICIENTS,
    EXPORT_LOCAL_POWER,
    EXPORT_PDF,
    EXPORT_PHASE,
    EXPORT_PNG,
    EXPORT_REAL_COEFFICIENTS,
    EXPORT_SVG,
    GENERATE_WAVELET_FIGURES,
    GLOBAL_SPECTRUM_LOG2_X_AXIS,
    GLOBAL_SPECTRUM_MAXIMUM_ANNOTATED_PEAKS,
    GLOBAL_SPECTRUM_PERIOD_UNIT,
    GLOBAL_SPECTRUM_SHOW_BACKGROUND,
    GLOBAL_SPECTRUM_SHOW_BANDS,
    GLOBAL_SPECTRUM_SHOW_PEAKS,
    GLOBAL_SPECTRUM_SHOW_SIGNIFICANCE,
    GLOBAL_SPECTRUM_SHADE_BANDS,
    LONGITUDE_UNWRAP_BEFORE_ANALYSIS,
    MAXIMUM_INTERPOLATION_GAP_DAYS,
    MAXIMUM_PERIOD_YEARS,
    MINIMUM_RECORDS,
    MINIMUM_SCALES_PER_BAND,
    MORLET_OMEGA0,
    POWER_FLOOR,
    REAL_PART_LEVELS,
    REAL_PART_PERCENTILE_LIMIT,
    REMOVE_LINEAR_TREND,
    RETAIN_CHECK_QC_RECORDS,
    SAMPLING_INTERVAL_DAYS,
    SCALE_SPACING_DJ,
    SCALOGRAM_PERIOD_UNIT,
    SCALOGRAM_POWER_LEVELS,
    SCALOGRAM_SHOW_BAND_LIMITS,
    SCALOGRAM_SHOW_CONE_OF_INFLUENCE,
    SCALOGRAM_SHOW_SIGNIFICANCE_CONTOUR,
    SHADE_SIGNIFICANT_BAND_POWER,
    SHOW_BAND_MEAN_LINE,
    SHOW_BAND_SIGNIFICANCE,
    SIGNIFICANCE_LEVEL,
    SMALLEST_SCALE_DAYS,
    STANDARDIZE_SERIES,
    WaveletBand,
    enabled_wavelet_bands,
    expected_figure_count_per_variable,
    validate_wavelet_configuration,
    wavelet_configuration_summary_lines,
)


PROGRAM_NAME = "PACIFIC WARM POOL CONTINUOUS WAVELET ANALYSIS"
PROGRAM_VERSION = "4.0.2"





@dataclass(frozen=True)
class Context:
    threshold_c: float
    input_csv: Path
    processed: Path
    tables: Path
    figures: Path
    reports: Path
    metadata: Path
    report: Path


@dataclass(frozen=True)
class Prepared:
    variable: str
    dates: pd.DatetimeIndex
    original: np.ndarray
    completed: np.ndarray
    detrended: np.ndarray
    normalized: np.ndarray
    interpolated: np.ndarray
    mean: float
    standard_deviation: float
    slope_per_day: float
    intercept: float
    ar1: float


@dataclass(frozen=True)
class BandResult:
    band: WaveletBand
    indices: np.ndarray
    power: np.ndarray
    significance: np.ndarray
    significance_ratio: np.ndarray
    inside_coi_fraction: np.ndarray
    mean_power: float
    median_power: float
    maximum_power: float
    maximum_date: pd.Timestamp
    variance: float
    significant_fraction: float
    overall_inside_coi_fraction: float


@dataclass(frozen=True)
class Result:
    variable: str
    coefficients: np.ndarray
    real: np.ndarray
    imaginary: np.ndarray
    phase: np.ndarray
    normalized_power: np.ndarray
    physical_power: np.ndarray
    scales_days: np.ndarray
    periods_days: np.ndarray
    frequencies_per_day: np.ndarray
    coi_days: np.ndarray
    local_significance: np.ndarray
    local_significance_ratio: np.ndarray
    background: np.ndarray
    global_power: np.ndarray
    global_significance: np.ndarray
    bands: tuple[BandResult, ...]












def build_context(threshold: float) -> Context:
    paths = get_threshold_paths(threshold)
    wave = paths.wavelet
    return Context(
        threshold_c=float(paths.threshold_c),
        input_csv=(
            paths.quality_control.processed_dir
            / "pwp_centroid_series_qc.csv"
        ),
        processed=wave.processed_dir,
        tables=wave.table_dir,
        figures=wave.figure_dir,
        reports=wave.report_dir,
        metadata=wave.report_dir / "pwp_wavelet_metadata.json",
        report=wave.report_dir / "pwp_wavelet_analysis_report.txt",
    )


def subdirectories(context: Context) -> dict[str, Path]:
    return {
        "transform": context.processed / "transform",
        "power": context.processed / "power",
        "real": context.processed / "real_part",
        "imaginary": context.processed / "imaginary_part",
        "phase": context.processed / "phase",
        "global": context.processed / "global_spectrum",
        "significance": context.processed / "significance",
        "bands": context.processed / "bands",
        "table_transform": context.tables / "transform_summary",
        "table_global": context.tables / "global_spectrum",
        "table_bands": context.tables / "bands",
        "table_significance": context.tables / "significance",
        "fig_complete": context.figures / "complete",
        "fig_scalogram": context.figures / "scalogram",
        "fig_real": context.figures / "real_part",
        "fig_global": context.figures / "global_spectrum",
        "fig_band_multi": context.figures / "band_multipanel",
        "fig_band_single": context.figures / "band_individual",
        "fig_band_overlay": context.figures / "band_overlay",
    }


def validate_and_create(contexts: tuple[Context, ...]) -> tuple[str, ...]:
    project_messages = validate_project_configuration()
    wavelet_messages = validate_wavelet_configuration()

    if Path(PROJECT_DIR).resolve() != ROOT.resolve():
        raise ValueError(
            f"Project-root mismatch: {PROJECT_DIR} versus {ROOT}"
        )

    for context in contexts:
        for directory in (
            context.processed,
            context.tables,
            context.figures,
            context.reports,
            *subdirectories(context).values(),
        ):
            directory.mkdir(parents=True, exist_ok=True)

    return tuple(project_messages + wavelet_messages)


def load_qc(context: Context) -> tuple[pd.DataFrame, dict[str, int]]:
    if not context.input_csv.is_file():
        raise FileNotFoundError(
            f"Program 11 input not found: {context.input_csv}"
        )

    data = pd.read_csv(context.input_csv, low_memory=False)
    required = ("date", "threshold_c", "qc_status", *VARIABLES)
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    data = data.copy()
    data["date"] = pd.to_datetime(data["date"], errors="raise")
    data["threshold_c"] = pd.to_numeric(
        data["threshold_c"], errors="raise"
    )
    data["qc_status"] = (
        data["qc_status"].astype(str).str.upper().str.strip()
    )
    for variable in VARIABLES:
        data[variable] = pd.to_numeric(data[variable], errors="coerce")

    data = data.sort_values("date").reset_index(drop=True)
    if data["date"].duplicated().any():
        raise ValueError("Duplicate dates in QC input.")

    counts = {
        status: int(np.count_nonzero(data["qc_status"] == status))
        for status in ("PASS", "CHECK", "FAIL")
    }

    retained_status = (
        ("PASS", "CHECK")
        if RETAIN_CHECK_QC_RECORDS
        else ("PASS",)
    )
    data = data.loc[data["qc_status"].isin(retained_status)].copy()

    threshold_values = data["threshold_c"].dropna().unique()
    if threshold_values.size != 1:
        raise ValueError("QC file must contain one threshold.")
    if not np.isclose(threshold_values[0], context.threshold_c):
        raise ValueError("QC threshold differs from active threshold.")

    return data, counts


def unwrap_longitude(values: np.ndarray) -> np.ndarray:
    output = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    output[valid] = np.rad2deg(
        np.unwrap(np.deg2rad(values[valid]))
    )
    return output


def complete_daily(
    dates: pd.DatetimeIndex,
    values: np.ndarray,
) -> tuple[pd.DatetimeIndex, np.ndarray, np.ndarray]:
    full_dates = pd.date_range(dates.min(), dates.max(), freq="D")
    series = pd.Series(values, index=dates).reindex(full_dates)
    originally_missing = series.isna()

    if ALLOW_SHORT_GAP_INTERPOLATION and originally_missing.any():
        series = series.interpolate(
            method="time",
            limit=MAXIMUM_INTERPOLATION_GAP_DAYS,
            limit_area="inside",
        )

    if series.isna().any():
        remaining = series.index[series.isna()]
        raise ValueError(
            "Long or boundary gaps remain after interpolation. "
            f"First: {remaining.min():%Y-%m-%d}; "
            f"last: {remaining.max():%Y-%m-%d}."
        )

    return (
        pd.DatetimeIndex(full_dates),
        series.to_numpy(dtype=float),
        originally_missing.to_numpy(dtype=bool),
    )


def estimate_ar1(values: np.ndarray) -> float:
    value = np.corrcoef(values[:-1], values[1:])[0, 1]
    if not np.isfinite(value):
        value = 0.0
    return float(np.clip(value, AR1_MINIMUM, AR1_MAXIMUM))


def prepare(data: pd.DataFrame, variable: str) -> Prepared:
    values = data[variable].to_numpy(dtype=float)
    if variable == "lon_360" and LONGITUDE_UNWRAP_BEFORE_ANALYSIS:
        values = unwrap_longitude(values)

    dates, completed, interpolated = complete_daily(
        pd.DatetimeIndex(data["date"]),
        values,
    )
    if completed.size < MINIMUM_RECORDS:
        raise ValueError(
            f"{variable}: fewer than {MINIMUM_RECORDS} records."
        )

    time = np.arange(completed.size, dtype=float)
    slope, intercept = np.polyfit(time, completed, 1)
    if REMOVE_LINEAR_TREND:
        detrended = completed - (slope * time + intercept)
    else:
        detrended = completed - np.mean(completed)

    standard_deviation = float(np.std(detrended, ddof=1))
    if standard_deviation <= 0 or not np.isfinite(standard_deviation):
        raise ValueError(f"{variable}: invalid standard deviation.")

    normalized = (
        detrended / standard_deviation
        if STANDARDIZE_SERIES
        else detrended.copy()
    )

    return Prepared(
        variable=variable,
        dates=dates,
        original=values,
        completed=completed,
        detrended=detrended,
        normalized=normalized,
        interpolated=interpolated,
        mean=float(np.mean(completed)),
        standard_deviation=standard_deviation,
        slope_per_day=float(slope),
        intercept=float(intercept),
        ar1=estimate_ar1(normalized),
    )


def scale_count(records: int) -> int:
    maximum_period = min(
        MAXIMUM_PERIOD_YEARS * DAYS_PER_YEAR,
        records * SAMPLING_INTERVAL_DAYS / 3.0,
    )
    return int(
        math.floor(
            math.log2(maximum_period / SMALLEST_SCALE_DAYS)
            / SCALE_SPACING_DJ
        )
    )


def ar1_background(
    periods_days: np.ndarray,
    ar1: float,
    variance: float,
) -> np.ndarray:
    omega = 2 * np.pi * SAMPLING_INTERVAL_DAYS / periods_days
    denominator = 1 + ar1 ** 2 - 2 * ar1 * np.cos(omega)
    return variance * (1 - ar1 ** 2) / np.maximum(
        denominator, POWER_FLOOR
    )


def band_indices(periods: np.ndarray, band: WaveletBand) -> np.ndarray:
    indices = np.flatnonzero(
        (periods >= band.minimum_period_days)
        & (periods <= band.maximum_period_days)
    )
    if indices.size < MINIMUM_SCALES_PER_BAND:
        raise ValueError(
            f"{band.key}: only {indices.size} scales inside the band."
        )
    return indices


def scale_average(
    normalized_power: np.ndarray,
    scales: np.ndarray,
    indices: np.ndarray,
    variance: float,
    mother: Any,
) -> np.ndarray:
    cdelta = float(mother.cdelta)
    if cdelta <= 0 or not np.isfinite(cdelta):
        raise ValueError("Invalid Morlet reconstruction factor.")
    weighted = np.sum(
        normalized_power[indices]
        / scales[indices, np.newaxis],
        axis=0,
    )
    return (
        variance
        * SCALE_SPACING_DJ
        * SAMPLING_INTERVAL_DAYS
        / cdelta
        * weighted
    )


def band_significance(
    band: WaveletBand,
    scales: np.ndarray,
    ar1: float,
    variance: float,
    mother: Any,
) -> float:
    significance, _ = wavelet.significance(
        variance,
        SAMPLING_INTERVAL_DAYS,
        scales,
        2,
        ar1,
        significance_level=SIGNIFICANCE_LEVEL,
        dof=[band.minimum_period_days, band.maximum_period_days],
        wavelet=mother,
    )
    return float(np.nanmean(np.asarray(significance, dtype=float)))


def run_cwt(prepared: Prepared) -> Result:
    mother = wavelet.Morlet(MORLET_OMEGA0)
    coefficients, scales, frequencies, coi, _, _ = wavelet.cwt(
        prepared.normalized,
        SAMPLING_INTERVAL_DAYS,
        SCALE_SPACING_DJ,
        SMALLEST_SCALE_DAYS,
        scale_count(prepared.normalized.size),
        mother,
    )

    coefficients = np.asarray(coefficients, dtype=complex)
    scales = np.asarray(scales, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    coi = np.asarray(coi, dtype=float)
    periods = 1.0 / frequencies

    normalized_power = np.abs(coefficients) ** 2
    variance = prepared.standard_deviation ** 2
    physical_power = normalized_power * variance

    local_norm, _ = wavelet.significance(
        1.0,
        SAMPLING_INTERVAL_DAYS,
        scales,
        0,
        prepared.ar1,
        significance_level=SIGNIFICANCE_LEVEL,
        wavelet=mother,
    )
    local_significance = np.asarray(local_norm, dtype=float) * variance
    local_ratio = physical_power / np.maximum(
        local_significance[:, np.newaxis], POWER_FLOOR
    )

    global_power = np.mean(physical_power, axis=1)
    dof = np.maximum(prepared.normalized.size - scales, 1.0)
    global_significance, _ = wavelet.significance(
        variance,
        SAMPLING_INTERVAL_DAYS,
        scales,
        1,
        prepared.ar1,
        significance_level=SIGNIFICANCE_LEVEL,
        dof=dof,
        wavelet=mother,
    )
    global_significance = np.asarray(global_significance, dtype=float)
    background = ar1_background(periods, prepared.ar1, variance)

    band_results: list[BandResult] = []
    for band in enabled_wavelet_bands():
        indices = band_indices(periods, band)
        power = scale_average(
            normalized_power, scales, indices, variance, mother
        )
        significance_scalar = band_significance(
            band, scales, prepared.ar1, variance, mother
        )
        significance = np.full(power.shape, significance_scalar)
        ratio = power / np.maximum(significance, POWER_FLOOR)
        inside = np.mean(
            periods[indices, np.newaxis] <= coi[np.newaxis, :],
            axis=0,
        )
        maximum_index = int(np.nanargmax(power))
        band_results.append(
            BandResult(
                band=band,
                indices=indices,
                power=power,
                significance=significance,
                significance_ratio=ratio,
                inside_coi_fraction=inside,
                mean_power=float(np.mean(power)),
                median_power=float(np.median(power)),
                maximum_power=float(power[maximum_index]),
                maximum_date=pd.Timestamp(prepared.dates[maximum_index]),
                variance=float(np.var(power, ddof=1)),
                significant_fraction=float(np.mean(ratio >= 1)),
                overall_inside_coi_fraction=float(np.mean(inside)),
            )
        )

    return Result(
        variable=prepared.variable,
        coefficients=coefficients,
        real=np.real(coefficients) * prepared.standard_deviation,
        imaginary=np.imag(coefficients) * prepared.standard_deviation,
        phase=np.angle(coefficients),
        normalized_power=normalized_power,
        physical_power=physical_power,
        scales_days=scales,
        periods_days=periods,
        frequencies_per_day=frequencies,
        coi_days=coi,
        local_significance=local_significance,
        local_significance_ratio=local_ratio,
        background=background,
        global_power=global_power,
        global_significance=global_significance,
        bands=tuple(band_results),
    )


def add_panel(axis: plt.Axes, label: str) -> None:
    axis.text(
        PANEL_LABEL_X,
        PANEL_LABEL_Y,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=PANEL_LABEL_FONT_SIZE,
        fontweight=PANEL_LABEL_FONT_WEIGHT,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.72,
            "pad": 2.0,
        },
        zorder=50,
    )


def date_axis(
    axis: plt.Axes,
    dates: pd.DatetimeIndex,
    labels: bool,
) -> None:
    axis.set_xlim(dates.min(), dates.max())
    axis.margins(x=0)
    axis.xaxis.set_major_locator(
        mdates.YearLocator(base=YEAR_TICK_INTERVAL)
    )
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    if not labels:
        axis.tick_params(axis="x", labelbottom=False)


def period_values(result: Result, unit: str) -> np.ndarray:
    if unit == "years":
        return result.periods_days / DAYS_PER_YEAR
    return result.periods_days


def period_label(unit: str) -> str:
    return (
        "Fourier period (years)"
        if unit == "years"
        else "Fourier period (days)"
    )


def tick_candidates(unit: str) -> np.ndarray:
    days = np.array([
        8, 16, 32, 64, 128, 256,
        DAYS_PER_YEAR,
        2 * DAYS_PER_YEAR,
        4 * DAYS_PER_YEAR,
        8 * DAYS_PER_YEAR,
        16 * DAYS_PER_YEAR,
    ], dtype=float)
    if unit == "years":
        return days / DAYS_PER_YEAR
    return days


def tick_text(value: float, unit: str) -> str:
    if unit == "years":
        if value < 0.75:
            return f"{int(round(value * DAYS_PER_YEAR))} d"
        return f"{value:g} yr"
    if value < 300:
        return f"{int(round(value))} d"
    return f"{value / DAYS_PER_YEAR:g} yr"


def configure_period_y(axis: plt.Axes, result: Result) -> None:
    values = period_values(result, SCALOGRAM_PERIOD_UNIT)
    candidates = tick_candidates(SCALOGRAM_PERIOD_UNIT)
    ticks = candidates[
        (candidates >= values.min())
        & (candidates <= values.max())
    ]
    axis.set_yticks(np.log2(ticks))
    axis.set_yticklabels(
        [tick_text(value, SCALOGRAM_PERIOD_UNIT) for value in ticks]
    )
    axis.set_ylim(np.log2(values.min()), np.log2(values.max()))
    axis.invert_yaxis()
    axis.set_ylabel(period_label(SCALOGRAM_PERIOD_UNIT))


def draw_coi(
    axis: plt.Axes,
    prepared: Prepared,
    result: Result,
) -> None:
    if not SCALOGRAM_SHOW_CONE_OF_INFLUENCE:
        return
    values = period_values(result, SCALOGRAM_PERIOD_UNIT)
    coi = (
        result.coi_days / DAYS_PER_YEAR
        if SCALOGRAM_PERIOD_UNIT == "years"
        else result.coi_days
    )
    coi = np.maximum(coi, values.min())
    axis.plot(prepared.dates, np.log2(coi), linewidth=1.0)
    axis.fill_between(
        prepared.dates,
        np.log2(coi),
        np.log2(values.max()),
        alpha=0.16,
        hatch="//",
    )


def draw_band_limits(axis: plt.Axes) -> None:
    if not SCALOGRAM_SHOW_BAND_LIMITS:
        return
    for band in enabled_wavelet_bands():
        if not band.show_on_scalogram:
            continue
        for days in (
            band.minimum_period_days,
            band.maximum_period_days,
        ):
            value = (
                days / DAYS_PER_YEAR
                if SCALOGRAM_PERIOD_UNIT == "years"
                else days
            )
            axis.axhline(
                np.log2(value),
                linewidth=0.55,
                linestyle="--",
            )


def scalogram(
    axis: plt.Axes,
    prepared: Prepared,
    result: Result,
) -> Any:
    coordinates = np.log2(
        period_values(result, SCALOGRAM_PERIOD_UNIT)
    )
    contour = axis.contourf(
        prepared.dates,
        coordinates,
        np.log2(np.maximum(result.physical_power, POWER_FLOOR)),
        levels=SCALOGRAM_POWER_LEVELS,
        extend="both",
    )
    if SCALOGRAM_SHOW_SIGNIFICANCE_CONTOUR:
        axis.contour(
            prepared.dates,
            coordinates,
            result.local_significance_ratio,
            levels=[1.0],
            linewidths=0.9,
        )
    draw_coi(axis, prepared, result)
    draw_band_limits(axis)
    configure_period_y(axis, result)
    date_axis(axis, prepared.dates, False)
    return contour


def real_scalogram(
    axis: plt.Axes,
    prepared: Prepared,
    result: Result,
) -> Any:
    limit = float(
        np.nanpercentile(
            np.abs(result.real),
            REAL_PART_PERCENTILE_LIMIT,
        )
    )
    if limit <= 0 or not np.isfinite(limit):
        limit = 1.0
    contour = axis.contourf(
        prepared.dates,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        result.real,
        levels=np.linspace(-limit, limit, REAL_PART_LEVELS),
        cmap="RdBu_r",
        extend="both",
    )
    draw_coi(axis, prepared, result)
    draw_band_limits(axis)
    configure_period_y(axis, result)
    date_axis(axis, prepared.dates, False)
    return contour


def normalize_band(values: np.ndarray) -> np.ndarray:
    if BAND_OVERLAY_NORMALIZATION == "absolute":
        return values
    if BAND_OVERLAY_NORMALIZATION == "mean":
        return values / max(float(np.mean(values)), POWER_FLOOR)
    if BAND_OVERLAY_NORMALIZATION == "maximum":
        return values / max(float(np.max(np.abs(values))), POWER_FLOOR)
    standard = float(np.std(values, ddof=1))
    if standard <= POWER_FLOOR:
        return np.zeros_like(values)
    return (values - np.mean(values)) / standard


def band_line(
    axis: plt.Axes,
    dates: pd.DatetimeIndex,
    band: BandResult,
    normalized: bool,
) -> None:
    values = normalize_band(band.power) if normalized else band.power
    axis.plot(
        dates,
        values,
        linewidth=1.0,
        label=band.band.short_label,
    )
    if SHOW_BAND_SIGNIFICANCE and not normalized:
        axis.plot(
            dates,
            band.significance,
            linestyle="--",
            linewidth=0.75,
            label="95% significance",
        )
        if SHADE_SIGNIFICANT_BAND_POWER:
            axis.fill_between(
                dates,
                0,
                values,
                where=band.significance_ratio >= 1,
                alpha=0.18,
            )
    if SHOW_BAND_MEAN_LINE:
        axis.axhline(np.mean(values), linestyle=":", linewidth=0.7)
    axis.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )


def global_periods(result: Result) -> np.ndarray:
    if GLOBAL_SPECTRUM_PERIOD_UNIT == "years":
        return result.periods_days / DAYS_PER_YEAR
    return result.periods_days


def configure_global_x(axis: plt.Axes, result: Result) -> None:
    periods = global_periods(result)
    axis.set_xlim(periods.min(), periods.max())
    if GLOBAL_SPECTRUM_LOG2_X_AXIS:
        axis.set_xscale("log", base=2)
    candidates = tick_candidates(GLOBAL_SPECTRUM_PERIOD_UNIT)
    ticks = candidates[
        (candidates >= periods.min())
        & (candidates <= periods.max())
    ]
    axis.set_xticks(ticks)
    axis.set_xticklabels(
        [tick_text(value, GLOBAL_SPECTRUM_PERIOD_UNIT) for value in ticks]
    )
    axis.set_xlabel(period_label(GLOBAL_SPECTRUM_PERIOD_UNIT))


def draw_global_bands(axis: plt.Axes) -> None:
    if not GLOBAL_SPECTRUM_SHOW_BANDS:
        return
    for band in enabled_wavelet_bands():
        if not band.show_on_global_spectrum:
            continue
        minimum = (
            band.minimum_period_days / DAYS_PER_YEAR
            if GLOBAL_SPECTRUM_PERIOD_UNIT == "years"
            else band.minimum_period_days
        )
        maximum = (
            band.maximum_period_days / DAYS_PER_YEAR
            if GLOBAL_SPECTRUM_PERIOD_UNIT == "years"
            else band.maximum_period_days
        )
        if GLOBAL_SPECTRUM_SHADE_BANDS:
            axis.axvspan(minimum, maximum, alpha=0.08)
        axis.axvline(minimum, linestyle="--", linewidth=0.7)
        axis.axvline(maximum, linestyle="--", linewidth=0.7)
        axis.text(
            math.sqrt(minimum * maximum),
            0.98,
            band.short_label,
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="top",
            rotation=90,
            fontsize=8,
        )


def global_peaks(result: Result) -> pd.DataFrame:
    prominence = 0.02 * float(np.nanmax(result.global_power))
    indices, properties = signal.find_peaks(
        result.global_power,
        prominence=prominence,
    )
    periods = global_periods(result)
    return pd.DataFrame({
        "period": periods[indices],
        "period_unit": GLOBAL_SPECTRUM_PERIOD_UNIT,
        "power": result.global_power[indices],
        "significance": result.global_significance[indices],
        "significant": (
            result.global_power[indices]
            >= result.global_significance[indices]
        ),
        "prominence": properties.get("prominences", np.array([])),
    }).sort_values(
        ["significant", "prominence", "power"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def global_spectrum(
    axis: plt.Axes,
    result: Result,
    annotate: bool,
) -> None:
    periods = global_periods(result)
    axis.plot(
        periods,
        result.global_power,
        linewidth=1.3,
        label="Global wavelet power",
    )
    if GLOBAL_SPECTRUM_SHOW_BACKGROUND:
        axis.plot(
            periods,
            result.background,
            linestyle=":",
            linewidth=1.0,
            label="AR(1) background",
        )
    if GLOBAL_SPECTRUM_SHOW_SIGNIFICANCE:
        axis.plot(
            periods,
            result.global_significance,
            linestyle="--",
            linewidth=1.0,
            label="95% significance",
        )
    draw_global_bands(axis)
    configure_global_x(axis, result)
    axis.set_ylabel("Global wavelet power")
    axis.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    axis.legend(frameon=False, loc="best")

    if annotate and GLOBAL_SPECTRUM_SHOW_PEAKS:
        for record in global_peaks(result).head(
            GLOBAL_SPECTRUM_MAXIMUM_ANNOTATED_PEAKS
        ).to_dict(orient="records"):
            axis.annotate(
                tick_text(
                    float(record["period"]),
                    GLOBAL_SPECTRUM_PERIOD_UNIT,
                ),
                xy=(record["period"], record["power"]),
                xytext=(4, 6),
                textcoords="offset points",
                fontsize=8,
            )


def save(figure: plt.Figure, base: Path) -> tuple[Path, ...]:
    files: list[Path] = []
    for enabled, suffix, dpi in (
        (EXPORT_PNG, ".png", FIGURE_DPI),
        (EXPORT_PDF, ".pdf", None),
        (EXPORT_SVG, ".svg", None),
    ):
        if not enabled:
            continue
        path = base.with_suffix(suffix)
        kwargs = {
            "bbox_inches": SAVE_BBOX,
            "pad_inches": SAVE_PAD_INCHES,
            "transparent": SAVE_TRANSPARENT,
        }
        if dpi is not None:
            kwargs["dpi"] = dpi
        figure.savefig(path, **kwargs)
        files.append(path)
    plt.close(figure)
    return tuple(files)


def traditional(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure = plt.figure(figsize=(17, 13.8))
    grid = figure.add_gridspec(
        4, 3,
        width_ratios=(1.0, 0.028, 0.24),
        height_ratios=(0.72, 1.45, 1.45, 1.0),
        left=0.075, right=0.975,
        bottom=0.075, top=0.945,
        hspace=0.16, wspace=0.045,
    )
    a = figure.add_subplot(grid[0, 0])
    b = figure.add_subplot(grid[1, 0], sharex=a)
    c = figure.add_subplot(grid[2, 0], sharex=a, sharey=b)
    d = figure.add_subplot(grid[3, 0], sharex=a)
    cb1 = figure.add_subplot(grid[1, 1])
    cb2 = figure.add_subplot(grid[2, 1])
    e = figure.add_subplot(grid[1:3, 2], sharey=b)

    a.plot(prepared.dates, prepared.detrended, linewidth=0.75)
    a.axhline(0, linewidth=0.7)
    a.set_ylabel(
        f"{VARIABLE_LABELS[result.variable]}\n"
        f"({VARIABLE_UNITS[result.variable]})"
    )
    a.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    date_axis(a, prepared.dates, False)
    add_panel(a, "A)")

    power = scalogram(b, prepared, result)
    add_panel(b, "B)")
    real = real_scalogram(c, prepared, result)
    add_panel(c, "C)")

    for band in result.bands:
        band_line(d, prepared.dates, band, True)
    d.set_ylabel(f"Band power\n({BAND_OVERLAY_NORMALIZATION})")
    d.legend(frameon=False, ncol=min(3, len(result.bands)))
    date_axis(d, prepared.dates, True)
    add_panel(d, "D)")

    figure.colorbar(power, cax=cb1).set_label("log₂ wavelet power")
    figure.colorbar(real, cax=cb2).set_label("Real coefficient")

    e.plot(
        result.global_power,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        label="Global power",
    )
    e.plot(
        result.global_significance,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        linestyle="--",
        label="95% significance",
    )
    e.plot(
        result.background,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        linestyle=":",
        label="AR(1) background",
    )
    e.tick_params(axis="y", labelleft=False)
    e.set_xlabel("Power")
    e.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    e.legend(frameon=False, loc="lower right")
    add_panel(e, "E)")

    figure.suptitle(
        f"Complex Morlet CWT — {VARIABLE_LABELS[result.variable]} "
        f"(SST ≥ {context.threshold_c:g} °C)",
        y=0.982,
    )
    return save(
        figure,
        directories["fig_complete"]
        / f"pwp_wavelet_complete_{result.variable}",
    )


def band_multipanel(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    count = len(result.bands)
    figure, axes = plt.subplots(
        count, 1,
        figsize=(15, max(3 * count, 5)),
        sharex=True,
        sharey=BAND_MULTIPANEL_SHARE_Y,
        squeeze=False,
    )
    for index, (axis, band) in enumerate(zip(axes[:, 0], result.bands)):
        band_line(axis, prepared.dates, band, False)
        axis.set_ylabel("Scale-averaged\npower")
        axis.set_title(
            f"{band.band.label}: {band.band.period_range_text}",
            loc="left",
            fontsize=10,
        )
        date_axis(axis, prepared.dates, index == count - 1)
        add_panel(axis, f"{chr(65 + index)})")
        if index == 0:
            axis.legend(frameon=False, loc="upper right")
    figure.suptitle(
        f"Scale-averaged wavelet bands — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.995,
    )
    return save(
        figure,
        directories["fig_band_multi"]
        / f"pwp_wavelet_bands_multipanel_{result.variable}",
    )


def individual_bands(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    files: list[Path] = []
    for index, band in enumerate(result.bands, start=1):
        figure, axis = plt.subplots(figsize=(15, 4.8))
        band_line(axis, prepared.dates, band, False)
        axis.set_xlabel("Year")
        axis.set_ylabel("Scale-averaged wavelet power")
        axis.legend(frameon=False, loc="upper right")
        date_axis(axis, prepared.dates, True)
        figure.suptitle(
            f"{band.band.label} ({band.band.period_range_text}) — "
            f"{VARIABLE_LABELS[result.variable]} — "
            f"SST ≥ {context.threshold_c:g} °C",
            y=0.985,
        )
        files.extend(
            save(
                figure,
                directories["fig_band_single"]
                / (
                    f"pwp_wavelet_band_{index:02d}_"
                    f"{slug(band.band.key)}_{result.variable}"
                ),
            )
        )
    return tuple(files)


def band_overlay(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure, axis = plt.subplots(figsize=(15, 5.8))
    for band in result.bands:
        band_line(axis, prepared.dates, band, True)
    axis.set_xlabel("Year")
    axis.set_ylabel(f"Band power ({BAND_OVERLAY_NORMALIZATION})")
    axis.legend(frameon=False, ncol=min(3, len(result.bands)))
    date_axis(axis, prepared.dates, True)
    figure.suptitle(
        f"User-defined wavelet bands — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.985,
    )
    return save(
        figure,
        directories["fig_band_overlay"]
        / f"pwp_wavelet_all_bands_overlay_{result.variable}",
    )


def integrated(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure = plt.figure(figsize=(17, 11.8))
    grid = figure.add_gridspec(
        3, 3,
        width_ratios=(1.0, 0.028, 0.24),
        height_ratios=(0.72, 1.6, 1.0),
        left=0.075, right=0.975,
        bottom=0.075, top=0.94,
        hspace=0.16, wspace=0.045,
    )
    a = figure.add_subplot(grid[0, 0])
    b = figure.add_subplot(grid[1, 0], sharex=a)
    c = figure.add_subplot(grid[2, 0], sharex=a)
    cb = figure.add_subplot(grid[1, 1])
    d = figure.add_subplot(grid[1, 2], sharey=b)

    a.plot(prepared.dates, prepared.detrended, linewidth=0.75)
    a.axhline(0, linewidth=0.7)
    a.set_ylabel(
        f"{VARIABLE_LABELS[result.variable]}\n"
        f"({VARIABLE_UNITS[result.variable]})"
    )
    a.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    date_axis(a, prepared.dates, False)
    add_panel(a, "A)")

    contour = scalogram(b, prepared, result)
    add_panel(b, "B)")
    figure.colorbar(contour, cax=cb).set_label("log₂ wavelet power")

    for band in result.bands:
        band_line(c, prepared.dates, band, True)
    c.set_ylabel(f"Band power\n({BAND_OVERLAY_NORMALIZATION})")
    c.legend(frameon=False, ncol=min(3, len(result.bands)))
    date_axis(c, prepared.dates, True)
    add_panel(c, "C)")

    d.plot(
        result.global_power,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        label="Global power",
    )
    d.plot(
        result.global_significance,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        linestyle="--",
        label="95% significance",
    )
    d.plot(
        result.background,
        np.log2(period_values(result, SCALOGRAM_PERIOD_UNIT)),
        linestyle=":",
        label="AR(1) background",
    )
    d.tick_params(axis="y", labelleft=False)
    d.set_xlabel("Power")
    d.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    d.legend(frameon=False, loc="lower right")
    add_panel(d, "D)")

    figure.suptitle(
        f"Wavelet decomposition and user bands — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.98,
    )
    return save(
        figure,
        directories["fig_complete"]
        / f"pwp_wavelet_series_scalogram_bands_global_{result.variable}",
    )


def series_scalogram(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure = plt.figure(figsize=(15.5, 8.8))
    grid = figure.add_gridspec(
        2, 2,
        width_ratios=(1.0, 0.028),
        height_ratios=(0.72, 1.75),
        left=0.08, right=0.96,
        bottom=0.09, top=0.93,
        hspace=0.12, wspace=0.04,
    )
    a = figure.add_subplot(grid[0, 0])
    b = figure.add_subplot(grid[1, 0], sharex=a)
    cb = figure.add_subplot(grid[1, 1])

    a.plot(prepared.dates, prepared.detrended, linewidth=0.75)
    a.axhline(0, linewidth=0.7)
    a.set_ylabel(
        f"{VARIABLE_LABELS[result.variable]}\n"
        f"({VARIABLE_UNITS[result.variable]})"
    )
    a.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )
    date_axis(a, prepared.dates, False)
    add_panel(a, "A)")

    contour = scalogram(b, prepared, result)
    date_axis(b, prepared.dates, True)
    b.set_xlabel("Year")
    add_panel(b, "B)")
    figure.colorbar(contour, cax=cb).set_label("log₂ wavelet power")

    figure.suptitle(
        f"Continuous wavelet power — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.98,
    )
    return save(
        figure,
        directories["fig_scalogram"]
        / f"pwp_wavelet_series_scalogram_{result.variable}",
    )


def isolated_global(
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure, axis = plt.subplots(figsize=(12.5, 6.8))
    global_spectrum(axis, result, True)
    figure.suptitle(
        f"Global Wavelet Spectrum — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.985,
    )
    return save(
        figure,
        directories["fig_global"]
        / f"pwp_wavelet_global_spectrum_{result.variable}",
    )


def isolated_real(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    figure = plt.figure(figsize=(15.5, 8))
    grid = figure.add_gridspec(
        1, 2,
        width_ratios=(1.0, 0.028),
        left=0.08, right=0.96,
        bottom=0.10, top=0.91,
        wspace=0.04,
    )
    axis = figure.add_subplot(grid[0, 0])
    cb = figure.add_subplot(grid[0, 1])
    contour = real_scalogram(axis, prepared, result)
    date_axis(axis, prepared.dates, True)
    axis.set_xlabel("Year")
    figure.colorbar(contour, cax=cb).set_label("Real coefficient")
    figure.suptitle(
        f"Real component of Complex Morlet CWT — "
        f"{VARIABLE_LABELS[result.variable]} — "
        f"SST ≥ {context.threshold_c:g} °C",
        y=0.97,
    )
    return save(
        figure,
        directories["fig_real"]
        / f"pwp_wavelet_real_part_{result.variable}",
    )


def figures(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[Path, ...]:
    files: list[Path] = []
    mapping = (
        ("traditional_complete", traditional),
        ("bands_multipanel", band_multipanel),
        ("individual_band_figures", individual_bands),
        ("all_bands_overlay", band_overlay),
        ("series_scalogram_bands_global", integrated),
        ("series_and_scalogram", series_scalogram),
        ("global_spectrum_only", isolated_global),
        ("real_part_only", isolated_real),
    )

    for key, function in mapping:
        if not GENERATE_WAVELET_FIGURES[key]:
            continue
        if key == "global_spectrum_only":
            files.extend(function(result, context, directories))
        else:
            files.extend(function(prepared, result, context, directories))

    formats = int(EXPORT_PNG) + int(EXPORT_PDF) + int(EXPORT_SVG)
    expected = expected_figure_count_per_variable() * formats
    if len(files) != expected:
        raise RuntimeError(
            f"Figure-count mismatch: expected {expected}, created {len(files)}."
        )
    return tuple(files)


def export(
    prepared: Prepared,
    result: Result,
    context: Context,
    directories: dict[str, Path],
) -> tuple[tuple[Path, ...], pd.DataFrame, pd.DataFrame]:
    files: list[Path] = []

    arrays: dict[str, np.ndarray] = {
        "dates": prepared.dates.to_numpy(dtype="datetime64[ns]"),
        "scales_days": result.scales_days,
        "fourier_period_days": result.periods_days,
        "frequencies_per_day": result.frequencies_per_day,
        "coi_days": result.coi_days,
        "interpolated_flag": prepared.interpolated,
        "local_significance": result.local_significance,
        "local_significance_ratio": result.local_significance_ratio,
    }
    if EXPORT_COMPLEX_COEFFICIENTS:
        arrays["coefficients_complex"] = result.coefficients
    if EXPORT_REAL_COEFFICIENTS:
        arrays["coefficients_real"] = result.real
    if EXPORT_IMAGINARY_COEFFICIENTS:
        arrays["coefficients_imaginary"] = result.imaginary
    if EXPORT_PHASE:
        arrays["phase_radians"] = result.phase
    if EXPORT_LOCAL_POWER:
        arrays["normalized_power"] = result.normalized_power
        arrays["physical_power"] = result.physical_power

    transform_file = (
        directories["transform"]
        / f"pwp_wavelet_transform_{result.variable}.npz"
    )
    np.savez_compressed(transform_file, **arrays)
    files.append(transform_file)

    global_table = pd.DataFrame({
        "scale_days": result.scales_days,
        "fourier_period_days": result.periods_days,
        "fourier_period_years": result.periods_days / DAYS_PER_YEAR,
        "global_wavelet_power": result.global_power,
        "ar1_background_power": result.background,
        "global_significance_power": result.global_significance,
        "significance_ratio": (
            result.global_power
            / np.maximum(result.global_significance, POWER_FLOOR)
        ),
    })
    global_file = (
        directories["table_global"]
        / f"pwp_wavelet_global_spectrum_{result.variable}.csv"
    )
    global_table.to_csv(
        global_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.append(global_file)

    peaks_file = (
        directories["table_global"]
        / f"pwp_wavelet_global_peaks_{result.variable}.csv"
    )
    global_peaks(result).to_csv(
        peaks_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.append(peaks_file)

    consolidated = pd.DataFrame({
        "date": prepared.dates,
        "threshold_c": context.threshold_c,
        "variable": result.variable,
    })
    summaries: list[dict[str, Any]] = []
    for band in result.bands:
        table = pd.DataFrame({
            "date": prepared.dates,
            "threshold_c": context.threshold_c,
            "variable": result.variable,
            "band_key": band.band.key,
            "band_label": band.band.label,
            "minimum_period_days": band.band.minimum_period_days,
            "maximum_period_days": band.band.maximum_period_days,
            "scale_averaged_power": band.power,
            "scale_averaged_significance": band.significance,
            "significance_ratio": band.significance_ratio,
            "significant_95pct": band.significance_ratio >= 1,
            "inside_coi_scale_fraction": band.inside_coi_fraction,
        })
        band_file = (
            directories["table_bands"]
            / f"pwp_wavelet_band_{slug(band.band.key)}_{result.variable}.csv"
        )
        table.to_csv(
            band_file,
            index=False,
            float_format="%.12e",
            lineterminator="\n",
        )
        files.append(band_file)
        consolidated[f"{band.band.key}_power"] = band.power
        consolidated[
            f"{band.band.key}_significance_ratio"
        ] = band.significance_ratio
        summaries.append({
            "threshold_c": context.threshold_c,
            "variable": result.variable,
            "band_key": band.band.key,
            "band_label": band.band.label,
            "minimum_period_days": band.band.minimum_period_days,
            "maximum_period_days": band.band.maximum_period_days,
            "scale_count": int(band.indices.size),
            "mean_power": band.mean_power,
            "median_power": band.median_power,
            "maximum_power": band.maximum_power,
            "maximum_power_date": band.maximum_date,
            "variance": band.variance,
            "significant_fraction": band.significant_fraction,
            "overall_inside_coi_fraction": (
                band.overall_inside_coi_fraction
            ),
        })

    consolidated_file = (
        directories["table_bands"]
        / f"pwp_wavelet_all_bands_{result.variable}.csv"
    )
    consolidated.to_csv(
        consolidated_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.append(consolidated_file)

    band_summary = pd.DataFrame(summaries)
    band_summary_file = (
        directories["table_bands"]
        / f"pwp_wavelet_band_summary_{result.variable}.csv"
    )
    band_summary.to_csv(
        band_summary_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.append(band_summary_file)

    transform_summary = pd.DataFrame({
        "threshold_c": [context.threshold_c],
        "variable": [result.variable],
        "records": [prepared.dates.size],
        "start_date": [prepared.dates.min()],
        "end_date": [prepared.dates.max()],
        "interpolated_records": [
            int(np.count_nonzero(prepared.interpolated))
        ],
        "series_mean": [prepared.mean],
        "series_standard_deviation": [prepared.standard_deviation],
        "linear_slope_per_day": [prepared.slope_per_day],
        "linear_slope_per_year": [
            prepared.slope_per_day * DAYS_PER_YEAR
        ],
        "ar1_coefficient": [prepared.ar1],
        "scale_count": [result.scales_days.size],
        "minimum_fourier_period_days": [result.periods_days.min()],
        "maximum_fourier_period_days": [result.periods_days.max()],
        "enabled_band_count": [len(result.bands)],
        "expected_figures_per_variable": [
            expected_figure_count_per_variable()
        ],
    })
    transform_summary_file = (
        directories["table_transform"]
        / f"pwp_wavelet_transform_summary_{result.variable}.csv"
    )
    transform_summary.to_csv(
        transform_summary_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.append(transform_summary_file)

    return tuple(files), transform_summary, band_summary


def write_outputs(
    context: Context,
    checksum: str,
    counts: dict[str, int],
    transform_summary: pd.DataFrame,
    band_summary: pd.DataFrame,
    files: list[Path],
) -> None:
    metadata = {
        "program": PROGRAM_NAME,
        "version": PROGRAM_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "threshold_c": context.threshold_c,
        "input_file": str(context.input_csv),
        "input_sha256": checksum,
        "qc_counts": counts,
        "wavelet": {
            "mother": "Complex Morlet",
            "omega0": MORLET_OMEGA0,
            "sampling_interval_days": SAMPLING_INTERVAL_DAYS,
            "scale_spacing_dj": SCALE_SPACING_DJ,
            "maximum_period_years": MAXIMUM_PERIOD_YEARS,
            "significance_level": SIGNIFICANCE_LEVEL,
            "background": "AR(1) red noise",
        },
        "bands": [
            {
                "key": band.key,
                "label": band.label,
                "minimum_period_days": band.minimum_period_days,
                "maximum_period_days": band.maximum_period_days,
            }
            for band in enabled_wavelet_bands()
        ],
        "global_spectrum_orientation": {
            "x": f"Fourier period ({GLOBAL_SPECTRUM_PERIOD_UNIT})",
            "y": "Global wavelet power",
            "log2_x": GLOBAL_SPECTRUM_LOG2_X_AXIS,
        },
        "figure_families": GENERATE_WAVELET_FIGURES,
        "transform_summary": transform_summary.to_dict(orient="records"),
        "band_summary": band_summary.to_dict(orient="records"),
        "created_files": [str(path) for path in files],
    }
    context.metadata.write_text(
        json.dumps(
            metadata,
            indent=2,
            ensure_ascii=False,
            default=str,
        ) + "\n",
        encoding="utf-8",
    )
    files.append(context.metadata)

    lines = [
        "PACIFIC WARM POOL CONTINUOUS WAVELET ANALYSIS REPORT",
        "=" * 78,
        "",
        f"Program                    : {PROGRAM_NAME}",
        f"Version                    : {PROGRAM_VERSION}",
        f"Threshold                  : {context.threshold_c:.1f} °C",
        f"Input file                 : {context.input_csv}",
        f"Input SHA-256              : {checksum}",
        f"PASS records               : {counts['PASS']:,}",
        f"CHECK records              : {counts['CHECK']:,}",
        f"FAIL records               : {counts['FAIL']:,}",
        "",
        "SCIENTIFIC CONFIGURATION",
        "-" * 78,
        *wavelet_configuration_summary_lines(),
        "",
        "INTERPRETIVE CAUTIONS",
        "-" * 78,
        (
            "Interpret local power jointly with the cone of influence and "
            "the AR(1) significance contour."
        ),
        (
            "Pointwise significance does not correct automatically for "
            "multiple time-scale testing."
        ),
        (
            "Overlapping user-defined bands are not independent summaries."
        ),
        (
            "Interdecadal results contain few independent cycles in the "
            "satellite-era record."
        ),
        "",
        "REFERENCES",
        "-" * 78,
        "Torrence and Compo (1998). A Practical Guide to Wavelet Analysis.",
        "Domingues et al. (2005). Wavelet technique in atmospheric sciences.",
        "Gu and Philander (1995). Annual and interannual tropical variability.",
        "Torrence and Webster (1999). Interdecadal ENSO-monsoon changes.",
        "",
        "FILES",
        "-" * 78,
        *[str(path) for path in files],
        "",
        "STATUS                     : SUCCESS",
        "=" * 78,
    ]
    context.report.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    files.append(context.report)


def process_context(context: Context) -> tuple[Path, ...]:
    section(f"PROCESSING THRESHOLD — {context.threshold_c:.1f} °C")
    data, counts = load_qc(context)
    checksum = sha256(context.input_csv)
    directories = subdirectories(context)

    item("QC input file", context.input_csv)
    item("Records", f"{len(data):,}")
    item("Start date", f"{data['date'].min():%Y-%m-%d}")
    item("End date", f"{data['date'].max():%Y-%m-%d}")

    transform_summaries: list[pd.DataFrame] = []
    band_summaries: list[pd.DataFrame] = []
    files: list[Path] = []

    for variable in VARIABLES:
        section(
            f"WAVELET TRANSFORM — {context.threshold_c:.1f} °C — {variable}"
        )
        prepared = prepare(data, variable)
        item(
            "Interpolated records",
            f"{np.count_nonzero(prepared.interpolated):,}",
        )
        item(
            "Standard deviation",
            f"{prepared.standard_deviation:.10f}",
        )
        item("AR(1) coefficient", f"{prepared.ar1:.10f}")

        result = run_cwt(prepared)
        export_files, transform_summary, band_summary = export(
            prepared, result, context, directories
        )
        figure_files = figures(
            prepared, result, context, directories
        )
        files.extend(export_files)
        files.extend(figure_files)
        transform_summaries.append(transform_summary)
        band_summaries.append(band_summary)

    transform_all = pd.concat(transform_summaries, ignore_index=True)
    band_all = pd.concat(band_summaries, ignore_index=True)

    transform_file = (
        directories["table_transform"]
        / "pwp_wavelet_transform_summary.csv"
    )
    band_file = (
        directories["table_bands"]
        / "pwp_wavelet_band_summary.csv"
    )
    transform_all.to_csv(
        transform_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    band_all.to_csv(
        band_file,
        index=False,
        float_format="%.12e",
        lineterminator="\n",
    )
    files.extend((transform_file, band_file))

    write_outputs(
        context,
        checksum,
        counts,
        transform_all,
        band_all,
        files,
    )
    return tuple(files)


def main() -> None:
    rule()
    print(PROGRAM_NAME)
    rule()

    contexts = tuple(
        build_context(threshold)
        for threshold in thresholds_to_run()
    )
    messages = validate_and_create(contexts)

    section("CONFIGURATION")
    item("Program version", PROGRAM_VERSION)
    item("Project root", PROJECT_DIR)
    item("Script", SCRIPT_FILE)
    item(
        "Execution mode",
        "all thresholds" if RUN_ALL_THRESHOLDS else "single threshold",
    )
    item(
        "Thresholds selected",
        ", ".join(f"{context.threshold_c:.1f} °C" for context in contexts),
    )
    for line in wavelet_configuration_summary_lines():
        print(line)

    section("VALIDATION")
    for message in messages:
        print(message)

    files: list[Path] = []
    for context in contexts:
        files.extend(process_context(context))

    section("FILES CREATED")
    for path in files:
        print(path)

    print()
    rule()
    print("PROGRAM 11 COMPLETED SUCCESSFULLY.")
    print(
        "Band limits and figure families were controlled by "
        "config/wavelet_config.py."
    )
    rule()


if __name__ == "__main__":
    main()
