# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

FILE
    config/config.py

VERSION
    4.1.0

PURPOSE
    Provide the authoritative, portable, typed, and validated configuration
    layer for the complete Pacific Warm Pool scientific-analysis framework.

ARCHITECTURAL ROLE
    This file is the single source of truth for:

        - project discovery;
        - project metadata and software versioning;
        - input and output roots;
        - threshold selection;
        - threshold-safe naming;
        - scientific-module registration;
        - threshold-specific directory construction;
        - cross-threshold directory construction;
        - project-level methodological and synthesis products;
        - raw and processed external climate-index locations;
        - Program 21 climate-index ingestion architecture;
        - dedicated PWP × Niño/SOI cross-wavelet paths;
        - spectral-synthesis paths;
        - spherical-Earth and OISST-grid constants;
        - publication-figure standards;
        - time-axis standards;
        - reusable validation and directory helpers.

PORTABILITY
    PROJECT_DIR is derived from the physical location of this file:

        <project_root>/config/config.py

    Therefore, the project can be moved, copied, cloned, or restored without
    editing machine-specific absolute paths.

SCIENTIFIC THRESHOLDS
    28.0 °C -> folder "28"
    28.5 °C -> folder "28.5"
    29.0 °C -> folder "29"

DIRECTORY PRINCIPLE
    Every threshold-dependent scientific module has a mirrored structure:

        data/processed/<threshold>/<module>/
        outputs/figures/<threshold>/<module>/
        outputs/tables/<threshold>/<module>/
        outputs/reports/<threshold>/<module>/
        outputs/logs/<threshold>/<module>/

    Every cross-threshold module has a mirrored structure:

        outputs/figures/threshold_comparison/<module>/
        outputs/tables/threshold_comparison/<module>/
        outputs/reports/threshold_comparison/<module>/

DESIGN PRINCIPLES
    1. No analysis program constructs threshold-dependent paths manually.
    2. No machine-specific project path is hard-coded.
    3. Module names are registered once and validated centrally.
    4. Directory creation is explicit, deterministic, and idempotent.
    5. Backward-compatible aliases are retained only where scientifically safe.
    6. Configuration validation does not modify scientific data.
    7. Structural errors fail loudly with actionable messages.

PYTHON
    Python 3.10+

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Final, Iterable, Iterator, Mapping, Sequence


# =============================================================================
# PROJECT IDENTITY AND PORTABLE DISCOVERY
# =============================================================================

CONFIG_FILE: Final[Path] = Path(__file__).resolve()
CONFIG_DIR: Final[Path] = CONFIG_FILE.parent
PROJECT_DIR: Final[Path] = CONFIG_DIR.parent
PROJECT_ROOT: Final[Path] = PROJECT_DIR

PROJECT_NAME: Final[str] = (
    "Pacific Warm Pool Scientific Analysis Pipeline"
)
PROJECT_SHORT_NAME: Final[str] = "PWP Scientific Pipeline"
PROJECT_VERSION: Final[str] = "4.1.0"
PROJECT_AUTHOR: Final[str] = "Fabio Vieira Machado"

SRC_DIR: Final[Path] = PROJECT_DIR / "src"
DOCUMENTATION_DIR: Final[Path] = PROJECT_DIR / "documentation"
DOCUMENTATION_VALIDATION_DIR: Final[Path] = (
    DOCUMENTATION_DIR / "validation"
)
DOCUMENTATION_PROGRAM21_DIR: Final[Path] = (
    DOCUMENTATION_DIR / "program_21"
)
DOCUMENTATION_PROGRAM19_NINO_SOI_DIR: Final[Path] = (
    DOCUMENTATION_DIR / "program_19_nino_soi"
)
TESTS_DIR: Final[Path] = PROJECT_DIR / "tests"


# =============================================================================
# PRIMARY DIRECTORY ROOTS
# =============================================================================

DATA_DIR: Final[Path] = PROJECT_DIR / "data"
MASK_DIR: Final[Path] = DATA_DIR / "masks"

RAW_DATA_DIR: Final[Path] = DATA_DIR / "raw"
RAW_DIR: Final[Path] = RAW_DATA_DIR
RAW_OISST_DIR: Final[Path] = RAW_DATA_DIR

PROCESSED_DIR: Final[Path] = DATA_DIR / "processed"
PREPARED_DIR: Final[Path] = DATA_DIR / "prepared"

# Canonical OSAF climate-index architecture (Program 21).
RAW_CLIMATE_INDICES_DIR: Final[Path] = (
    RAW_DATA_DIR / "climate_indices"
)
RAW_NINO_DIR: Final[Path] = (
    RAW_CLIMATE_INDICES_DIR / "nino"
)
RAW_SOI_DIR: Final[Path] = (
    RAW_CLIMATE_INDICES_DIR / "soi"
)

PROCESSED_CLIMATE_INDICES_DIR: Final[Path] = (
    PROCESSED_DIR / "climate_indices"
)
PROCESSED_NINO_DIR: Final[Path] = (
    PROCESSED_CLIMATE_INDICES_DIR / "nino"
)
PROCESSED_SOI_DIR: Final[Path] = (
    PROCESSED_CLIMATE_INDICES_DIR / "soi"
)

# Retained only for compatibility with older archived scripts.
LEGACY_EXTERNAL_INDEX_DIR: Final[Path] = (
    DATA_DIR / "external_indices"
)
EXTERNAL_INDEX_DIR: Final[Path] = LEGACY_EXTERNAL_INDEX_DIR

OUTPUT_DIR: Final[Path] = PROJECT_DIR / "outputs"
FIGURE_DIR: Final[Path] = OUTPUT_DIR / "figures"
TABLE_DIR: Final[Path] = OUTPUT_DIR / "tables"
REPORT_DIR: Final[Path] = OUTPUT_DIR / "reports"
LOG_DIR: Final[Path] = OUTPUT_DIR / "logs"
MODEL_DIR: Final[Path] = OUTPUT_DIR / "models"
FORECAST_DIR: Final[Path] = OUTPUT_DIR / "forecasts"

# Program-21 project-level outputs.
CLIMATE_INDEX_TABLE_DIR: Final[Path] = (
    TABLE_DIR / "climate_indices"
)
CLIMATE_INDEX_REPORT_DIR: Final[Path] = (
    REPORT_DIR / "climate_indices"
)

NINO_INGESTION_SUMMARY_FILE: Final[Path] = (
    CLIMATE_INDEX_TABLE_DIR / "nino_ingestion_summary.csv"
)
SOI_INGESTION_SUMMARY_FILE: Final[Path] = (
    CLIMATE_INDEX_TABLE_DIR / "soi_ingestion_summary.csv"
)

NINO_VALIDATION_REPORT_FILE: Final[Path] = (
    CLIMATE_INDEX_REPORT_DIR / "nino_validation.txt"
)
SOI_VALIDATION_REPORT_FILE: Final[Path] = (
    CLIMATE_INDEX_REPORT_DIR / "soi_validation.txt"
)
CLIMATE_INDICES_INGESTION_REPORT_FILE: Final[Path] = (
    CLIMATE_INDEX_REPORT_DIR
    / "climate_indices_ingestion_report.txt"
)

ARCHIVE_DIR: Final[Path] = PROJECT_DIR / "archive"
TEMP_DIR: Final[Path] = PROJECT_DIR / "temp"


# =============================================================================
# THRESHOLD-INDEPENDENT SCIENTIFIC INPUTS
# =============================================================================

MASK_FILE: Final[Path] = MASK_DIR / "IndonesiaNew0.msk"

PACIFIC_MASK_OISST_FILE: Final[Path] = (
    PROCESSED_DIR / "pacific_mask_oisst.npy"
)
PACIFIC_MASK_FILE: Final[Path] = PACIFIC_MASK_OISST_FILE

GRID_LAT_FILE: Final[Path] = PROCESSED_DIR / "grid_lat.npy"
GRID_LON_FILE: Final[Path] = PROCESSED_DIR / "grid_lon.npy"

LATITUDE_GRID_FILE: Final[Path] = GRID_LAT_FILE
LONGITUDE_GRID_FILE: Final[Path] = GRID_LON_FILE

# Legacy aliases retained for early-stage scripts.
MASK_NPY: Final[Path] = PROCESSED_DIR / "pacific_mask.npy"
PWP_OUTPUT: Final[Path] = PROCESSED_DIR / "pwp.txt"

# =============================================================================
# EXTERNAL CLIMATE INDICES — PROGRAM 21 DATA CONTRACT
# =============================================================================

# Immutable raw sources.
RAW_NINO_FILE: Final[Path] = (
    RAW_NINO_DIR / "SSTOI_NINOS.txt"
)
RAW_SOI_FILE: Final[Path] = (
    RAW_SOI_DIR / "monthly_soi.txt"
)

# Canonical Program-21 processed products.
WEEKLY_NINO_INDICES_FILE: Final[Path] = (
    PROCESSED_NINO_DIR / "weekly_nino_indices.csv"
)
WEEKLY_NINO_INDICES_QC_FILE: Final[Path] = (
    PROCESSED_NINO_DIR / "weekly_nino_indices_qc.csv"
)
WEEKLY_NINO_INDICES_METADATA_FILE: Final[Path] = (
    PROCESSED_NINO_DIR / "weekly_nino_indices_metadata.json"
)

MONTHLY_SOI_FILE: Final[Path] = (
    PROCESSED_SOI_DIR / "monthly_soi.csv"
)
MONTHLY_SOI_QC_FILE: Final[Path] = (
    PROCESSED_SOI_DIR / "monthly_soi_qc.csv"
)
MONTHLY_SOI_METADATA_FILE: Final[Path] = (
    PROCESSED_SOI_DIR / "monthly_soi_metadata.json"
)

# Historical monthly-Niño alias retained only for archived scripts. Program 21
# currently provides authoritative weekly SST/SSTA Niño products.
MONTHLY_NINO_INDICES_FILE: Final[Path] = (
    LEGACY_EXTERNAL_INDEX_DIR / "monthly_nino_indices.csv"
)


# =============================================================================
# OISST GRID AND EARTH CONSTANTS
# =============================================================================

LAT_MIN: Final[float] = -89.875
LAT_MAX: Final[float] = 89.875
LON_MIN: Final[float] = 0.125
LON_MAX: Final[float] = 359.875

DLAT: Final[float] = 0.25
DLON: Final[float] = 0.25

NLAT: Final[int] = 720
NLON: Final[int] = 1440

EARTH_RADIUS_KM: Final[float] = 6371.0088
EARTH_RADIUS: Final[float] = EARTH_RADIUS_KM

MASK_ROWS: Final[int] = 810
MASK_COLS: Final[int] = 80
MASK_OCEAN: Final[int] = 1
MASK_LAND: Final[int] = 0


# =============================================================================
# THRESHOLD CONFIGURATION
# =============================================================================

RUN_ALL_THRESHOLDS: bool = True

PWP_SST_THRESHOLD_C: float = 29.0

PWP_ALLOWED_THRESHOLDS_C: Final[tuple[float, ...]] = (
    28.0,
    28.5,
    29.0,
)

SST_THRESHOLD: float = PWP_SST_THRESHOLD_C
PWP_THRESHOLDS: Final[tuple[float, ...]] = PWP_ALLOWED_THRESHOLDS_C


def normalize_threshold(
    threshold_c: float,
) -> float:
    """Validate and return one approved PWP SST threshold."""

    try:
        threshold = float(threshold_c)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "PWP SST threshold must be numeric."
        ) from error

    for allowed in PWP_ALLOWED_THRESHOLDS_C:
        if abs(threshold - allowed) <= 1.0e-10:
            return float(allowed)

    allowed_text = ", ".join(
        f"{value:g}"
        for value in PWP_ALLOWED_THRESHOLDS_C
    )

    raise ValueError(
        "Unsupported PWP SST threshold.\n"
        f"Selected : {threshold:g} °C\n"
        f"Allowed  : {allowed_text} °C"
    )


def threshold_folder_name(
    threshold_c: float,
) -> str:
    """Convert an approved threshold to a readable folder name."""

    threshold = normalize_threshold(threshold_c)

    return (
        f"{threshold:.1f}"
        .rstrip("0")
        .rstrip(".")
    )


def threshold_tag(
    threshold_c: float,
) -> str:
    """Convert an approved threshold to a filename-safe tag."""

    return threshold_folder_name(
        threshold_c
    ).replace(
        ".",
        "_",
    )


def thresholds_to_run() -> tuple[float, ...]:
    """Return thresholds selected by the active execution mode."""

    if RUN_ALL_THRESHOLDS:
        return tuple(
            normalize_threshold(value)
            for value in PWP_ALLOWED_THRESHOLDS_C
        )

    return (
        normalize_threshold(PWP_SST_THRESHOLD_C),
    )


PWP_THRESHOLD_FOLDER: str = threshold_folder_name(
    PWP_SST_THRESHOLD_C
)
PWP_THRESHOLD_TAG: str = threshold_tag(
    PWP_SST_THRESHOLD_C
)


# =============================================================================
# SCIENTIFIC MODULE REGISTRY
# =============================================================================

PIPELINE_MODULES: Final[tuple[str, ...]] = (
    "centroid",
    "methodology",
    "quality_control",
    "descriptive_statistics",
    "distribution",
    "trend",
    "seasonality",
    "stl",
    "stationarity",
    "autocorrelation",
    "spectral",
    "wavelet",
    "real_wavelet",
    "cross_wavelet",
    "cross_wavelet_nino_soi",
    "spectral_synthesis",
    "models",
    "forecast",
    "comparison",
)

# Modules that produce threshold-specific processed numerical data.
PROCESSED_MODULES: Final[tuple[str, ...]] = (
    "centroid",
    "quality_control",
    "descriptive_statistics",
    "distribution",
    "trend",
    "seasonality",
    "stl",
    "stationarity",
    "autocorrelation",
    "spectral",
    "wavelet",
    "real_wavelet",
    "cross_wavelet",
    "cross_wavelet_nino_soi",
    "spectral_synthesis",
    "models",
    "forecast",
)

# Modules expected to produce cross-threshold products.
THRESHOLD_COMPARISON_MODULES: Final[tuple[str, ...]] = (
    "descriptive_statistics",
    "distribution",
    "trend",
    "seasonality",
    "stl",
    "stationarity",
    "autocorrelation",
    "spectral",
    "wavelet",
    "real_wavelet",
    "cross_wavelet",
    "cross_wavelet_nino_soi",
    "spectral_synthesis",
    "models",
    "forecast",
    "comparison",
)


# =============================================================================
# PATH DATA MODELS
# =============================================================================

@dataclass(frozen=True)
class ModulePaths:
    """Mirrored directories for one threshold-specific module."""

    name: str
    processed_dir: Path
    figure_dir: Path
    table_dir: Path
    report_dir: Path
    log_dir: Path

    def all_directories(self) -> tuple[Path, ...]:
        """Return all directories belonging to this module."""

        return (
            self.processed_dir,
            self.figure_dir,
            self.table_dir,
            self.report_dir,
            self.log_dir,
        )


@dataclass(frozen=True)
class ComparisonModulePaths:
    """Project-level threshold-comparison directories for one module."""

    name: str
    figure_dir: Path
    table_dir: Path
    report_dir: Path

    def all_directories(self) -> tuple[Path, ...]:
        """Return all comparison directories for this module."""

        return (
            self.figure_dir,
            self.table_dir,
            self.report_dir,
        )


@dataclass(frozen=True)
class ThresholdPaths:
    """Complete threshold-specific path architecture."""

    threshold_c: float
    folder_name: str
    tag: str

    processed_dir: Path
    figure_dir: Path
    table_dir: Path
    report_dir: Path
    log_dir: Path

    centroid: ModulePaths
    methodology: ModulePaths
    quality_control: ModulePaths
    descriptive_statistics: ModulePaths
    distribution: ModulePaths
    trend: ModulePaths
    seasonality: ModulePaths
    stl: ModulePaths
    stationarity: ModulePaths
    autocorrelation: ModulePaths
    spectral: ModulePaths
    wavelet: ModulePaths
    real_wavelet: ModulePaths
    cross_wavelet: ModulePaths
    cross_wavelet_nino_soi: ModulePaths
    spectral_synthesis: ModulePaths
    models: ModulePaths
    forecast: ModulePaths
    comparison: ModulePaths

    centroid_series_csv: Path
    centroid_qc_csv: Path

    processing_summary_txt: Path
    scientific_methodology_txt: Path
    processing_metadata_json: Path

    quality_control_report_txt: Path
    stl_report_txt: Path
    spectral_report_txt: Path
    wavelet_report_txt: Path
    real_wavelet_report_txt: Path
    cross_wavelet_report_txt: Path
    cross_wavelet_nino_soi_report_txt: Path
    spectral_synthesis_report_txt: Path
    distribution_report_txt: Path

    def module(
        self,
        module_name: str,
    ) -> ModulePaths:
        """Return one registered module path object."""

        normalized = str(module_name).strip().lower()

        if normalized not in PIPELINE_MODULES:
            raise KeyError(
                "Unknown pipeline module.\n"
                f"Requested : {module_name}\n"
                f"Available : {', '.join(PIPELINE_MODULES)}"
            )

        return getattr(self, normalized)

    def modules(self) -> tuple[ModulePaths, ...]:
        """Return all module path objects in registry order."""

        return tuple(
            self.module(name)
            for name in PIPELINE_MODULES
        )

    def root_directories(self) -> tuple[Path, ...]:
        """Return threshold-level root directories."""

        return (
            self.processed_dir,
            self.figure_dir,
            self.table_dir,
            self.report_dir,
            self.log_dir,
        )

    def all_directories(self) -> tuple[Path, ...]:
        """Return threshold roots and all mirrored module directories."""

        directories: list[Path] = list(
            self.root_directories()
        )

        for module_paths in self.modules():
            directories.extend(
                module_paths.all_directories()
            )

        return tuple(dict.fromkeys(directories))


def _build_module_paths(
    module_name: str,
    processed_root: Path,
    figure_root: Path,
    table_root: Path,
    report_root: Path,
    log_root: Path,
) -> ModulePaths:
    """Build the mirrored path object for one module."""

    normalized = str(module_name).strip().lower()

    if normalized not in PIPELINE_MODULES:
        raise ValueError(
            f"Unsupported pipeline module: {module_name}"
        )

    return ModulePaths(
        name=normalized,
        processed_dir=processed_root / normalized,
        figure_dir=figure_root / normalized,
        table_dir=table_root / normalized,
        report_dir=report_root / normalized,
        log_dir=log_root / normalized,
    )


def get_threshold_paths(
    threshold_c: float,
) -> ThresholdPaths:
    """Build all paths associated with one approved threshold."""

    threshold = normalize_threshold(threshold_c)
    folder_name = threshold_folder_name(threshold)
    tag = threshold_tag(threshold)

    processed_dir = PROCESSED_DIR / folder_name
    figure_dir = FIGURE_DIR / folder_name
    table_dir = TABLE_DIR / folder_name
    report_dir = REPORT_DIR / folder_name
    log_dir = LOG_DIR / folder_name

    module_paths = {
        module_name: _build_module_paths(
            module_name=module_name,
            processed_root=processed_dir,
            figure_root=figure_dir,
            table_root=table_dir,
            report_root=report_dir,
            log_root=log_dir,
        )
        for module_name in PIPELINE_MODULES
    }

    centroid = module_paths["centroid"]
    quality_control = module_paths["quality_control"]
    distribution = module_paths["distribution"]
    stl = module_paths["stl"]
    spectral = module_paths["spectral"]
    wavelet = module_paths["wavelet"]
    real_wavelet = module_paths["real_wavelet"]
    cross_wavelet = module_paths["cross_wavelet"]
    cross_wavelet_nino_soi = module_paths[
        "cross_wavelet_nino_soi"
    ]
    spectral_synthesis = module_paths[
        "spectral_synthesis"
    ]

    return ThresholdPaths(
        threshold_c=threshold,
        folder_name=folder_name,
        tag=tag,
        processed_dir=processed_dir,
        figure_dir=figure_dir,
        table_dir=table_dir,
        report_dir=report_dir,
        log_dir=log_dir,
        centroid=centroid,
        methodology=module_paths["methodology"],
        quality_control=quality_control,
        descriptive_statistics=module_paths[
            "descriptive_statistics"
        ],
        distribution=distribution,
        trend=module_paths["trend"],
        seasonality=module_paths["seasonality"],
        stl=stl,
        stationarity=module_paths["stationarity"],
        autocorrelation=module_paths["autocorrelation"],
        spectral=spectral,
        wavelet=wavelet,
        real_wavelet=real_wavelet,
        cross_wavelet=cross_wavelet,
        cross_wavelet_nino_soi=cross_wavelet_nino_soi,
        spectral_synthesis=spectral_synthesis,
        models=module_paths["models"],
        forecast=module_paths["forecast"],
        comparison=module_paths["comparison"],
        centroid_series_csv=(
            centroid.processed_dir
            / "pwp_centroid_series.csv"
        ),
        centroid_qc_csv=(
            quality_control.processed_dir
            / "pwp_centroid_series_qc.csv"
        ),
        processing_summary_txt=(
            centroid.report_dir
            / "pwp_processing_summary.txt"
        ),
        scientific_methodology_txt=(
            centroid.report_dir
            / "pwp_methodology_scientific_report.txt"
        ),
        processing_metadata_json=(
            centroid.report_dir
            / "pwp_processing_metadata.json"
        ),
        quality_control_report_txt=(
            quality_control.report_dir
            / "pwp_quality_control_report.txt"
        ),
        stl_report_txt=(
            stl.report_dir
            / "pwp_stl_decomposition_report.txt"
        ),
        spectral_report_txt=(
            spectral.report_dir
            / "pwp_spectral_analysis_report.txt"
        ),
        wavelet_report_txt=(
            wavelet.report_dir
            / "pwp_wavelet_analysis_report.txt"
        ),
        real_wavelet_report_txt=(
            real_wavelet.report_dir
            / "pwp_wavelet_real_part_report.txt"
        ),
        cross_wavelet_report_txt=(
            cross_wavelet.report_dir
            / "pwp_cross_wavelet_report.txt"
        ),
        cross_wavelet_nino_soi_report_txt=(
            cross_wavelet_nino_soi.report_dir
            / "pwp_cross_wavelet_nino_soi_report.txt"
        ),
        spectral_synthesis_report_txt=(
            spectral_synthesis.report_dir
            / "scientific_summary.txt"
        ),
        distribution_report_txt=(
            distribution.report_dir
            / "pwp_distribution_analysis_report.txt"
        ),
    )


def get_active_threshold_paths() -> ThresholdPaths:
    """Return the active single-threshold path object."""

    return get_threshold_paths(
        PWP_SST_THRESHOLD_C
    )


def get_module_paths(
    threshold_c: float,
    module_name: str,
) -> ModulePaths:
    """Return one module path object for one threshold."""

    return get_threshold_paths(
        threshold_c
    ).module(
        module_name
    )


# =============================================================================
# CROSS-THRESHOLD PATHS
# =============================================================================

THRESHOLD_COMPARISON_TAG: Final[str] = "threshold_comparison"

THRESHOLD_COMPARISON_FIGURE_DIR: Final[Path] = (
    FIGURE_DIR / THRESHOLD_COMPARISON_TAG
)
THRESHOLD_COMPARISON_TABLE_DIR: Final[Path] = (
    TABLE_DIR / THRESHOLD_COMPARISON_TAG
)
THRESHOLD_COMPARISON_REPORT_DIR: Final[Path] = (
    REPORT_DIR / THRESHOLD_COMPARISON_TAG
)


def get_comparison_module_paths(
    module_name: str,
) -> ComparisonModulePaths:
    """Return cross-threshold directories for one registered module."""

    normalized = str(module_name).strip().lower()

    if normalized not in THRESHOLD_COMPARISON_MODULES:
        raise KeyError(
            "Unknown threshold-comparison module.\n"
            f"Requested : {module_name}\n"
            "Available : "
            f"{', '.join(THRESHOLD_COMPARISON_MODULES)}"
        )

    return ComparisonModulePaths(
        name=normalized,
        figure_dir=(
            THRESHOLD_COMPARISON_FIGURE_DIR
            / normalized
        ),
        table_dir=(
            THRESHOLD_COMPARISON_TABLE_DIR
            / normalized
        ),
        report_dir=(
            THRESHOLD_COMPARISON_REPORT_DIR
            / normalized
        ),
    )


def comparison_module_paths() -> tuple[
    ComparisonModulePaths,
    ...,
]:
    """Return all cross-threshold module path objects."""

    return tuple(
        get_comparison_module_paths(module_name)
        for module_name in THRESHOLD_COMPARISON_MODULES
    )


# =============================================================================
# PROJECT-LEVEL MODULES
# =============================================================================

METHODOLOGICAL_DOMAIN_TAG: Final[str] = "methodological_domain"
PIPELINE_SYNTHESIS_TAG: Final[str] = "pipeline_synthesis"

METHODOLOGICAL_DOMAIN_PROCESSED_DIR: Final[Path] = (
    PROCESSED_DIR / METHODOLOGICAL_DOMAIN_TAG
)
METHODOLOGICAL_DOMAIN_FIGURE_DIR: Final[Path] = (
    FIGURE_DIR / METHODOLOGICAL_DOMAIN_TAG
)
METHODOLOGICAL_DOMAIN_TABLE_DIR: Final[Path] = (
    TABLE_DIR / METHODOLOGICAL_DOMAIN_TAG
)
METHODOLOGICAL_DOMAIN_REPORT_DIR: Final[Path] = (
    REPORT_DIR / METHODOLOGICAL_DOMAIN_TAG
)
METHODOLOGICAL_DOMAIN_LOG_DIR: Final[Path] = (
    LOG_DIR / METHODOLOGICAL_DOMAIN_TAG
)

PIPELINE_SYNTHESIS_FIGURE_DIR: Final[Path] = (
    FIGURE_DIR / PIPELINE_SYNTHESIS_TAG
)
PIPELINE_SYNTHESIS_TABLE_DIR: Final[Path] = (
    TABLE_DIR / PIPELINE_SYNTHESIS_TAG
)
PIPELINE_SYNTHESIS_REPORT_DIR: Final[Path] = (
    REPORT_DIR / PIPELINE_SYNTHESIS_TAG
)
PIPELINE_SYNTHESIS_LOG_DIR: Final[Path] = (
    LOG_DIR / PIPELINE_SYNTHESIS_TAG
)


# =============================================================================
# PUBLICATION FIGURE STANDARDS
# =============================================================================

FIGURE_DPI: Final[int] = 300

SAVE_BBOX: Final[str] = "tight"
SAVE_PAD_INCHES: Final[float] = 0.05
SAVE_TRANSPARENT: Final[bool] = False

GRID_LINESTYLE: Final[str] = "--"
GRID_LINEWIDTH: Final[float] = 0.55
GRID_ALPHA: Final[float] = 0.35

PANEL_LABEL_FONT_SIZE: Final[int] = 12
PANEL_LABEL_FONT_WEIGHT: Final[str] = "bold"
PANEL_LABEL_X: Final[float] = 0.012
PANEL_LABEL_Y: Final[float] = 0.975

YEAR_TICK_INTERVAL: Final[int] = 4
YEAR_TICK_ANCHOR: Final[int] = 1982
TIME_AXIS_LABEL: Final[str] = "Year"

DEFAULT_FIGURE_WIDTH_INCHES: Final[float] = 14.0
DEFAULT_FIGURE_HEIGHT_INCHES: Final[float] = 8.0


# =============================================================================
# WAVELET AND TEMPORAL-STANDARD CONSTANTS
# =============================================================================

DAILY_SAMPLING_INTERVAL_DAYS: Final[float] = 1.0
WEEKLY_SAMPLING_INTERVAL_DAYS: Final[float] = 7.0

MORLET_OMEGA0: Final[float] = 6.0
WAVELET_DJ: Final[float] = 1.0 / 12.0
WAVELET_SIGNIFICANCE_LEVEL: Final[float] = 0.95

ANNUAL_BAND_DAYS: Final[tuple[float, float]] = (
    330.0,
    400.0,
)
INTERANNUAL_BAND_DAYS: Final[tuple[float, float]] = (
    2.0 * 365.2425,
    8.0 * 365.2425,
)
INTERDECADAL_BAND_DAYS: Final[tuple[float, float]] = (
    8.0 * 365.2425,
    16.0 * 365.2425,
)

WEEKLY_PWP_FREQUENCY: Final[str] = "W-WED"
MONTHLY_PWP_FREQUENCY: Final[str] = "ME"


# =============================================================================
# DIRECTORY CREATION HELPERS
# =============================================================================

def base_directories() -> tuple[Path, ...]:
    """Return all threshold-independent base directories."""

    return (
        CONFIG_DIR,
        SRC_DIR,
        DOCUMENTATION_DIR,
        DOCUMENTATION_VALIDATION_DIR,
        DOCUMENTATION_PROGRAM21_DIR,
        DOCUMENTATION_PROGRAM19_NINO_SOI_DIR,
        TESTS_DIR,
        DATA_DIR,
        MASK_DIR,
        RAW_DATA_DIR,
        RAW_CLIMATE_INDICES_DIR,
        RAW_NINO_DIR,
        RAW_SOI_DIR,
        PROCESSED_DIR,
        PROCESSED_CLIMATE_INDICES_DIR,
        PROCESSED_NINO_DIR,
        PROCESSED_SOI_DIR,
        PREPARED_DIR,
        LEGACY_EXTERNAL_INDEX_DIR,
        OUTPUT_DIR,
        FIGURE_DIR,
        TABLE_DIR,
        REPORT_DIR,
        LOG_DIR,
        MODEL_DIR,
        FORECAST_DIR,
        CLIMATE_INDEX_TABLE_DIR,
        CLIMATE_INDEX_REPORT_DIR,
        ARCHIVE_DIR,
        TEMP_DIR,
        THRESHOLD_COMPARISON_FIGURE_DIR,
        THRESHOLD_COMPARISON_TABLE_DIR,
        THRESHOLD_COMPARISON_REPORT_DIR,
        METHODOLOGICAL_DOMAIN_PROCESSED_DIR,
        METHODOLOGICAL_DOMAIN_FIGURE_DIR,
        METHODOLOGICAL_DOMAIN_TABLE_DIR,
        METHODOLOGICAL_DOMAIN_REPORT_DIR,
        METHODOLOGICAL_DOMAIN_LOG_DIR,
        PIPELINE_SYNTHESIS_FIGURE_DIR,
        PIPELINE_SYNTHESIS_TABLE_DIR,
        PIPELINE_SYNTHESIS_REPORT_DIR,
        PIPELINE_SYNTHESIS_LOG_DIR,
    )


def all_project_directories() -> tuple[Path, ...]:
    """Return the complete deterministic project-directory specification."""

    directories: list[Path] = list(
        base_directories()
    )

    for threshold in PWP_ALLOWED_THRESHOLDS_C:
        directories.extend(
            get_threshold_paths(
                threshold
            ).all_directories()
        )

    for comparison_paths in comparison_module_paths():
        directories.extend(
            comparison_paths.all_directories()
        )

    return tuple(
        dict.fromkeys(directories)
    )


def ensure_directories(
    directories: Iterable[Path],
) -> tuple[Path, ...]:
    """Create directories idempotently and return normalized unique paths."""

    normalized = tuple(
        dict.fromkeys(
            Path(directory)
            for directory in directories
        )
    )

    for directory in normalized:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    return normalized


def ensure_base_directories() -> tuple[Path, ...]:
    """Create threshold-independent base directories."""

    return ensure_directories(
        base_directories()
    )


def ensure_threshold_directories(
    threshold_c: float,
) -> ThresholdPaths:
    """Create all threshold-specific directories for one threshold."""

    paths = get_threshold_paths(
        threshold_c
    )

    ensure_directories(
        paths.all_directories()
    )

    return paths


def ensure_comparison_directories() -> tuple[Path, ...]:
    """Create all registered cross-threshold directories."""

    directories: list[Path] = [
        THRESHOLD_COMPARISON_FIGURE_DIR,
        THRESHOLD_COMPARISON_TABLE_DIR,
        THRESHOLD_COMPARISON_REPORT_DIR,
    ]

    for paths in comparison_module_paths():
        directories.extend(
            paths.all_directories()
        )

    return ensure_directories(
        directories
    )


def ensure_complete_project_structure() -> tuple[Path, ...]:
    """Create the complete directory architecture without deleting data."""

    return ensure_directories(
        all_project_directories()
    )


# =============================================================================
# VALIDATION
# =============================================================================

def _is_relative_to(
    path: Path,
    parent: Path,
) -> bool:
    """Compatibility-safe Path.is_relative_to implementation."""

    try:
        path.resolve().relative_to(
            parent.resolve()
        )
        return True
    except ValueError:
        return False


def validate_project_configuration() -> tuple[str, ...]:
    """
    Validate internal configuration consistency.

    Returns
    -------
    tuple[str, ...]
        Human-readable validation messages.

    Raises
    ------
    ValueError
        If one or more structural errors are detected.
    """

    errors: list[str] = []
    messages: list[str] = []

    expected_config_dir = PROJECT_DIR / "config"
    expected_config_file = expected_config_dir / "config.py"

    if CONFIG_DIR.resolve() != expected_config_dir.resolve():
        errors.append(
            "CONFIG_DIR does not equal PROJECT_DIR/config."
        )

    if CONFIG_FILE.resolve() != expected_config_file.resolve():
        errors.append(
            "CONFIG_FILE does not equal PROJECT_DIR/config/config.py."
        )

    canonical_climate_paths = (
        RAW_CLIMATE_INDICES_DIR,
        RAW_NINO_DIR,
        RAW_SOI_DIR,
        PROCESSED_CLIMATE_INDICES_DIR,
        PROCESSED_NINO_DIR,
        PROCESSED_SOI_DIR,
        CLIMATE_INDEX_TABLE_DIR,
        CLIMATE_INDEX_REPORT_DIR,
    )

    for climate_path in canonical_climate_paths:
        if not _is_relative_to(
            climate_path,
            PROJECT_DIR,
        ):
            errors.append(
                "Configured climate-index path lies outside PROJECT_DIR: "
                f"{climate_path}"
            )

    if WEEKLY_NINO_INDICES_FILE.parent != PROCESSED_NINO_DIR:
        errors.append(
            "WEEKLY_NINO_INDICES_FILE is not located in PROCESSED_NINO_DIR."
        )

    if MONTHLY_SOI_FILE.parent != PROCESSED_SOI_DIR:
        errors.append(
            "MONTHLY_SOI_FILE is not located in PROCESSED_SOI_DIR."
        )

    if not PWP_ALLOWED_THRESHOLDS_C:
        errors.append(
            "PWP_ALLOWED_THRESHOLDS_C cannot be empty."
        )

    if len(
        set(PWP_ALLOWED_THRESHOLDS_C)
    ) != len(
        PWP_ALLOWED_THRESHOLDS_C
    ):
        errors.append(
            "PWP_ALLOWED_THRESHOLDS_C contains duplicate values."
        )

    try:
        normalize_threshold(
            PWP_SST_THRESHOLD_C
        )
    except ValueError as error:
        errors.append(str(error))

    if len(
        set(PIPELINE_MODULES)
    ) != len(
        PIPELINE_MODULES
    ):
        errors.append(
            "PIPELINE_MODULES contains duplicate names."
        )

    if len(
        set(THRESHOLD_COMPARISON_MODULES)
    ) != len(
        THRESHOLD_COMPARISON_MODULES
    ):
        errors.append(
            "THRESHOLD_COMPARISON_MODULES contains duplicate names."
        )

    unknown_comparison_modules = (
        set(THRESHOLD_COMPARISON_MODULES)
        - set(PIPELINE_MODULES)
    )

    if unknown_comparison_modules:
        errors.append(
            "Threshold-comparison modules are not registered in "
            f"PIPELINE_MODULES: {sorted(unknown_comparison_modules)}"
        )

    for directory in all_project_directories():
        if not _is_relative_to(
            directory,
            PROJECT_DIR,
        ):
            errors.append(
                "Configured directory lies outside PROJECT_DIR: "
                f"{directory}"
            )

    for threshold in PWP_ALLOWED_THRESHOLDS_C:
        paths = get_threshold_paths(
            threshold
        )

        if paths.folder_name != threshold_folder_name(
            threshold
        ):
            errors.append(
                f"Threshold folder-name mismatch for {threshold}."
            )

        if paths.tag != threshold_tag(
            threshold
        ):
            errors.append(
                f"Threshold tag mismatch for {threshold}."
            )

        for module_name in PIPELINE_MODULES:
            module_paths = paths.module(
                module_name
            )

            if module_paths.name != module_name:
                errors.append(
                    "Module-path name mismatch: "
                    f"{module_name} / {module_paths.name}"
                )

    if errors:
        raise ValueError(
            "Project configuration validation failed:\n"
            + "\n".join(
                f"  - {message}"
                for message in errors
            )
        )

    messages.extend(
        (
            f"Project root validated: {PROJECT_DIR}",
            f"Thresholds validated: {PWP_ALLOWED_THRESHOLDS_C}",
            f"Pipeline modules validated: {len(PIPELINE_MODULES)}",
            (
                "Threshold-comparison modules validated: "
                f"{len(THRESHOLD_COMPARISON_MODULES)}"
            ),
            (
                "Configured directories validated: "
                f"{len(all_project_directories())}"
            ),
        )
    )

    return tuple(messages)


# =============================================================================
# HUMAN-READABLE CONFIGURATION SUMMARY
# =============================================================================

def configuration_summary_lines() -> tuple[str, ...]:
    """Return a complete human-readable configuration summary."""

    selected_thresholds = ", ".join(
        f"{value:.1f} °C"
        for value in thresholds_to_run()
    )

    return (
        f"Project name               : {PROJECT_NAME}",
        f"Project version            : {PROJECT_VERSION}",
        f"Configuration file         : {CONFIG_FILE}",
        f"Project root               : {PROJECT_DIR}",
        f"Source directory           : {SRC_DIR}",
        f"Raw OISST directory        : {RAW_OISST_DIR}",
        f"Raw climate-index root     : {RAW_CLIMATE_INDICES_DIR}",
        f"Raw Niño file              : {RAW_NINO_FILE}",
        f"Raw SOI file               : {RAW_SOI_FILE}",
        f"Processed-data directory   : {PROCESSED_DIR}",
        f"Processed climate indices  : {PROCESSED_CLIMATE_INDICES_DIR}",
        f"Weekly Niño product        : {WEEKLY_NINO_INDICES_FILE}",
        f"Monthly SOI product        : {MONTHLY_SOI_FILE}",
        f"Figure directory           : {FIGURE_DIR}",
        f"Table directory            : {TABLE_DIR}",
        f"Report directory           : {REPORT_DIR}",
        f"Log directory              : {LOG_DIR}",
        (
            "Execution mode             : "
            f"{'all thresholds' if RUN_ALL_THRESHOLDS else 'single threshold'}"
        ),
        f"Selected thresholds        : {selected_thresholds}",
        (
            "Registered modules         : "
            f"{', '.join(PIPELINE_MODULES)}"
        ),
        (
            "Comparison modules         : "
            f"{', '.join(THRESHOLD_COMPARISON_MODULES)}"
        ),
        f"Weekly Niño input          : {WEEKLY_NINO_INDICES_FILE}",
        f"Monthly SOI input          : {MONTHLY_SOI_FILE}",
    )


def configuration_as_dict() -> dict[str, object]:
    """Return selected public configuration values as a serializable mapping."""

    return {
        "project_name": PROJECT_NAME,
        "project_version": PROJECT_VERSION,
        "project_root": str(PROJECT_DIR),
        "configuration_file": str(CONFIG_FILE),
        "run_all_thresholds": RUN_ALL_THRESHOLDS,
        "active_threshold_c": PWP_SST_THRESHOLD_C,
        "allowed_thresholds_c": list(
            PWP_ALLOWED_THRESHOLDS_C
        ),
        "pipeline_modules": list(
            PIPELINE_MODULES
        ),
        "threshold_comparison_modules": list(
            THRESHOLD_COMPARISON_MODULES
        ),
        "climate_indices": {
            "raw_root": str(
                RAW_CLIMATE_INDICES_DIR
            ),
            "raw_nino_file": str(
                RAW_NINO_FILE
            ),
            "raw_soi_file": str(
                RAW_SOI_FILE
            ),
            "processed_root": str(
                PROCESSED_CLIMATE_INDICES_DIR
            ),
            "weekly_nino_file": str(
                WEEKLY_NINO_INDICES_FILE
            ),
            "monthly_soi_file": str(
                MONTHLY_SOI_FILE
            ),
            "table_dir": str(
                CLIMATE_INDEX_TABLE_DIR
            ),
            "report_dir": str(
                CLIMATE_INDEX_REPORT_DIR
            ),
        },
        "raw_oisst_directory": str(
            RAW_OISST_DIR
        ),
        "weekly_nino_indices_file": str(
            WEEKLY_NINO_INDICES_FILE
        ),
        "monthly_soi_file": str(
            MONTHLY_SOI_FILE
        ),
    }


# =============================================================================
# DIRECT EXECUTION
# =============================================================================

def main() -> None:
    """Validate and print the authoritative project configuration."""

    messages = validate_project_configuration()

    print("=" * 78)
    print("PACIFIC WARM POOL PROJECT CONFIGURATION")
    print("=" * 78)

    for line in configuration_summary_lines():
        print(line)

    print()
    print("VALIDATION")
    print("-" * 78)

    for message in messages:
        print(message)

    print()
    print("=" * 78)
    print("CONFIGURATION VALIDATION COMPLETED SUCCESSFULLY.")
    print("=" * 78)


if __name__ == "__main__":
    main()