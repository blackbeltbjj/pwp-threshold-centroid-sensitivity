# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

PROGRAM
    05_calculate_pwp_centroid.py

VERSION
    3.2.0

PURPOSE
    Calculate daily Pacific Warm Pool area, area-weighted spherical centroid,
    internal SST statistics, spatial-concentration diagnostics, quality flags,
    and complete provenance for one or more approved SST thresholds.

DESIGN
    This program was written as a clean implementation for the portable project
    structure controlled by config/config.py.

    It does not depend on legacy centroid scripts or legacy fallback variables.

THRESHOLD MODES
    Single-threshold mode
        RUN_ALL_THRESHOLDS = False

        Only PWP_SST_THRESHOLD_C is processed.

    Threshold-sensitivity mode
        RUN_ALL_THRESHOLDS = True

        The approved thresholds are processed in the same OISST-reading cycle:

            28.0 °C
            28.5 °C
            29.0 °C

EFFICIENCY
    Every daily SST field is read and converted to degrees Celsius once.
    All selected thresholds are then evaluated from that same in-memory field.

SCIENTIFIC DEFINITION
    For threshold T, the daily PWP consists of cells satisfying:

        validated Pacific Ocean mask
        AND finite SST
        AND SST >= T

GRID-CELL AREA
    Exact spherical latitude-longitude cell area:

        A_i = R² Δλ_i [sin(phi_north,i) - sin(phi_south,i)]

SPHERICAL CENTROID
    Each selected grid-cell centre is transformed to a Cartesian unit vector:

        x_i = cos(phi_i) cos(lambda_i)
        y_i = cos(phi_i) sin(lambda_i)
        z_i = sin(phi_i)

    Area-weighted sums:

        X = Σ A_i x_i
        Y = Σ A_i y_i
        Z = Σ A_i z_i

    Centroid:

        lambda_c = atan2(Y, X)
        phi_c    = atan2(Z, sqrt(X² + Y²))

    This avoids arithmetic averaging of circular longitude and remains valid
    across the 0°/360° discontinuity.

INPUTS
    Configuration:
        config/config.py

    Daily absolute NOAA OISST:
        data/raw/**/*.nc
        data/raw/**/*.nc4
        data/raw/**/*.cdf

    Threshold-independent processed inputs:
        data/processed/pacific_mask_oisst.npy
        data/processed/grid_lat.npy
        data/processed/grid_lon.npy

OUTPUTS FOR EACH THRESHOLD
    Daily centroid series:
        data/processed/<threshold>/centroid/
            pwp_centroid_series.csv

    Technical report:
        outputs/reports/<threshold>/centroid/
            pwp_processing_summary.txt

    Technical-scientific methodology:
        outputs/reports/<threshold>/centroid/
            pwp_methodology_scientific_report.txt

    Machine-readable provenance:
        outputs/reports/<threshold>/centroid/
            pwp_processing_metadata.json

PATH POLICY
    All threshold-dependent paths are obtained from:

        get_threshold_paths(threshold).centroid

    The program never constructs threshold directories independently.

CROSS-THRESHOLD OUTPUTS
    When more than one threshold is selected:

        outputs/tables/threshold_comparison/centroid/
            pwp_centroid_threshold_comparison_summary.csv

        outputs/reports/threshold_comparison/centroid/
            pwp_centroid_threshold_comparison_report.txt

OUTPUT COLUMNS
    date
    lon_360
    lon_180
    lat
    area_km2
    pwp_cell_count
    valid_ocean_cell_count
    mask_cell_count
    ocean_data_coverage_pct
    mean_pwp_sst_c
    max_pwp_sst_c
    area_weighted_sst_std_c
    threshold_c
    threshold_source
    threshold_fallback_used
    centroid_resultant_length
    centroid_angular_dispersion_deg
    quality_flag
    source_file

QUALITY CONTROL
    - Absolute-SST files only.
    - SST-anomaly products are excluded by path indicators and rejected by
      variable metadata.
    - Grid shape, coordinates, and orientation are validated.
    - SST units are interpreted explicitly.
    - Implausible finite SST values are rejected.
    - Duplicate dates are rejected.
    - Missing calendar dates are recorded.
    - Low ocean-data coverage is flagged.
    - Threshold consistency is validated before export.

DIRECTORY POLICY
    Run 00_create_project_structure.py before this program.
    The program requires the threshold directories to exist and does not build
    the complete project tree itself.

DEPENDENCIES
    numpy
    pandas
    xarray

PYTHON
    Python 3.10+

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

import hashlib
import json
import sys
import warnings
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import xarray as xr


# =============================================================================
# PROJECT CONFIGURATION
# =============================================================================

SCRIPT_FILE = Path(__file__).resolve()
PROJECT_ROOT_FROM_SCRIPT = SCRIPT_FILE.parents[1]

if str(PROJECT_ROOT_FROM_SCRIPT) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_ROOT_FROM_SCRIPT),
    )

from config.config import (  # noqa: E402
    EARTH_RADIUS_KM,
    GRID_LAT_FILE,
    GRID_LON_FILE,
    PACIFIC_MASK_OISST_FILE,
    PROJECT_DIR,
    PWP_ALLOWED_THRESHOLDS_C,
    PWP_SST_THRESHOLD_C,
    RAW_DIR,
    RUN_ALL_THRESHOLDS,
    THRESHOLD_COMPARISON_REPORT_DIR,
    THRESHOLD_COMPARISON_TABLE_DIR,
    get_threshold_paths,
    normalize_threshold,
    threshold_folder_name,
    thresholds_to_run,
    validate_project_configuration,
)


# =============================================================================
# PROGRAM CONSTANTS
# =============================================================================

PROGRAM_NAME = "PACIFIC WARM POOL CENTROID AND AREA ANALYSIS"
PROGRAM_VERSION = "3.2.0"

MINIMUM_DAILY_COVERAGE_PCT = 99.0
COORDINATE_TOLERANCE = 1.0e-6

PLAUSIBLE_SST_MIN_C = -3.0
PLAUSIBLE_SST_MAX_C = 45.0

NETCDF_SUFFIXES = {
    ".nc",
    ".nc4",
    ".cdf",
}

SST_VARIABLE_CANDIDATES = (
    "sst",
    "sea_surface_temperature",
    "analysed_sst",
    "tos",
)

TIME_COORDINATE_CANDIDATES = (
    "time",
    "date",
    "day",
)

LATITUDE_COORDINATE_CANDIDATES = (
    "lat",
    "latitude",
    "y",
)

LONGITUDE_COORDINATE_CANDIDATES = (
    "lon",
    "longitude",
    "x",
)

ANOMALY_PATH_INDICATORS = (
    "anom",
    "anomaly",
    "ssta",
)

ABSOLUTE_SST_METADATA_REJECTION_TERMS = (
    "anomaly",
    "departure",
)

OUTPUT_COMPARISON_TABLE = (
    THRESHOLD_COMPARISON_TABLE_DIR
    / "centroid"
    / "pwp_centroid_threshold_comparison_summary.csv"
)

OUTPUT_COMPARISON_REPORT = (
    THRESHOLD_COMPARISON_REPORT_DIR
    / "centroid"
    / "pwp_centroid_threshold_comparison_report.txt"
)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class ThresholdContext:
    """Threshold value and its configured output paths."""

    threshold_c: float
    folder_name: str
    centroid_csv: Path
    technical_report: Path
    scientific_report: Path
    metadata_json: Path


@dataclass(frozen=True)
class InputFileMetadata:
    """Metadata recorded for one NetCDF input file."""

    source_file: str
    sst_variable: str
    time_coordinate: str
    latitude_coordinate: str
    longitude_coordinate: str
    original_units: str
    unit_interpretation: str
    records_per_threshold: int
    thresholds_processed_c: tuple[float, ...]


@dataclass(frozen=True)
class GridGeometry:
    """Precomputed exact cell areas and Cartesian unit vectors."""

    cell_area_km2: np.ndarray
    x_unit: np.ndarray
    y_unit: np.ndarray
    z_unit: np.ndarray
    global_area_km2: float


# =============================================================================
# CONSOLE UTILITIES
# =============================================================================

def print_rule(
    character: str = "=",
    width: int = 78,
) -> None:
    """Print a horizontal rule."""

    print(
        character
        * width
    )


def print_header(
    title: str,
) -> None:
    """Print the main program header."""

    print_rule()
    print(
        title
    )
    print_rule()


def print_section(
    title: str,
) -> None:
    """Print a terminal section heading."""

    print()
    print(
        title
    )
    print_rule(
        "-",
    )


def print_key_value(
    key: str,
    value: object,
) -> None:
    """Print an aligned key-value pair."""

    print(
        f"{key:<28s}: {value}"
    )


# =============================================================================
# CONFIGURATION AND PATHS
# =============================================================================

def selected_threshold_contexts() -> tuple[
    ThresholdContext,
    ...,
]:
    """Build validated threshold-specific output contexts."""

    contexts: list[
        ThresholdContext
    ] = []

    for threshold_c in thresholds_to_run():
        threshold = normalize_threshold(
            threshold_c
        )

        paths = get_threshold_paths(
            threshold
        )

        centroid_paths = paths.centroid

        contexts.append(
            ThresholdContext(
                threshold_c=threshold,
                folder_name=paths.folder_name,
                centroid_csv=(
                    centroid_paths.processed_dir
                    / "pwp_centroid_series.csv"
                ),
                technical_report=(
                    centroid_paths.report_dir
                    / "pwp_processing_summary.txt"
                ),
                scientific_report=(
                    centroid_paths.report_dir
                    / "pwp_methodology_scientific_report.txt"
                ),
                metadata_json=(
                    centroid_paths.report_dir
                    / "pwp_processing_metadata.json"
                ),
            )
        )

    return tuple(
        contexts
    )


def validate_configuration(
    contexts: tuple[
        ThresholdContext,
        ...,
    ],
) -> None:
    """Validate the central configuration and execution mode."""

    validate_project_configuration()

    if PROJECT_ROOT_FROM_SCRIPT.resolve() != Path(
        PROJECT_DIR
    ).resolve():
        raise ValueError(
            "Project-root mismatch.\n"
            f"Script-derived root : {PROJECT_ROOT_FROM_SCRIPT}\n"
            f"Configured root     : {PROJECT_DIR}"
        )

    if not contexts:
        raise ValueError(
            "No PWP thresholds were selected."
        )

    for context in contexts:
        configured_paths = get_threshold_paths(
            context.threshold_c
        )

        configured_centroid = (
            configured_paths.centroid
        )

        if (
            context.centroid_csv.parent.resolve()
            != configured_centroid.processed_dir.resolve()
        ):
            raise ValueError(
                "Centroid-series output is outside the configured centroid "
                "processed-data directory.\n"
                f"Configured : {configured_centroid.processed_dir}\n"
                f"Selected   : {context.centroid_csv.parent}"
            )

        for report_file in (
            context.technical_report,
            context.scientific_report,
            context.metadata_json,
        ):
            if (
                report_file.parent.resolve()
                != configured_centroid.report_dir.resolve()
            ):
                raise ValueError(
                    "Centroid report output is outside the configured "
                    "centroid report directory.\n"
                    f"Configured : {configured_centroid.report_dir}\n"
                    f"Selected   : {report_file.parent}"
                )

    selected = tuple(
        context.threshold_c
        for context in contexts
    )

    if len(
        selected
    ) != len(
        set(
            selected
        )
    ):
        raise ValueError(
            "Selected thresholds must be unique."
        )

    if RUN_ALL_THRESHOLDS:
        expected = tuple(
            float(
                value
            )
            for value in PWP_ALLOWED_THRESHOLDS_C
        )

        if selected != expected:
            raise ValueError(
                "All-threshold mode must process the approved thresholds "
                "in configured order."
            )

    if not np.isfinite(
        float(
            EARTH_RADIUS_KM
        )
    ):
        raise ValueError(
            "EARTH_RADIUS_KM must be finite."
        )

    if not (
        0.0
        < float(
            EARTH_RADIUS_KM
        )
        < 10_000.0
    ):
        raise ValueError(
            f"Invalid Earth radius: {EARTH_RADIUS_KM}"
        )

    if not (
        0.0
        <= MINIMUM_DAILY_COVERAGE_PCT
        <= 100.0
    ):
        raise ValueError(
            "MINIMUM_DAILY_COVERAGE_PCT must lie in [0, 100]."
        )


def validate_required_inputs_and_outputs(
    contexts: tuple[
        ThresholdContext,
        ...,
    ],
) -> None:
    """Validate required files and directories before expensive processing."""

    required_files = (
        PACIFIC_MASK_OISST_FILE,
        GRID_LAT_FILE,
        GRID_LON_FILE,
    )

    missing_files = [
        path
        for path in required_files
        if not path.is_file()
    ]

    if missing_files:
        raise FileNotFoundError(
            "Required threshold-independent input files are missing:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing_files
            )
        )

    if not RAW_DIR.is_dir():
        raise FileNotFoundError(
            "Raw OISST directory not found:\n"
            f"{RAW_DIR}"
        )

    required_directories: list[
        Path
    ] = []

    for context in contexts:
        required_directories.extend(
            [
                context.centroid_csv.parent,
                context.technical_report.parent,
                context.scientific_report.parent,
                context.metadata_json.parent,
            ]
        )

    if len(
        contexts
    ) > 1:
        required_directories.extend(
            [
                OUTPUT_COMPARISON_TABLE.parent,
                OUTPUT_COMPARISON_REPORT.parent,
            ]
        )

    missing_directories = sorted(
        {
            path
            for path in required_directories
            if not path.is_dir()
        },
        key=str,
    )

    if missing_directories:
        raise FileNotFoundError(
            "Required output directories are missing. "
            "Run 00_create_project_structure.py first:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing_directories
            )
        )


# =============================================================================
# FILE DISCOVERY AND PROVENANCE
# =============================================================================

def normalized_path_text(
    path: Path,
) -> str:
    """Return a lower-case slash-normalized path string."""

    return str(
        path
    ).lower().replace(
        "\\",
        "/",
    )


def is_anomaly_path(
    file_path: Path,
) -> bool:
    """Classify anomaly products using path indicators."""

    path_text = normalized_path_text(
        file_path
    )

    return any(
        indicator in path_text
        for indicator in ANOMALY_PATH_INDICATORS
    )


def discover_netcdf_files() -> tuple[
    list[Path],
    list[Path],
]:
    """Discover absolute-SST NetCDF files and separately list exclusions."""

    discovered = sorted(
        {
            path
            for path in RAW_DIR.rglob(
                "*"
            )
            if (
                path.is_file()
                and path.suffix.lower()
                in NETCDF_SUFFIXES
            )
        },
        key=lambda path: normalized_path_text(
            path
        ),
    )

    if not discovered:
        raise FileNotFoundError(
            "No NetCDF files were found under:\n"
            f"{RAW_DIR}"
        )

    absolute_files: list[
        Path
    ] = []

    anomaly_files: list[
        Path
    ] = []

    for file_path in discovered:
        if is_anomaly_path(
            file_path
        ):
            anomaly_files.append(
                file_path
            )
        else:
            absolute_files.append(
                file_path
            )

    if not absolute_files:
        raise FileNotFoundError(
            "No candidate absolute-SST files remain after anomaly exclusion."
        )

    return (
        absolute_files,
        anomaly_files,
    )


def calculate_sha256(
    file_path: Path,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Calculate a bounded-memory SHA-256 checksum."""

    digest = hashlib.sha256()

    with file_path.open(
        "rb"
    ) as file:
        while True:
            chunk = file.read(
                chunk_size
            )

            if not chunk:
                break

            digest.update(
                chunk
            )

    return digest.hexdigest()


# =============================================================================
# REFERENCE GRID AND MASK
# =============================================================================

def load_reference_grid_and_mask() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Load and validate the reference OISST grid and Pacific mask."""

    latitude = np.asarray(
        np.load(
            GRID_LAT_FILE
        ),
        dtype=np.float64,
    ).squeeze()

    longitude = np.asarray(
        np.load(
            GRID_LON_FILE
        ),
        dtype=np.float64,
    ).squeeze()

    mask_raw = np.asarray(
        np.load(
            PACIFIC_MASK_OISST_FILE
        )
    )

    if latitude.ndim != 1:
        raise ValueError(
            f"Latitude grid must be one-dimensional: {latitude.shape}"
        )

    if longitude.ndim != 1:
        raise ValueError(
            f"Longitude grid must be one-dimensional: {longitude.shape}"
        )

    expected_shape = (
        latitude.size,
        longitude.size,
    )

    if mask_raw.shape != expected_shape:
        raise ValueError(
            "Pacific mask shape does not match the reference grid.\n"
            f"Mask     : {mask_raw.shape}\n"
            f"Expected : {expected_shape}"
        )

    if not np.all(
        np.isfinite(
            latitude
        )
    ):
        raise ValueError(
            "Latitude grid contains nonfinite values."
        )

    if not np.all(
        np.isfinite(
            longitude
        )
    ):
        raise ValueError(
            "Longitude grid contains nonfinite values."
        )

    if not np.all(
        np.diff(
            latitude
        )
        > 0.0
    ):
        raise ValueError(
            "Reference latitude must be strictly increasing."
        )

    longitude_360 = np.mod(
        longitude,
        360.0,
    )

    if not np.all(
        np.diff(
            longitude_360
        )
        > 0.0
    ):
        raise ValueError(
            "Reference longitude must be strictly increasing in [0, 360)."
        )

    mask = (
        np.isfinite(
            mask_raw
        )
        & (
            mask_raw
            > 0
        )
    )

    if not np.any(
        mask
    ):
        raise ValueError(
            "Pacific mask contains no active cells."
        )

    return (
        latitude,
        longitude_360,
        mask,
    )


# =============================================================================
# SPHERICAL GEOMETRY
# =============================================================================

def coordinate_edges(
    centres: np.ndarray,
    periodic: bool,
) -> np.ndarray:
    """Calculate cell edges from regularly ordered cell centres."""

    centres = np.asarray(
        centres,
        dtype=np.float64,
    )

    if centres.ndim != 1:
        raise ValueError(
            "Coordinate centres must be one-dimensional."
        )

    if centres.size < 2:
        raise ValueError(
            "At least two coordinate centres are required."
        )

    midpoints = (
        centres[:-1]
        + centres[1:]
    ) / 2.0

    edges = np.empty(
        centres.size + 1,
        dtype=np.float64,
    )

    edges[1:-1] = midpoints

    if periodic:
        first_spacing = (
            centres[1]
            - centres[0]
        )
        last_spacing = (
            centres[-1]
            - centres[-2]
        )

        edges[0] = (
            centres[0]
            - first_spacing / 2.0
        )
        edges[-1] = (
            centres[-1]
            + last_spacing / 2.0
        )
    else:
        edges[0] = (
            centres[0]
            - (
                centres[1]
                - centres[0]
            ) / 2.0
        )
        edges[-1] = (
            centres[-1]
            + (
                centres[-1]
                - centres[-2]
            ) / 2.0
        )

    return edges


def build_spherical_geometry(
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> GridGeometry:
    """Precompute exact cell areas and unit vectors for the reference grid."""

    latitude_edges = coordinate_edges(
        latitude,
        periodic=False,
    )

    longitude_edges = coordinate_edges(
        longitude,
        periodic=True,
    )

    latitude_edges = np.clip(
        latitude_edges,
        -90.0,
        90.0,
    )

    delta_lambda = np.deg2rad(
        np.diff(
            longitude_edges
        )
    )

    latitude_factor = (
        np.sin(
            np.deg2rad(
                latitude_edges[1:]
            )
        )
        - np.sin(
            np.deg2rad(
                latitude_edges[:-1]
            )
        )
    )

    cell_area_km2 = (
        float(
            EARTH_RADIUS_KM
        )
        ** 2
        * latitude_factor[
            :,
            np.newaxis,
        ]
        * delta_lambda[
            np.newaxis,
            :,
        ]
    )

    if np.any(
        cell_area_km2
        <= 0.0
    ):
        raise ValueError(
            "Calculated cell areas must be strictly positive."
        )

    lon_grid, lat_grid = np.meshgrid(
        np.deg2rad(
            longitude
        ),
        np.deg2rad(
            latitude
        ),
    )

    cos_lat = np.cos(
        lat_grid
    )

    x_unit = (
        cos_lat
        * np.cos(
            lon_grid
        )
    )

    y_unit = (
        cos_lat
        * np.sin(
            lon_grid
        )
    )

    z_unit = np.sin(
        lat_grid
    )

    global_area_km2 = float(
        np.sum(
            cell_area_km2,
            dtype=np.float64,
        )
    )

    expected_global_area = (
        4.0
        * np.pi
        * float(
            EARTH_RADIUS_KM
        )
        ** 2
    )

    relative_error = abs(
        global_area_km2
        - expected_global_area
    ) / expected_global_area

    if relative_error > 1.0e-8:
        raise ValueError(
            "Global spherical area validation failed.\n"
            f"Calculated     : {global_area_km2:.6f} km²\n"
            f"Expected       : {expected_global_area:.6f} km²\n"
            f"Relative error : {relative_error:.3e}"
        )

    return GridGeometry(
        cell_area_km2=cell_area_km2,
        x_unit=x_unit,
        y_unit=y_unit,
        z_unit=z_unit,
        global_area_km2=global_area_km2,
    )


# =============================================================================
# NETCDF IDENTIFICATION
# =============================================================================

def first_existing_name(
    candidates: Iterable[
        str
    ],
    available_names: Iterable[
        str
    ],
    label: str,
) -> str:
    """Return the first candidate present in the available names."""

    available = set(
        available_names
    )

    for candidate in candidates:
        if candidate in available:
            return candidate

    raise KeyError(
        f"Unable to identify {label}. "
        f"Candidates: {tuple(candidates)}; "
        f"available: {sorted(available)}"
    )


def identify_sst_variable(
    dataset: xr.Dataset,
) -> str:
    """Identify the absolute SST variable."""

    for candidate in SST_VARIABLE_CANDIDATES:
        if candidate in dataset.data_vars:
            return candidate

    eligible: list[
        str
    ] = []

    for name, variable in dataset.data_vars.items():
        standard_name = str(
            variable.attrs.get(
                "standard_name",
                "",
            )
        ).lower()

        long_name = str(
            variable.attrs.get(
                "long_name",
                "",
            )
        ).lower()

        metadata = (
            standard_name
            + " "
            + long_name
        )

        if (
            "sea_surface_temperature"
            in metadata
            or "sea surface temperature"
            in metadata
        ):
            eligible.append(
                name
            )

    if len(
        eligible
    ) == 1:
        return eligible[0]

    raise KeyError(
        "Unable to identify a unique SST variable. "
        f"Data variables: {list(dataset.data_vars)}"
    )


def reject_anomaly_variable(
    dataset: xr.Dataset,
    sst_name: str,
    file_path: Path,
) -> None:
    """Reject variables whose metadata indicates an anomaly product."""

    variable = dataset[
        sst_name
    ]

    metadata_text = " ".join(
        str(
            variable.attrs.get(
                key,
                "",
            )
        ).lower()
        for key in (
            "long_name",
            "standard_name",
            "description",
            "comment",
            "title",
        )
    )

    if any(
        term in metadata_text
        for term in ABSOLUTE_SST_METADATA_REJECTION_TERMS
    ):
        raise ValueError(
            "The selected SST variable appears to be an anomaly product.\n"
            f"File     : {file_path}\n"
            f"Variable : {sst_name}"
        )


def standardize_sst(
    dataset: xr.Dataset,
    reference_latitude: np.ndarray,
    reference_longitude: np.ndarray,
) -> tuple[
    xr.DataArray,
    str,
    str,
    str,
    str,
]:
    """Identify dimensions and align SST to the reference grid."""

    sst_name = identify_sst_variable(
        dataset
    )

    sst = dataset[
        sst_name
    ]

    names = set(
        dataset.coords
    ) | set(
        dataset.dims
    )

    time_name = first_existing_name(
        TIME_COORDINATE_CANDIDATES,
        names,
        "time coordinate",
    )

    latitude_name = first_existing_name(
        LATITUDE_COORDINATE_CANDIDATES,
        names,
        "latitude coordinate",
    )

    longitude_name = first_existing_name(
        LONGITUDE_COORDINATE_CANDIDATES,
        names,
        "longitude coordinate",
    )

    required_dims = {
        time_name,
        latitude_name,
        longitude_name,
    }

    if not required_dims.issubset(
        set(
            sst.dims
        )
    ):
        raise ValueError(
            "SST variable does not contain the required dimensions.\n"
            f"SST dimensions : {sst.dims}\n"
            f"Required       : {sorted(required_dims)}"
        )

    extra_dims = [
        dim
        for dim in sst.dims
        if dim not in required_dims
    ]

    for dim in extra_dims:
        if sst.sizes[
            dim
        ] != 1:
            raise ValueError(
                "Unsupported non-singleton SST dimension.\n"
                f"Dimension : {dim}\n"
                f"Size      : {sst.sizes[dim]}"
            )

        sst = sst.isel(
            {
                dim: 0
            },
            drop=True,
        )

    sst = sst.transpose(
        time_name,
        latitude_name,
        longitude_name,
    )

    source_latitude = np.asarray(
        sst[
            latitude_name
        ].values,
        dtype=np.float64,
    )

    source_longitude = np.mod(
        np.asarray(
            sst[
                longitude_name
            ].values,
            dtype=np.float64,
        ),
        360.0,
    )

    if source_latitude.size != reference_latitude.size:
        raise ValueError(
            "Latitude dimension size differs from the reference grid."
        )

    if source_longitude.size != reference_longitude.size:
        raise ValueError(
            "Longitude dimension size differs from the reference grid."
        )

    if np.allclose(
        source_latitude,
        reference_latitude,
        rtol=0.0,
        atol=COORDINATE_TOLERANCE,
    ):
        pass
    elif np.allclose(
        source_latitude[::-1],
        reference_latitude,
        rtol=0.0,
        atol=COORDINATE_TOLERANCE,
    ):
        sst = sst.isel(
            {
                latitude_name: slice(
                    None,
                    None,
                    -1,
                )
            }
        )
    else:
        raise ValueError(
            "Source latitude values do not match the reference grid."
        )

    longitude_order = np.argsort(
        source_longitude
    )

    sorted_longitude = source_longitude[
        longitude_order
    ]

    if not np.allclose(
        sorted_longitude,
        reference_longitude,
        rtol=0.0,
        atol=COORDINATE_TOLERANCE,
    ):
        raise ValueError(
            "Source longitude values do not match the reference grid."
        )

    if not np.array_equal(
        longitude_order,
        np.arange(
            source_longitude.size
        ),
    ):
        sst = sst.isel(
            {
                longitude_name: longitude_order
            }
        )

    return (
        sst,
        sst_name,
        time_name,
        latitude_name,
        longitude_name,
    )


# =============================================================================
# SST UNITS AND PHYSICAL VALIDATION
# =============================================================================

def convert_sst_to_celsius(
    values: np.ndarray,
    units_attribute: object,
) -> tuple[
    np.ndarray,
    str,
]:
    """Convert one daily SST field to degrees Celsius."""

    values = np.asarray(
        values,
        dtype=np.float64,
    )

    units = str(
        units_attribute
        or ""
    ).strip().lower()

    kelvin_indicators = (
        "kelvin",
        "degree_k",
        "degrees_k",
        "degk",
        " k",
    )

    celsius_indicators = (
        "celsius",
        "degree_c",
        "degrees_c",
        "degc",
        "°c",
    )

    if (
        units == "k"
        or any(
            indicator in units
            for indicator in kelvin_indicators
        )
    ):
        converted = (
            values
            - 273.15
        )
        interpretation = "kelvin_to_celsius"

    elif (
        units in {
            "c",
            "deg c",
        }
        or any(
            indicator in units
            for indicator in celsius_indicators
        )
    ):
        converted = values
        interpretation = "celsius_from_metadata"

    else:
        finite_values = values[
            np.isfinite(
                values
            )
        ]

        if finite_values.size == 0:
            return (
                values,
                "undetermined_empty_field",
            )

        median_value = float(
            np.median(
                finite_values
            )
        )

        if median_value > 100.0:
            converted = (
                values
                - 273.15
            )
            interpretation = "kelvin_inferred_from_values"
        else:
            converted = values
            interpretation = "celsius_inferred_from_values"

    finite_converted = converted[
        np.isfinite(
            converted
        )
    ]

    if finite_converted.size:
        minimum = float(
            np.min(
                finite_converted
            )
        )

        maximum = float(
            np.max(
                finite_converted
            )
        )

        if (
            minimum
            < PLAUSIBLE_SST_MIN_C
            or maximum
            > PLAUSIBLE_SST_MAX_C
        ):
            raise ValueError(
                "Finite SST values lie outside the accepted physical range.\n"
                f"Minimum : {minimum:.6f} °C\n"
                f"Maximum : {maximum:.6f} °C"
            )

    return (
        converted,
        interpretation,
    )


# =============================================================================
# DAILY SCIENTIFIC CALCULATIONS
# =============================================================================

def area_weighted_mean(
    values: np.ndarray,
    selection: np.ndarray,
    weights: np.ndarray,
) -> float:
    """Calculate an area-weighted mean."""

    selected_weights = weights[
        selection
    ]

    denominator = float(
        np.sum(
            selected_weights,
            dtype=np.float64,
        )
    )

    if denominator <= 0.0:
        return np.nan

    numerator = float(
        np.sum(
            selected_weights
            * values[
                selection
            ],
            dtype=np.float64,
        )
    )

    return numerator / denominator


def area_weighted_standard_deviation(
    values: np.ndarray,
    selection: np.ndarray,
    weights: np.ndarray,
    mean_value: float,
) -> float:
    """Calculate an area-weighted population standard deviation."""

    if not np.isfinite(
        mean_value
    ):
        return np.nan

    selected_weights = weights[
        selection
    ]

    denominator = float(
        np.sum(
            selected_weights,
            dtype=np.float64,
        )
    )

    if denominator <= 0.0:
        return np.nan

    variance = float(
        np.sum(
            selected_weights
            * (
                values[
                    selection
                ]
                - mean_value
            )
            ** 2,
            dtype=np.float64,
        )
        / denominator
    )

    return float(
        np.sqrt(
            max(
                variance,
                0.0,
            )
        )
    )


def spherical_centroid(
    selection: np.ndarray,
    geometry: GridGeometry,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
]:
    """Calculate total area, spherical centroid, and concentration diagnostics."""

    if not np.any(
        selection
    ):
        return (
            np.nan,
            np.nan,
            0.0,
            np.nan,
            np.nan,
        )

    selected_area = geometry.cell_area_km2[
        selection
    ]

    total_area_km2 = float(
        np.sum(
            selected_area,
            dtype=np.float64,
        )
    )

    x_sum = float(
        np.sum(
            selected_area
            * geometry.x_unit[
                selection
            ],
            dtype=np.float64,
        )
    )

    y_sum = float(
        np.sum(
            selected_area
            * geometry.y_unit[
                selection
            ],
            dtype=np.float64,
        )
    )

    z_sum = float(
        np.sum(
            selected_area
            * geometry.z_unit[
                selection
            ],
            dtype=np.float64,
        )
    )

    resultant_norm = float(
        np.sqrt(
            x_sum
            ** 2
            + y_sum
            ** 2
            + z_sum
            ** 2
        )
    )

    numerical_tolerance = (
        np.finfo(
            np.float64
        ).eps
        * max(
            total_area_km2,
            1.0,
        )
    )

    if resultant_norm <= numerical_tolerance:
        raise ValueError(
            "The spherical centroid is numerically undefined because the "
            "area-weighted resultant vector is effectively zero."
        )

    longitude_360 = float(
        np.mod(
            np.rad2deg(
                np.arctan2(
                    y_sum,
                    x_sum,
                )
            ),
            360.0,
        )
    )

    latitude = float(
        np.rad2deg(
            np.arctan2(
                z_sum,
                np.hypot(
                    x_sum,
                    y_sum,
                ),
            )
        )
    )

    resultant_length = float(
        np.clip(
            resultant_norm
            / total_area_km2,
            0.0,
            1.0,
        )
    )

    angular_dispersion_deg = float(
        np.rad2deg(
            np.arccos(
                resultant_length
            )
        )
    )

    return (
        longitude_360,
        latitude,
        total_area_km2,
        resultant_length,
        angular_dispersion_deg,
    )


def longitude_360_to_180(
    longitude_360: float,
) -> float:
    """Convert [0, 360) longitude to [-180, 180)."""

    if not np.isfinite(
        longitude_360
    ):
        return np.nan

    return float(
        (
            longitude_360
            + 180.0
        )
        % 360.0
        - 180.0
    )


def process_daily_threshold(
    sst_c: np.ndarray,
    current_date: pd.Timestamp,
    source_file: Path,
    mask: np.ndarray,
    geometry: GridGeometry,
    threshold_c: float,
    threshold_source: str,
) -> dict[str, Any]:
    """Calculate all diagnostics for one date and one threshold."""

    if sst_c.shape != mask.shape:
        raise ValueError(
            "Daily SST shape differs from the Pacific mask.\n"
            f"SST  : {sst_c.shape}\n"
            f"Mask : {mask.shape}"
        )

    finite_sst = np.isfinite(
        sst_c
    )

    valid_ocean = (
        mask
        & finite_sst
    )

    pwp_selection = (
        valid_ocean
        & (
            sst_c
            >= threshold_c
        )
    )

    mask_cell_count = int(
        np.count_nonzero(
            mask
        )
    )

    valid_ocean_cell_count = int(
        np.count_nonzero(
            valid_ocean
        )
    )

    pwp_cell_count = int(
        np.count_nonzero(
            pwp_selection
        )
    )

    coverage_pct = (
        100.0
        * valid_ocean_cell_count
        / mask_cell_count
    )

    (
        longitude_360,
        latitude,
        area_km2,
        resultant_length,
        angular_dispersion_deg,
    ) = spherical_centroid(
        selection=pwp_selection,
        geometry=geometry,
    )

    mean_pwp_sst_c = area_weighted_mean(
        values=sst_c,
        selection=pwp_selection,
        weights=geometry.cell_area_km2,
    )

    weighted_sst_std_c = area_weighted_standard_deviation(
        values=sst_c,
        selection=pwp_selection,
        weights=geometry.cell_area_km2,
        mean_value=mean_pwp_sst_c,
    )

    max_pwp_sst_c = (
        float(
            np.max(
                sst_c[
                    pwp_selection
                ]
            )
        )
        if pwp_cell_count
        else np.nan
    )

    quality_flags: list[
        str
    ] = []

    if coverage_pct < MINIMUM_DAILY_COVERAGE_PCT:
        quality_flags.append(
            "LOW_COVERAGE"
        )

    if pwp_cell_count == 0:
        quality_flags.append(
            "NO_PWP_CELLS"
        )

    quality_flag = (
        "OK"
        if not quality_flags
        else "|".join(
            quality_flags
        )
    )

    return {
        "date": current_date.normalize(),
        "lon_360": longitude_360,
        "lon_180": longitude_360_to_180(
            longitude_360
        ),
        "lat": latitude,
        "area_km2": area_km2,
        "pwp_cell_count": pwp_cell_count,
        "valid_ocean_cell_count": valid_ocean_cell_count,
        "mask_cell_count": mask_cell_count,
        "ocean_data_coverage_pct": coverage_pct,
        "mean_pwp_sst_c": mean_pwp_sst_c,
        "max_pwp_sst_c": max_pwp_sst_c,
        "area_weighted_sst_std_c": weighted_sst_std_c,
        "threshold_c": threshold_c,
        "threshold_source": threshold_source,
        "threshold_fallback_used": False,
        "centroid_resultant_length": resultant_length,
        "centroid_angular_dispersion_deg": angular_dispersion_deg,
        "quality_flag": quality_flag,
        "source_file": str(
            source_file
        ),
    }


# =============================================================================
# NETCDF PROCESSING
# =============================================================================

def process_netcdf_file(
    file_path: Path,
    reference_latitude: np.ndarray,
    reference_longitude: np.ndarray,
    mask: np.ndarray,
    geometry: GridGeometry,
    contexts: tuple[
        ThresholdContext,
        ...,
    ],
) -> tuple[
    dict[
        float,
        list[
            dict[str, Any]
        ],
    ],
    InputFileMetadata,
]:
    """Process one NetCDF file for all selected thresholds."""

    results = {
        context.threshold_c: []
        for context in contexts
    }

    with xr.open_dataset(
        file_path,
        decode_cf=True,
        mask_and_scale=True,
        decode_times=True,
        cache=False,
    ) as dataset:

        (
            sst,
            sst_name,
            time_name,
            latitude_name,
            longitude_name,
        ) = standardize_sst(
            dataset=dataset,
            reference_latitude=reference_latitude,
            reference_longitude=reference_longitude,
        )

        reject_anomaly_variable(
            dataset=dataset,
            sst_name=sst_name,
            file_path=file_path,
        )

        units = sst.attrs.get(
            "units",
            "",
        )

        time_values = pd.to_datetime(
            sst[
                time_name
            ].values
        )

        if pd.isna(
            time_values
        ).any():
            raise ValueError(
                f"Invalid decoded time values in {file_path.name}"
            )

        if time_values.size != sst.sizes[
            time_name
        ]:
            raise ValueError(
                "Decoded time coordinate size differs from the SST time axis."
            )

        unit_interpretations: set[
            str
        ] = set()

        threshold_source = (
            "config.config::PWP_ALLOWED_THRESHOLDS_C"
            if RUN_ALL_THRESHOLDS
            else "config.config::PWP_SST_THRESHOLD_C"
        )

        for time_index, time_value in enumerate(
            time_values
        ):
            daily_values = np.asarray(
                sst.isel(
                    {
                        time_name: time_index
                    }
                ).values,
                dtype=np.float64,
            )

            (
                daily_sst_c,
                unit_interpretation,
            ) = convert_sst_to_celsius(
                values=daily_values,
                units_attribute=units,
            )

            unit_interpretations.add(
                unit_interpretation
            )

            current_date = pd.Timestamp(
                time_value
            )

            for context in contexts:
                results[
                    context.threshold_c
                ].append(
                    process_daily_threshold(
                        sst_c=daily_sst_c,
                        current_date=current_date,
                        source_file=file_path,
                        mask=mask,
                        geometry=geometry,
                        threshold_c=context.threshold_c,
                        threshold_source=threshold_source,
                    )
                )

        metadata = InputFileMetadata(
            source_file=str(
                file_path
            ),
            sst_variable=sst_name,
            time_coordinate=time_name,
            latitude_coordinate=latitude_name,
            longitude_coordinate=longitude_name,
            original_units=str(
                units
            ),
            unit_interpretation=";".join(
                sorted(
                    unit_interpretations
                )
            ),
            records_per_threshold=int(
                time_values.size
            ),
            thresholds_processed_c=tuple(
                context.threshold_c
                for context in contexts
            ),
        )

    return (
        results,
        metadata,
    )


# =============================================================================
# RESULT VALIDATION
# =============================================================================

def build_validated_dataframe(
    records: list[
        dict[str, Any]
    ],
    expected_threshold_c: float,
) -> pd.DataFrame:
    """Build, sort, and rigorously validate one threshold time series."""

    if not records:
        raise ValueError(
            "No daily records were generated."
        )

    dataframe = pd.DataFrame(
        records
    )

    required_columns = (
        "date",
        "lon_360",
        "lon_180",
        "lat",
        "area_km2",
        "pwp_cell_count",
        "valid_ocean_cell_count",
        "mask_cell_count",
        "ocean_data_coverage_pct",
        "mean_pwp_sst_c",
        "max_pwp_sst_c",
        "area_weighted_sst_std_c",
        "threshold_c",
        "threshold_source",
        "threshold_fallback_used",
        "centroid_resultant_length",
        "centroid_angular_dispersion_deg",
        "quality_flag",
        "source_file",
    )

    missing_columns = [
        column
        for column in required_columns
        if column not in dataframe.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Generated output is missing columns: {missing_columns}"
        )

    dataframe[
        "date"
    ] = pd.to_datetime(
        dataframe[
            "date"
        ],
        errors="raise",
    )

    dataframe = dataframe.sort_values(
        "date"
    ).reset_index(
        drop=True
    )

    duplicated = dataframe[
        "date"
    ].duplicated(
        keep=False
    )

    if duplicated.any():
        duplicate_dates = dataframe.loc[
            duplicated,
            "date",
        ].dt.strftime(
            "%Y-%m-%d"
        ).tolist()

        raise ValueError(
            "Duplicate dates were generated:\n"
            + "\n".join(
                duplicate_dates[
                    :20
                ]
            )
        )

    threshold_values = pd.to_numeric(
        dataframe[
            "threshold_c"
        ],
        errors="coerce",
    ).dropna().unique()

    if threshold_values.size != 1:
        raise ValueError(
            "A threshold-specific output contains multiple threshold values."
        )

    if not np.isclose(
        float(
            threshold_values[0]
        ),
        expected_threshold_c,
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError(
            "Stored threshold differs from the expected threshold."
        )

    if not dataframe[
        "ocean_data_coverage_pct"
    ].between(
        0.0,
        100.0,
        inclusive="both",
    ).all():
        raise ValueError(
            "Ocean-data coverage must lie in [0, 100]."
        )

    valid_resultant = dataframe[
        "centroid_resultant_length"
    ].dropna()

    if not valid_resultant.between(
        0.0,
        1.0,
        inclusive="both",
    ).all():
        raise ValueError(
            "Centroid resultant length must lie in [0, 1]."
        )

    if (
        dataframe[
            "area_km2"
        ]
        < 0.0
    ).any():
        raise ValueError(
            "PWP area cannot be negative."
        )

    complete_dates = pd.date_range(
        start=dataframe[
            "date"
        ].min(),
        end=dataframe[
            "date"
        ].max(),
        freq="D",
    )

    missing_dates = complete_dates.difference(
        pd.DatetimeIndex(
            dataframe[
                "date"
            ]
        )
    )

    dataframe.attrs[
        "missing_dates"
    ] = missing_dates

    return dataframe


# =============================================================================
# STATISTICAL UTILITIES
# =============================================================================

def circular_mean_degrees(
    values: np.ndarray,
) -> float:
    """Calculate a circular mean in [0, 360)."""

    values = np.asarray(
        values,
        dtype=np.float64,
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if values.size == 0:
        return np.nan

    radians = np.deg2rad(
        values
    )

    return float(
        np.mod(
            np.rad2deg(
                np.arctan2(
                    np.mean(
                        np.sin(
                            radians
                        )
                    ),
                    np.mean(
                        np.cos(
                            radians
                        )
                    ),
                )
            ),
            360.0,
        )
    )


def dataframe_summary(
    dataframe: pd.DataFrame,
) -> dict[str, Any]:
    """Build one JSON-safe summary of a threshold time series."""

    missing_dates = dataframe.attrs.get(
        "missing_dates",
        pd.DatetimeIndex([]),
    )

    quality_counts = {
        str(
            key
        ): int(
            value
        )
        for key, value in dataframe[
            "quality_flag"
        ].value_counts(
            dropna=False
        ).items()
    }

    return {
        "records": int(
            len(
                dataframe
            )
        ),
        "start_date": dataframe[
            "date"
        ].min().strftime(
            "%Y-%m-%d"
        ),
        "end_date": dataframe[
            "date"
        ].max().strftime(
            "%Y-%m-%d"
        ),
        "missing_calendar_dates": [
            date.strftime(
                "%Y-%m-%d"
            )
            for date in missing_dates
        ],
        "quality_flag_counts": quality_counts,
        "mean_area_km2": float(
            dataframe[
                "area_km2"
            ].mean()
        ),
        "median_area_km2": float(
            dataframe[
                "area_km2"
            ].median()
        ),
        "minimum_area_km2": float(
            dataframe[
                "area_km2"
            ].min()
        ),
        "maximum_area_km2": float(
            dataframe[
                "area_km2"
            ].max()
        ),
        "mean_centroid_longitude_360": circular_mean_degrees(
            dataframe[
                "lon_360"
            ].to_numpy(
                dtype=np.float64
            )
        ),
        "mean_centroid_latitude": float(
            dataframe[
                "lat"
            ].mean()
        ),
        "mean_internal_sst_c": float(
            dataframe[
                "mean_pwp_sst_c"
            ].mean()
        ),
        "mean_weighted_sst_std_c": float(
            dataframe[
                "area_weighted_sst_std_c"
            ].mean()
        ),
        "mean_resultant_length": float(
            dataframe[
                "centroid_resultant_length"
            ].mean()
        ),
        "mean_angular_dispersion_deg": float(
            dataframe[
                "centroid_angular_dispersion_deg"
            ].mean()
        ),
    }


# =============================================================================
# EXPORTS AND REPORTS
# =============================================================================

def export_dataframe(
    dataframe: pd.DataFrame,
    output_file: Path,
) -> None:
    """Export one threshold time series using stable formatting."""

    export_frame = dataframe.copy()

    export_frame[
        "date"
    ] = export_frame[
        "date"
    ].dt.strftime(
        "%Y-%m-%d"
    )

    export_frame.to_csv(
        output_file,
        index=False,
        encoding="utf-8",
        float_format="%.10f",
        lineterminator="\n",
    )


def build_metadata(
    context: ThresholdContext,
    dataframe: pd.DataFrame,
    input_files: list[Path],
    excluded_anomaly_files: list[Path],
    file_metadata: list[
        InputFileMetadata
    ],
    latitude: np.ndarray,
    longitude: np.ndarray,
    mask: np.ndarray,
    geometry: GridGeometry,
    mask_checksum: str,
    contexts: tuple[
        ThresholdContext,
        ...,
    ],
) -> dict[str, Any]:
    """Build machine-readable scientific provenance."""

    threshold_source = (
        "config.config::PWP_ALLOWED_THRESHOLDS_C"
        if RUN_ALL_THRESHOLDS
        else "config.config::PWP_SST_THRESHOLD_C"
    )

    return {
        "program": {
            "name": PROGRAM_NAME,
            "version": PROGRAM_VERSION,
            "script": str(
                SCRIPT_FILE
            ),
            "generated_utc": datetime.now(
                timezone.utc
            ).isoformat(),
        },
        "project": {
            "project_root": str(
                PROJECT_DIR
            ),
            "raw_directory": str(
                RAW_DIR
            ),
            "pipeline_module": "centroid",
            "module_processed_directory": str(
                context.centroid_csv.parent
            ),
            "module_report_directory": str(
                context.technical_report.parent
            ),
            "output_csv": str(
                context.centroid_csv
            ),
        },
        "execution": {
            "mode": (
                "all thresholds"
                if RUN_ALL_THRESHOLDS
                else "single threshold"
            ),
            "thresholds_selected_c": [
                item.threshold_c
                for item in contexts
            ],
        },
        "threshold": {
            "value_c": context.threshold_c,
            "folder_name": context.folder_name,
            "source": threshold_source,
            "fallback_used": False,
            "fallback_policy": (
                "No silent fallback. config/config.py is mandatory."
            ),
            "definition": (
                "finite SST cells inside the validated Pacific mask "
                f"with SST >= {context.threshold_c:.3f} °C"
            ),
        },
        "geometry": {
            "earth_radius_km": float(
                EARTH_RADIUS_KM
            ),
            "cell_area_method": (
                "exact spherical latitude-longitude cell area from cell edges"
            ),
            "centroid_method": (
                "exact-area-weighted three-dimensional Cartesian unit-vector "
                "resultant on a sphere"
            ),
            "global_grid_area_km2": geometry.global_area_km2,
            "grid_shape": [
                int(
                    latitude.size
                ),
                int(
                    longitude.size
                ),
            ],
            "latitude_min_deg": float(
                latitude.min()
            ),
            "latitude_max_deg": float(
                latitude.max()
            ),
            "longitude_min_deg_e": float(
                longitude.min()
            ),
            "longitude_max_deg_e": float(
                longitude.max()
            ),
            "pacific_mask_cells": int(
                np.count_nonzero(
                    mask
                )
            ),
            "pacific_mask_area_km2": float(
                np.sum(
                    geometry.cell_area_km2[
                        mask
                    ],
                    dtype=np.float64,
                )
            ),
            "pacific_mask_sha256": mask_checksum,
        },
        "processing": {
            "absolute_sst_files": len(
                input_files
            ),
            "anomaly_files_excluded": len(
                excluded_anomaly_files
            ),
            "minimum_daily_coverage_pct": (
                MINIMUM_DAILY_COVERAGE_PCT
            ),
            "plausible_sst_range_c": [
                PLAUSIBLE_SST_MIN_C,
                PLAUSIBLE_SST_MAX_C,
            ],
            "coordinate_tolerance": COORDINATE_TOLERANCE,
            "input_files": [
                asdict(
                    item
                )
                for item in file_metadata
            ],
        },
        "output_series": dataframe_summary(
            dataframe
        ),
        "output_columns": list(
            dataframe.columns
        ),
    }


def write_json(
    metadata: dict[str, Any],
    output_file: Path,
) -> None:
    """Write UTF-8 JSON with deterministic indentation."""

    output_file.write_text(
        json.dumps(
            metadata,
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_technical_report(
    context: ThresholdContext,
    dataframe: pd.DataFrame,
    metadata: dict[str, Any],
    input_files: list[Path],
    excluded_anomaly_files: list[Path],
) -> None:
    """Write the detailed processing report."""

    summary = dataframe_summary(
        dataframe
    )

    lines = [
        "PACIFIC WARM POOL CENTROID PROCESSING REPORT",
        "=" * 78,
        "",
        "1. PROGRAM",
        "-" * 78,
        f"Program                    : {PROGRAM_NAME}",
        f"Version                    : {PROGRAM_VERSION}",
        f"Script                     : {SCRIPT_FILE}",
        f"Generated UTC              : {datetime.now(timezone.utc).isoformat()}",
        "",
        "2. CONFIGURATION",
        "-" * 78,
        f"Project root               : {PROJECT_DIR}",
        (
            "Execution mode             : "
            f"{'all thresholds' if RUN_ALL_THRESHOLDS else 'single threshold'}"
        ),
        f"Threshold                  : {context.threshold_c:.3f} °C",
        f"Threshold folder           : {context.folder_name}",
        f"Threshold fallback used    : False",
        f"Earth radius               : {float(EARTH_RADIUS_KM):.4f} km",
        "",
        "3. INPUTS",
        "-" * 78,
        f"Raw-data directory         : {RAW_DIR}",
        f"Absolute-SST files         : {len(input_files):,}",
        f"Anomaly files excluded     : {len(excluded_anomaly_files):,}",
        f"Pacific mask               : {PACIFIC_MASK_OISST_FILE}",
        f"Latitude grid              : {GRID_LAT_FILE}",
        f"Longitude grid             : {GRID_LON_FILE}",
        "",
        "4. OUTPUT SERIES",
        "-" * 78,
        f"Output CSV                 : {context.centroid_csv}",
        f"Records                    : {summary['records']:,}",
        f"Start date                 : {summary['start_date']}",
        f"End date                   : {summary['end_date']}",
        (
            "Missing calendar dates     : "
            f"{len(summary['missing_calendar_dates']):,}"
        ),
        f"Mean area                  : {summary['mean_area_km2']:,.3f} km²",
        f"Median area                : {summary['median_area_km2']:,.3f} km²",
        f"Minimum area               : {summary['minimum_area_km2']:,.3f} km²",
        f"Maximum area               : {summary['maximum_area_km2']:,.3f} km²",
        (
            "Mean centroid longitude    : "
            f"{summary['mean_centroid_longitude_360']:.8f} °E"
        ),
        (
            "Mean centroid latitude     : "
            f"{summary['mean_centroid_latitude']:.8f}°"
        ),
        (
            "Mean internal SST          : "
            f"{summary['mean_internal_sst_c']:.8f} °C"
        ),
        (
            "Mean weighted SST SD       : "
            f"{summary['mean_weighted_sst_std_c']:.8f} °C"
        ),
        (
            "Mean resultant length      : "
            f"{summary['mean_resultant_length']:.10f}"
        ),
        (
            "Mean angular dispersion    : "
            f"{summary['mean_angular_dispersion_deg']:.8f}°"
        ),
        "",
        "5. QUALITY FLAGS",
        "-" * 78,
    ]

    for flag, count in summary[
        "quality_flag_counts"
    ].items():
        lines.append(
            f"{flag:<27s}: {count:,}"
        )

    lines.extend(
        [
            "",
            "6. METHOD",
            "-" * 78,
            (
                "Grid-cell areas were calculated exactly on a sphere from "
                "latitude and longitude cell edges."
            ),
            (
                "The centroid was calculated from exact-area-weighted "
                "three-dimensional Cartesian unit vectors and transformed "
                "back to geographic longitude and latitude."
            ),
            (
                "The normalized resultant-vector magnitude and corresponding "
                "angular dispersion were stored as geometric concentration "
                "diagnostics; they are not uncertainty estimates."
            ),
            "",
            "7. OUTPUTS",
            "-" * 78,
            str(
                context.centroid_csv
            ),
            str(
                context.technical_report
            ),
            str(
                context.scientific_report
            ),
            str(
                context.metadata_json
            ),
            "",
            "8. STATUS",
            "-" * 78,
            "Status                     : SUCCESS",
            "",
            "=" * 78,
            "END OF REPORT",
            "=" * 78,
        ]
    )

    context.technical_report.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )


def write_scientific_report(
    context: ThresholdContext,
    dataframe: pd.DataFrame,
) -> None:
    """Write a technical-scientific methodology record."""

    summary = dataframe_summary(
        dataframe
    )

    lines = [
        "PACIFIC WARM POOL CENTROID — TECHNICAL-SCIENTIFIC METHOD",
        "=" * 78,
        "",
        "1. OPERATIONAL DEFINITION",
        "-" * 78,
        (
            "For each day, the Pacific Warm Pool was defined as all finite "
            "sea-surface-temperature grid cells inside the validated Pacific "
            f"Ocean mask satisfying SST >= {context.threshold_c:.3f} °C."
        ),
        "",
        "2. GRID-CELL AREA",
        "-" * 78,
        (
            "The area of each regular latitude-longitude cell was evaluated "
            "on a sphere from its angular bounds:"
        ),
        "A_i = R² Δλ_i [sin(phi_north,i) - sin(phi_south,i)].",
        (
            f"The adopted mean Earth radius was {float(EARTH_RADIUS_KM):.4f} km."
        ),
        "",
        "3. SPHERICAL CENTROID",
        "-" * 78,
        (
            "Each selected cell centre was represented by a Cartesian unit "
            "vector. These vectors were weighted by exact grid-cell area, "
            "summed, and transformed back to longitude and latitude."
        ),
        "X = Σ A_i cos(phi_i) cos(lambda_i)",
        "Y = Σ A_i cos(phi_i) sin(lambda_i)",
        "Z = Σ A_i sin(phi_i)",
        "lambda_c = atan2(Y, X)",
        "phi_c = atan2(Z, sqrt(X² + Y²))",
        (
            "This formulation treats longitude as circular and therefore "
            "avoids discontinuity errors at 0°/360°."
        ),
        "",
        "4. INTERNAL SST STATISTICS",
        "-" * 78,
        (
            "The mean SST and spatial SST standard deviation inside the PWP "
            "were weighted by exact cell area."
        ),
        "",
        "5. GEOMETRIC CONCENTRATION",
        "-" * 78,
        (
            "The normalized magnitude of the weighted Cartesian resultant "
            "was stored as centroid_resultant_length."
        ),
        (
            "centroid_angular_dispersion_deg was calculated as "
            "arccos(resultant_length) in degrees."
        ),
        (
            "These quantities describe spatial concentration and are not "
            "formal confidence intervals."
        ),
        "",
        "6. QUALITY CONTROL",
        "-" * 78,
        (
            "Daily valid-ocean coverage was calculated relative to the number "
            "of active cells in the Pacific mask."
        ),
        (
            "Records below "
            f"{MINIMUM_DAILY_COVERAGE_PCT:.3f}% coverage were flagged."
        ),
        (
            "Absolute-SST inputs, coordinate consistency, unit conversion, "
            "plausible SST limits, duplicate dates, and threshold consistency "
            "were validated."
        ),
        "",
        "7. SERIES SUMMARY",
        "-" * 78,
        f"Records                    : {summary['records']:,}",
        f"Start date                 : {summary['start_date']}",
        f"End date                   : {summary['end_date']}",
        f"Mean PWP area              : {summary['mean_area_km2']:,.3f} km²",
        (
            "Mean centroid longitude    : "
            f"{summary['mean_centroid_longitude_360']:.8f} °E"
        ),
        (
            "Mean centroid latitude     : "
            f"{summary['mean_centroid_latitude']:.8f}°"
        ),
        "",
        "8. REPRODUCIBILITY",
        "-" * 78,
        f"Daily series               : {context.centroid_csv}",
        f"Machine-readable metadata  : {context.metadata_json}",
        "",
        "=" * 78,
        "END OF REPORT",
        "=" * 78,
    ]

    context.scientific_report.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )


# =============================================================================
# CROSS-THRESHOLD COMPARISON
# =============================================================================

def build_comparison_summary(
    dataframes: dict[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    """Build one compact comparison table across thresholds."""

    records: list[
        dict[str, Any]
    ] = []

    for threshold_c, dataframe in sorted(
        dataframes.items()
    ):
        summary = dataframe_summary(
            dataframe
        )

        records.append(
            {
                "threshold_c": threshold_c,
                "threshold_folder": threshold_folder_name(
                    threshold_c
                ),
                "records": summary[
                    "records"
                ],
                "start_date": summary[
                    "start_date"
                ],
                "end_date": summary[
                    "end_date"
                ],
                "mean_area_km2": summary[
                    "mean_area_km2"
                ],
                "median_area_km2": summary[
                    "median_area_km2"
                ],
                "minimum_area_km2": summary[
                    "minimum_area_km2"
                ],
                "maximum_area_km2": summary[
                    "maximum_area_km2"
                ],
                "mean_centroid_longitude_360": summary[
                    "mean_centroid_longitude_360"
                ],
                "mean_centroid_latitude": summary[
                    "mean_centroid_latitude"
                ],
                "mean_internal_sst_c": summary[
                    "mean_internal_sst_c"
                ],
                "mean_weighted_sst_std_c": summary[
                    "mean_weighted_sst_std_c"
                ],
                "mean_resultant_length": summary[
                    "mean_resultant_length"
                ],
                "mean_angular_dispersion_deg": summary[
                    "mean_angular_dispersion_deg"
                ],
            }
        )

    return pd.DataFrame(
        records
    ).sort_values(
        "threshold_c"
    ).reset_index(
        drop=True
    )


def export_comparison(
    dataframes: dict[
        float,
        pd.DataFrame,
    ],
) -> tuple[
    Path,
    Path,
] | None:
    """Export cross-threshold summary products."""

    if len(
        dataframes
    ) <= 1:
        return None

    comparison = build_comparison_summary(
        dataframes
    )

    comparison.to_csv(
        OUTPUT_COMPARISON_TABLE,
        index=False,
        encoding="utf-8",
        float_format="%.10f",
        lineterminator="\n",
    )

    lines = [
        "PACIFIC WARM POOL THRESHOLD COMPARISON — CENTROID AND AREA",
        "=" * 78,
        "",
        f"Generated UTC              : {datetime.now(timezone.utc).isoformat()}",
        "",
        "SUMMARY",
        "-" * 78,
    ]

    for record in comparison.to_dict(
        orient="records"
    ):
        lines.extend(
            [
                f"Threshold                  : {record['threshold_c']:.1f} °C",
                f"Folder                     : {record['threshold_folder']}",
                f"Records                    : {record['records']:,}",
                f"Mean area                  : {record['mean_area_km2']:,.3f} km²",
                f"Median area                : {record['median_area_km2']:,.3f} km²",
                (
                    "Mean centroid longitude    : "
                    f"{record['mean_centroid_longitude_360']:.8f} °E"
                ),
                (
                    "Mean centroid latitude     : "
                    f"{record['mean_centroid_latitude']:.8f}°"
                ),
                f"Mean internal SST          : {record['mean_internal_sst_c']:.8f} °C",
                "",
            ]
        )

    lines.extend(
        [
            "FILES",
            "-" * 78,
            str(
                OUTPUT_COMPARISON_TABLE
            ),
            str(
                OUTPUT_COMPARISON_REPORT
            ),
            "",
            "STATUS",
            "-" * 78,
            "Status                     : SUCCESS",
            "",
            "=" * 78,
            "END OF REPORT",
            "=" * 78,
        ]
    )

    OUTPUT_COMPARISON_REPORT.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )

    return (
        OUTPUT_COMPARISON_TABLE,
        OUTPUT_COMPARISON_REPORT,
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    """Run the complete centroid and area processing workflow."""

    print_header(
        PROGRAM_NAME
    )

    contexts = selected_threshold_contexts()

    validate_configuration(
        contexts
    )

    validate_required_inputs_and_outputs(
        contexts
    )

    print_section(
        "CONFIGURATION, INPUTS, AND OUTPUTS"
    )

    print_key_value(
        "Program version",
        PROGRAM_VERSION,
    )

    print_key_value(
        "Project root",
        PROJECT_DIR,
    )

    print_key_value(
        "Script",
        SCRIPT_FILE,
    )

    print_key_value(
        "Raw OISST directory",
        RAW_DIR,
    )

    print_key_value(
        "Pacific mask",
        PACIFIC_MASK_OISST_FILE,
    )

    print_key_value(
        "Latitude grid",
        GRID_LAT_FILE,
    )

    print_key_value(
        "Longitude grid",
        GRID_LON_FILE,
    )

    print_key_value(
        "Execution mode",
        (
            "all thresholds"
            if RUN_ALL_THRESHOLDS
            else "single threshold"
        ),
    )

    print_key_value(
        "Thresholds selected",
        ", ".join(
            f"{context.threshold_c:.1f} °C"
            for context in contexts
        ),
    )

    for context in contexts:
        print_key_value(
            f"Centroid data [{context.threshold_c:.1f} °C]",
            context.centroid_csv,
        )

        print_key_value(
            f"Centroid reports [{context.threshold_c:.1f} °C]",
            context.technical_report.parent,
        )

    latitude, longitude, mask = (
        load_reference_grid_and_mask()
    )

    geometry = build_spherical_geometry(
        latitude=latitude,
        longitude=longitude,
    )

    input_files, anomaly_files = (
        discover_netcdf_files()
    )

    mask_checksum = calculate_sha256(
        PACIFIC_MASK_OISST_FILE
    )

    print_section(
        "REFERENCE GRID AND MASK"
    )

    print_key_value(
        "Grid shape",
        (
            latitude.size,
            longitude.size,
        ),
    )

    print_key_value(
        "Latitude range",
        (
            f"{latitude.min():.3f} "
            f"to {latitude.max():.3f}°"
        ),
    )

    print_key_value(
        "Longitude range",
        (
            f"{longitude.min():.3f} "
            f"to {longitude.max():.3f}°E"
        ),
    )

    print_key_value(
        "Pacific-mask cells",
        f"{np.count_nonzero(mask):,}",
    )

    print_key_value(
        "Global grid area",
        f"{geometry.global_area_km2:,.3f} km²",
    )

    print_key_value(
        "Absolute-SST files",
        f"{len(input_files):,}",
    )

    print_key_value(
        "Anomaly files excluded",
        f"{len(anomaly_files):,}",
    )

    all_records: dict[
        float,
        list[
            dict[str, Any]
        ],
    ] = {
        context.threshold_c: []
        for context in contexts
    }

    file_metadata: list[
        InputFileMetadata
    ] = []

    print_section(
        "PROCESSING ABSOLUTE OISST FILES"
    )

    print(
        "Each daily OISST field is read once and evaluated for all selected "
        "thresholds."
    )

    for index, file_path in enumerate(
        input_files,
        start=1,
    ):
        try:
            display_path = file_path.relative_to(
                RAW_DIR
            )
        except ValueError:
            display_path = file_path

        print(
            f"[{index:>5}/{len(input_files):<5}] "
            f"{display_path}"
        )

        try:
            records_by_threshold, metadata = process_netcdf_file(
                file_path=file_path,
                reference_latitude=latitude,
                reference_longitude=longitude,
                mask=mask,
                geometry=geometry,
                contexts=contexts,
            )
        except Exception as error:
            raise RuntimeError(
                "Processing failed for NetCDF file:\n"
                f"{file_path}\n\n"
                f"Original error: {error}"
            ) from error

        for threshold_c, records in records_by_threshold.items():
            all_records[
                threshold_c
            ].extend(
                records
            )

        file_metadata.append(
            metadata
        )

    dataframes: dict[
        float,
        pd.DataFrame,
    ] = {}

    created_files: list[
        Path
    ] = []

    for context in contexts:
        print_section(
            f"VALIDATING AND EXPORTING — {context.threshold_c:.1f} °C"
        )

        dataframe = build_validated_dataframe(
            records=all_records[
                context.threshold_c
            ],
            expected_threshold_c=context.threshold_c,
        )

        export_dataframe(
            dataframe=dataframe,
            output_file=context.centroid_csv,
        )

        metadata = build_metadata(
            context=context,
            dataframe=dataframe,
            input_files=input_files,
            excluded_anomaly_files=anomaly_files,
            file_metadata=file_metadata,
            latitude=latitude,
            longitude=longitude,
            mask=mask,
            geometry=geometry,
            mask_checksum=mask_checksum,
            contexts=contexts,
        )

        write_json(
            metadata=metadata,
            output_file=context.metadata_json,
        )

        write_technical_report(
            context=context,
            dataframe=dataframe,
            metadata=metadata,
            input_files=input_files,
            excluded_anomaly_files=anomaly_files,
        )

        write_scientific_report(
            context=context,
            dataframe=dataframe,
        )

        dataframes[
            context.threshold_c
        ] = dataframe

        created_files.extend(
            [
                context.centroid_csv,
                context.technical_report,
                context.scientific_report,
                context.metadata_json,
            ]
        )

        summary = dataframe_summary(
            dataframe
        )

        print_key_value(
            "Daily records",
            f"{summary['records']:,}",
        )

        print_key_value(
            "Start date",
            summary[
                "start_date"
            ],
        )

        print_key_value(
            "End date",
            summary[
                "end_date"
            ],
        )

        print_key_value(
            "Missing calendar dates",
            f"{len(summary['missing_calendar_dates']):,}",
        )

        print_key_value(
            "Mean PWP area",
            f"{summary['mean_area_km2']:,.3f} km²",
        )

        print_key_value(
            "Mean centroid longitude",
            (
                f"{summary['mean_centroid_longitude_360']:.8f} °E"
            ),
        )

        print_key_value(
            "Mean centroid latitude",
            (
                f"{summary['mean_centroid_latitude']:.8f}°"
            ),
        )

    comparison_outputs = export_comparison(
        dataframes
    )

    if comparison_outputs is not None:
        created_files.extend(
            comparison_outputs
        )

    print_section(
        "FILES CREATED"
    )

    for path in created_files:
        print(
            path
        )

    print()

    print_rule()

    print(
        "PROGRAM 05 COMPLETED SUCCESSFULLY."
    )

    print(
        "Daily PWP centroid and area series were generated for: "
        + ", ".join(
            f"{context.threshold_c:.1f} °C"
            for context in contexts
        )
    )

    print(
        "Each daily OISST field was read once for all selected thresholds."
    )

    print_rule()


if __name__ == "__main__":
    main()