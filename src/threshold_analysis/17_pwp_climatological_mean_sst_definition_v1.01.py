#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
===============================================================================
PROJECT
    Ocean Spectral Analysis Framework (OSAF)
    Pacific Warm Pool (PWP) Scientific Analysis Pipeline

PROGRAM
    17_pwp_climatological_mean_sst_definition_v1.01.py

TITLE
    Long-Term Mean SST Background, Pacific Mask, and Thermal PWP Definitions

VERSION
    1.01

PURPOSE
    Create a publication-quality five-panel (A-E) climatological methodological
    figure showing the long-term SST background within which the 28.0, 28.5,
    and 29.0 °C Pacific Warm Pool definitions exist.

    This is a CLIMATOLOGICAL BACKGROUND figure. It is deliberately distinct
    from the planned grid-cell occurrence/persistence analysis.

    The five panels are:

        A) Long-term mean daily OISST field
        B) Fixed Pacific Ocean mask
        C) Long-term mean SST restricted to the Pacific mask
        D) Mean Pacific SST with 28.0, 28.5, and 29.0 °C isotherms
        E) Nested thermal domains obtained by applying the three thresholds
           to the long-term mean SST field, with the long-term mean spherical
           centroids from the frozen daily centroid products marked.

IMPORTANT SCIENTIFIC DISTINCTION
    Panel E does NOT show climatological PWP occurrence frequency.

    In this program:

        mean_SST_i = mean_t [ SST_i(t) ]

    and the Panel-E background domain is diagnosed from:

        mean_SST_i >= T

    This answers:

        "What is the long-term mean SST structure within which the three
         thermal definitions exist?"

    It does NOT answer:

        "On what fraction of valid days does grid cell i satisfy SST >= T?"

    The latter requires:

        F_i(T) =
            sum_t 1[SST_i(t) >= T]
            --------------------------------
            sum_t 1[SST_i(t) is finite]

    and must be calculated directly from the DAILY threshold condition in a
    separate occurrence/persistence analysis. The two quantities must never
    be substituted for one another.

SCIENTIFIC ROLE IN PAPER 1
    Proposed visual sequence:

        Methods Figure
            one-day PWP construction

        Results Figure 1a / contextual figure
            this long-term mean SST background A-E figure

        Results Figure 1b / occurrence figure
            climatological fraction of valid days satisfying SST >= T

        Results Figure 2
            distribution of all daily spherical centroid positions

    Together these distinguish:
        mean temperature
        threshold occurrence
        threshold-defined thermal domain
        centroid position

AUTHORITATIVE TEMPORAL RECORD
    1981-09-01 to 2026-07-29

    N = 16,403 target daily observations.

    The exact target dates are read from the canonical frozen Program-05
    centroid product for 28.0 °C:

        data/processed/28/centroid/pwp_centroid_series.csv

    The 28.5 and 29.0 °C canonical centroid series must contain the same date
    sequence.

SST INPUT
    NOAA Optimum Interpolation Sea Surface Temperature (OISST)

    Daily absolute SST
    Native 0.25° latitude × 0.25° longitude grid

    Raw source directory:
        data/raw/

    The program prefers yearly absolute-SST files named similarly to:

        sst.day.mean.1981.nc
        ...
        sst.day.mean.2026.nc

    and excludes probable SST-anomaly products.

PACIFIC MASK AND GRID INPUTS
    Reused from canonical Program 17 configuration:

        data/processed/pacific_mask_oisst.npy
        data/processed/grid_lat.npy
        data/processed/grid_lon.npy

UPSTREAM PROGRAMS
    Program 01 / mask construction
    Program 04 / mask-grid alignment, where applicable
    Program 05 / frozen daily spherical PWP centroid products
    Program 06 / final QC policy
    Program 17 / canonical methodological-domain utilities
    Program 34 v1.01 / long-term mean spherical centroid summary

CENTROID PROVENANCE
    Mean centroid markers are NOT calculated from the thresholded mean SST
    field.

    They are the long-term mean spherical centroids of the canonical DAILY
    centroid populations. The daily lon/lat positions are converted to 3-D
    unit-sphere vectors, averaged, and converted back to spherical coordinates.

    Frozen values expected from Program 34 v1.01:

        28.0 °C : approximately 167.832°E, +2.020°
        28.5 °C : approximately 162.981°E, +1.537°
        29.0 °C : approximately 158.759°E, +1.084°

    The program recalculates these markers directly from the canonical
    Program-05 centroid CSV files, so it does not depend on a derived summary
    CSV being present.

THERMAL THRESHOLDS
    28.0 °C
    28.5 °C
    29.0 °C

LONG-TERM MEAN SST
    At each grid cell i:

                    sum over valid target dates SST_i(t)
        SSTbar_i = --------------------------------------
                    number of finite SST observations_i

    Missing values are excluded cell by cell.

    No threshold is applied before the temporal mean is calculated.

QUALITY CONTROL
    - exact target dates come from canonical Program-05 products;
    - the three centroid files must have identical target dates;
    - all expected 16,403 dates must be encountered exactly once;
    - duplicate SST dates from overlapping NetCDF products are ignored after
      first successful processing and reported;
    - coordinate grids must match the frozen Program-17 mask/reference grid;
    - anomaly files are excluded;
    - no interpolation or morphological operation is applied.

OUTPUTS
    Processed climatological arrays:
        data/processed/climatological_mean_sst_definition/
            pwp_long_term_mean_sst_definition_1981-09-01_2026-07-29.npz

    Tables:
        outputs/tables/threshold_comparison/
            pwp_climatological_mean_sst_definition/
                pwp_climatological_mean_sst_definition_summary.csv

    Figures:
        outputs/figures/threshold_comparison/
            pwp_climatological_mean_sst_definition/
                pwp_climatological_mean_sst_definition_A-E.png
                pwp_climatological_mean_sst_definition_A-E.pdf

    Reports:
        outputs/reports/threshold_comparison/
            pwp_climatological_mean_sst_definition/
                PROGRAM17_PWP_CLIMATOLOGICAL_MEAN_SST_DEFINITION.txt
                PROGRAM17_PWP_CLIMATOLOGICAL_MEAN_SST_DEFINITION.json

FIGURE INTERPRETATION
    A) Long-term mean OISST
       Background tropical-Pacific SST climatology.

    B) Fixed Pacific mask
       The geographical ocean domain within which PWP membership is allowed.

    C) Masked long-term mean SST
       Background SST climatology restricted to the analysed Pacific domain.

    D) Mean SST + 28/28.5/29 °C isotherms
       Direct visual comparison of the three thermal boundaries in the
       time-mean SST field.

    E) Mean-SST threshold regions + daily-population mean spherical centroids
       Visualizes the nested mean-SST thermal domains and where the mean
       spherical centroid of each DAILY threshold-defined population lies.

PUBLICATION WARNING
    Do not label Panel E as "frequency", "persistence", or "probability of
    belonging to the PWP". It is a thresholding of the time-mean SST field and
    serves only as long-term thermal-background context.

DEPENDENCIES
    Python 3.10+
    numpy
    pandas
    xarray
    matplotlib

EXECUTION
    From project root, for example:

        python src/threshold_analysis/17_pwp_climatological_mean_sst_definition_v1.01.py

===============================================================================
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import xarray as xr


# =============================================================================
# PROJECT DISCOVERY
# =============================================================================

SCRIPT_FILE = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_FILE.parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# LOAD CANONICAL PROGRAM 17 AS AN UPSTREAM UTILITY MODULE
# =============================================================================

PROGRAM17_CANDIDATES = (
    PROJECT_ROOT / "src_a" / "17_methodological_domain_figure.py",
    PROJECT_ROOT / "src" / "17_methodological_domain_figure.py",
)

PROGRAM17_SOURCE = next(
    (
        path
        for path in PROGRAM17_CANDIDATES
        if path.is_file()
    ),
    None,
)

if PROGRAM17_SOURCE is None:
    raise FileNotFoundError(
        "Canonical Program-17 source was not found.\n"
        "Checked:\n"
        + "\n".join(
            f"  - {path}"
            for path in PROGRAM17_CANDIDATES
        )
    )

program17_spec = importlib.util.spec_from_file_location(
    "osaf_program17_climatology_source",
    PROGRAM17_SOURCE,
)

if (
    program17_spec is None
    or program17_spec.loader is None
):
    raise RuntimeError(
        f"Unable to create import specification for {PROGRAM17_SOURCE}"
    )

P17 = importlib.util.module_from_spec(
    program17_spec
)

sys.modules[
    "osaf_program17_climatology_source"
] = P17

program17_spec.loader.exec_module(
    P17
)


# =============================================================================
# PROGRAM CONSTANTS
# =============================================================================

PROGRAM_NAME = (
    "OSAF PROGRAM 17 EXTENSION — LONG-TERM MEAN SST / "
    "PWP THERMAL-DEFINITION FIGURE"
)
PROGRAM_VERSION = "1.01"

START_DATE = pd.Timestamp("1981-09-01")
END_DATE = pd.Timestamp("2026-07-29")
EXPECTED_DAYS = 16_403

THRESHOLDS_C = (
    28.0,
    28.5,
    29.0,
)

PRIMARY_THRESHOLD_C = 28.0

MAP_LONGITUDE_MIN = float(
    getattr(
        P17,
        "MAP_LONGITUDE_MIN",
        100.0,
    )
)

MAP_LONGITUDE_MAX = float(
    getattr(
        P17,
        "MAP_LONGITUDE_MAX",
        290.0,
    )
)

MAP_LATITUDE_MIN = float(
    getattr(
        P17,
        "MAP_LATITUDE_MIN",
        -35.0,
    )
)

MAP_LATITUDE_MAX = float(
    getattr(
        P17,
        "MAP_LATITUDE_MAX",
        35.0,
    )
)

FIGURE_DPI = int(
    getattr(
        P17,
        "FIGURE_DPI",
        300,
    )
)

SAVE_BBOX = getattr(
    P17,
    "SAVE_BBOX",
    "tight",
)

SAVE_PAD_INCHES = float(
    getattr(
        P17,
        "SAVE_PAD_INCHES",
        0.05,
    )
)

SAVE_TRANSPARENT = bool(
    getattr(
        P17,
        "SAVE_TRANSPARENT",
        False,
    )
)

GRID_LINESTYLE = getattr(
    P17,
    "GRID_LINESTYLE",
    "--",
)

GRID_LINEWIDTH = float(
    getattr(
        P17,
        "GRID_LINEWIDTH",
        0.5,
    )
)

GRID_ALPHA = float(
    getattr(
        P17,
        "GRID_ALPHA",
        0.35,
    )
)

RAW_OISST_DIR = Path(
    P17.RAW_OISST_DIR
)

PACIFIC_MASK_FILE = Path(
    P17.PACIFIC_MASK_FILE
)

LATITUDE_GRID_FILE = Path(
    P17.LATITUDE_GRID_FILE
)

LONGITUDE_GRID_FILE = Path(
    P17.LONGITUDE_GRID_FILE
)

CENTROID_FILES = {
    28.0:
    PROJECT_ROOT
    / "data"
    / "processed"
    / "28"
    / "centroid"
    / "pwp_centroid_series.csv",

    28.5:
    PROJECT_ROOT
    / "data"
    / "processed"
    / "28.5"
    / "centroid"
    / "pwp_centroid_series.csv",

    29.0:
    PROJECT_ROOT
    / "data"
    / "processed"
    / "29"
    / "centroid"
    / "pwp_centroid_series.csv",
}

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "climatological_mean_sst_definition"
)

TABLE_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "tables"
    / "threshold_comparison"
    / "pwp_climatological_mean_sst_definition"
)

FIGURE_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "figures"
    / "threshold_comparison"
    / "pwp_climatological_mean_sst_definition"
)

REPORT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reports"
    / "threshold_comparison"
    / "pwp_climatological_mean_sst_definition"
)

NPZ_FILE = (
    PROCESSED_DIR
    / (
        "pwp_long_term_mean_sst_definition_"
        f"{START_DATE:%Y-%m-%d}_{END_DATE:%Y-%m-%d}.npz"
    )
)

SUMMARY_CSV = (
    TABLE_DIR
    / "pwp_climatological_mean_sst_definition_summary.csv"
)

FIGURE_PNG = (
    FIGURE_DIR
    / "pwp_climatological_mean_sst_definition_A-E.png"
)

FIGURE_PDF = (
    FIGURE_DIR
    / "pwp_climatological_mean_sst_definition_A-E.pdf"
)

REPORT_TXT = (
    REPORT_DIR
    / "PROGRAM17_PWP_CLIMATOLOGICAL_MEAN_SST_DEFINITION.txt"
)

REPORT_JSON = (
    REPORT_DIR
    / "PROGRAM17_PWP_CLIMATOLOGICAL_MEAN_SST_DEFINITION.json"
)

CHUNK_DAYS = 8


# =============================================================================
# TERMINAL UTILITIES
# =============================================================================

def rule(
    character: str = "=",
    width: int = 78,
) -> None:
    print(
        character
        * width
    )


def section(
    title: str,
) -> None:
    print()
    print(
        title
    )
    rule(
        "-"
    )


def item(
    label: str,
    value: object,
) -> None:
    print(
        f"{label:<46s}: {value}"
    )


# =============================================================================
# GENERIC UTILITIES
# =============================================================================

def sha256(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open(
        "rb"
    ) as stream:
        for block in iter(
            lambda: stream.read(
                1024
                * 1024
            ),
            b"",
        ):
            digest.update(
                block
            )

    return digest.hexdigest()


def identify_name(
    available: Iterable[str],
    candidates: tuple[str, ...],
    description: str,
) -> str:
    return P17.identify_name(
        available,
        candidates,
        description,
    )


def normalize_longitude_grid(
    longitude: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    return P17.normalize_longitude_grid(
        longitude
    )


def longitude_360_to_180(
    longitude_360: float,
) -> float:
    return (
        (
            float(
                longitude_360
            )
            + 180.0
        )
        % 360.0
    ) - 180.0


# =============================================================================
# CONFIGURATION / CANONICAL DATE VALIDATION
# =============================================================================

def validate_configuration() -> None:
    required = (
        PROGRAM17_SOURCE,
        RAW_OISST_DIR,
        PACIFIC_MASK_FILE,
        LATITUDE_GRID_FILE,
        LONGITUDE_GRID_FILE,
        *CENTROID_FILES.values(),
    )

    missing = [
        path
        for path in required
        if not path.exists()
    ]

    if missing:
        raise FileNotFoundError(
            "Required inputs are missing:\n"
            + "\n".join(
                f"  - {path}"
                for path in missing
            )
        )

    for directory in (
        PROCESSED_DIR,
        TABLE_DIR,
        FIGURE_DIR,
        REPORT_DIR,
    ):
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def load_common_target_dates() -> pd.DatetimeIndex:
    date_sets: dict[
        float,
        pd.DatetimeIndex,
    ] = {}

    for threshold_c, source in CENTROID_FILES.items():
        data = pd.read_csv(
            source,
            usecols=[
                "date",
                "threshold_c",
            ],
        )

        data[
            "date"
        ] = pd.to_datetime(
            data[
                "date"
            ],
            errors="raise",
        ).dt.normalize()

        if data[
            "date"
        ].duplicated().any():
            raise ValueError(
                f"Duplicate dates in canonical centroid file: {source}"
            )

        if len(
            data
        ) != EXPECTED_DAYS:
            raise ValueError(
                f"Expected {EXPECTED_DAYS:,} rows in {source}; "
                f"found {len(data):,}."
            )

        threshold_values = pd.to_numeric(
            data[
                "threshold_c"
            ],
            errors="coerce",
        ).dropna().to_numpy(
            dtype=float
        )

        if not np.allclose(
            threshold_values,
            threshold_c,
        ):
            raise ValueError(
                f"Threshold-column mismatch in {source}"
            )

        dates = pd.DatetimeIndex(
            data[
                "date"
            ]
        )

        if (
            dates.min()
            != START_DATE
            or dates.max()
            != END_DATE
        ):
            raise ValueError(
                "Canonical date range mismatch for "
                f"{threshold_c:.1f} °C:\n"
                f"  expected {START_DATE:%Y-%m-%d} to {END_DATE:%Y-%m-%d}\n"
                f"  found    {dates.min():%Y-%m-%d} to {dates.max():%Y-%m-%d}"
            )

        date_sets[
            threshold_c
        ] = dates

    reference = date_sets[
        PRIMARY_THRESHOLD_C
    ]

    for threshold_c in THRESHOLDS_C:
        current = date_sets[
            threshold_c
        ]

        if not reference.equals(
            current
        ):
            raise ValueError(
                "Canonical Program-05 date sequences are not identical "
                f"between 28.0 and {threshold_c:.1f} °C."
            )

    return reference


# =============================================================================
# REFERENCE GRID / MASK
# =============================================================================

def load_reference_grid_and_mask() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    mask = np.load(
        PACIFIC_MASK_FILE
    ).astype(
        bool
    )

    latitude = np.load(
        LATITUDE_GRID_FILE
    ).astype(
        np.float64
    ).squeeze()

    longitude = np.load(
        LONGITUDE_GRID_FILE
    ).astype(
        np.float64
    ).squeeze()

    longitude, order = normalize_longitude_grid(
        longitude
    )

    mask = mask[
        :,
        order
    ]

    if latitude[
        0
    ] > latitude[
        -1
    ]:
        latitude = latitude[
            ::-1
        ]

        mask = mask[
            ::-1,
            :
        ]

    if mask.shape != (
        latitude.size,
        longitude.size,
    ):
        raise ValueError(
            "Reference mask/grid mismatch.\n"
            f"Mask      : {mask.shape}\n"
            f"Latitude  : {latitude.size}\n"
            f"Longitude : {longitude.size}"
        )

    return (
        latitude,
        longitude,
        mask,
    )


# =============================================================================
# RAW OISST FILE DISCOVERY
# =============================================================================

def is_probable_anomaly_file(
    path: Path,
) -> bool:
    return P17.is_probable_anomaly_file(
        path
    )


def discover_absolute_sst_files() -> tuple[
    Path,
    ...,
]:
    """
    Discover raw absolute-SST NetCDF files.

    Prefer canonical annual sst.day.mean.YYYY.nc files when the full target
    year set can be represented that way. Otherwise fall back to all
    non-anomaly NetCDF files and rely on exact target-date filtering.
    """

    all_files = sorted(
        [
            path
            for path in RAW_OISST_DIR.rglob(
                "*"
            )
            if (
                path.is_file()
                and path.suffix.lower()
                in getattr(
                    P17,
                    "ABSOLUTE_SST_FILE_SUFFIXES",
                    (
                        ".nc",
                        ".nc4",
                        ".cdf",
                    ),
                )
                and not is_probable_anomaly_file(
                    path
                )
            )
        ],
        key=lambda path: str(
            path
        ),
    )

    if not all_files:
        raise FileNotFoundError(
            f"No absolute-SST NetCDF files found under {RAW_OISST_DIR}"
        )

    years = range(
        START_DATE.year,
        END_DATE.year
        + 1,
    )

    annual_by_year: dict[
        int,
        list[
            Path
        ],
    ] = {}

    for year in years:
        exact_name = (
            f"sst.day.mean.{year}.nc"
        )

        matches = [
            path
            for path in all_files
            if path.name.lower()
            == exact_name.lower()
        ]

        if matches:
            annual_by_year[
                year
            ] = matches

    if len(
        annual_by_year
    ) == len(
        tuple(
            years
        )
    ):
        selected: list[
            Path
        ] = []

        for year in range(
            START_DATE.year,
            END_DATE.year
            + 1,
        ):
            matches = annual_by_year[
                year
            ]

            if len(
                matches
            ) > 1:
                raise ValueError(
                    "More than one exact annual OISST file was found for "
                    f"{year}:\n"
                    + "\n".join(
                        f"  - {path}"
                        for path in matches
                    )
                )

            selected.append(
                matches[
                    0
                ]
            )

        return tuple(
            selected
        )

    return tuple(
        all_files
    )


# =============================================================================
# STREAMED LONG-TERM MEAN SST
# =============================================================================

def calculate_long_term_mean_sst(
    target_dates: pd.DatetimeIndex,
    reference_latitude: np.ndarray,
    reference_longitude: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[
        Path
    ],
    int,
]:
    """
    Stream raw OISST in small time chunks and calculate an exact target-date
    cell-wise mean.

    Every target date is accepted only once.
    """

    candidate_files = discover_absolute_sst_files()

    target_date_set = set(
        target_dates
    )

    processed_dates: set[
        pd.Timestamp
    ] = set()

    duplicate_dates_skipped = 0

    sum_sst = np.zeros(
        (
            reference_latitude.size,
            reference_longitude.size,
        ),
        dtype=np.float64,
    )

    valid_count = np.zeros(
        (
            reference_latitude.size,
            reference_longitude.size,
        ),
        dtype=np.uint32,
    )

    used_files: list[
        Path
    ] = []

    section(
        "## CALCULATING LONG-TERM MEAN DAILY SST"
    )

    for file_index, file_path in enumerate(
        candidate_files,
        start=1,
    ):
        with xr.open_dataset(
            file_path,
            decode_times=True,
            mask_and_scale=True,
        ) as dataset:
            try:
                sst_name = identify_name(
                    dataset.data_vars,
                    P17.SST_VARIABLE_CANDIDATES,
                    "absolute SST variable",
                )
            except Exception:
                # Some unrelated NetCDF files may exist under data/raw.
                continue

            sst_data = dataset[
                sst_name
            ]

            latitude_name = identify_name(
                (
                    list(
                        dataset.coords
                    )
                    + list(
                        sst_data.dims
                    )
                ),
                P17.LATITUDE_COORDINATE_CANDIDATES,
                "latitude coordinate",
            )

            longitude_name = identify_name(
                (
                    list(
                        dataset.coords
                    )
                    + list(
                        sst_data.dims
                    )
                ),
                P17.LONGITUDE_COORDINATE_CANDIDATES,
                "longitude coordinate",
            )

            time_name = None

            for candidate in P17.TIME_COORDINATE_CANDIDATES:
                for available in (
                    list(
                        dataset.coords
                    )
                    + list(
                        sst_data.dims
                    )
                ):
                    if (
                        available.lower()
                        == candidate.lower()
                    ):
                        time_name = available
                        break

                if time_name is not None:
                    break

            if time_name is None:
                # The climatological calculation requires explicit dates.
                continue

            raw_times = pd.DatetimeIndex(
                pd.to_datetime(
                    dataset[
                        time_name
                    ].values
                )
            ).normalize()

            relevant_indices: list[
                int
            ] = []

            relevant_dates: list[
                pd.Timestamp
            ] = []

            for index, date in enumerate(
                raw_times
            ):
                normalized = pd.Timestamp(
                    date
                ).normalize()

                if normalized not in target_date_set:
                    continue

                if normalized in processed_dates:
                    duplicate_dates_skipped += 1
                    continue

                relevant_indices.append(
                    index
                )

                relevant_dates.append(
                    normalized
                )

            if not relevant_indices:
                continue

            latitude = np.asarray(
                dataset[
                    latitude_name
                ].values,
                dtype=np.float64,
            ).squeeze()

            longitude_raw = np.asarray(
                dataset[
                    longitude_name
                ].values,
                dtype=np.float64,
            ).squeeze()

            longitude, longitude_order = normalize_longitude_grid(
                longitude_raw
            )

            latitude_reverse = bool(
                latitude[
                    0
                ]
                > latitude[
                    -1
                ]
            )

            if latitude_reverse:
                latitude = latitude[
                    ::-1
                ]

            if not np.allclose(
                latitude,
                reference_latitude,
                rtol=0.0,
                atol=1.0e-6,
            ):
                raise ValueError(
                    f"Latitude grid mismatch in {file_path}"
                )

            if not np.allclose(
                longitude,
                reference_longitude,
                rtol=0.0,
                atol=1.0e-6,
            ):
                raise ValueError(
                    f"Longitude grid mismatch in {file_path}"
                )

            used_files.append(
                file_path
            )

            for chunk_start in range(
                0,
                len(
                    relevant_indices
                ),
                CHUNK_DAYS,
            ):
                chunk_indices = relevant_indices[
                    chunk_start:
                    chunk_start
                    + CHUNK_DAYS
                ]

                chunk_dates = relevant_dates[
                    chunk_start:
                    chunk_start
                    + CHUNK_DAYS
                ]

                chunk = (
                    sst_data
                    .isel(
                        {
                            time_name:
                            chunk_indices
                        }
                    )
                    .transpose(
                        time_name,
                        latitude_name,
                        longitude_name,
                    )
                )

                values = np.asarray(
                    chunk.values,
                    dtype=np.float64,
                )

                values = values[
                    :,
                    :,
                    longitude_order
                ]

                if latitude_reverse:
                    values = values[
                        :,
                        ::-1,
                        :
                    ]

                finite = np.isfinite(
                    values
                )

                sum_sst += np.nansum(
                    values,
                    axis=0,
                )

                valid_count += np.sum(
                    finite,
                    axis=0,
                    dtype=np.uint32,
                )

                processed_dates.update(
                    chunk_dates
                )

        if (
            file_index
            % 5
            == 0
            or file_index
            == len(
                candidate_files
            )
        ):
            print(
                f"Files scanned: {file_index:4d}/{len(candidate_files):4d} | "
                f"unique target dates accumulated: {len(processed_dates):,}"
            )

    missing_dates = sorted(
        target_date_set
        - processed_dates
    )

    if missing_dates:
        preview = "\n".join(
            f"  - {date:%Y-%m-%d}"
            for date in missing_dates[
                :20
            ]
        )

        raise RuntimeError(
            f"{len(missing_dates):,} canonical target dates were not found "
            "in the raw absolute-SST archive.\n"
            f"First missing dates:\n{preview}"
        )

    if len(
        processed_dates
    ) != EXPECTED_DAYS:
        raise RuntimeError(
            "Unexpected processed-date count.\n"
            f"Expected : {EXPECTED_DAYS:,}\n"
            f"Found    : {len(processed_dates):,}"
        )

    mean_sst = np.full(
        sum_sst.shape,
        np.nan,
        dtype=np.float64,
    )

    valid_cells = (
        valid_count
        > 0
    )

    mean_sst[
        valid_cells
    ] = (
        sum_sst[
            valid_cells
        ]
        / valid_count[
            valid_cells
        ]
    )

    return (
        mean_sst,
        valid_count,
        used_files,
        duplicate_dates_skipped,
    )


# =============================================================================
# SPHERICAL MEAN OF DAILY CANONICAL CENTROID POPULATIONS
# =============================================================================

def spherical_mean_centroid(
    longitude_360_deg: np.ndarray,
    latitude_deg: np.ndarray,
) -> tuple[
    float,
    float,
    float,
]:
    lon = np.deg2rad(
        np.asarray(
            longitude_360_deg,
            dtype=np.float64,
        )
    )

    lat = np.deg2rad(
        np.asarray(
            latitude_deg,
            dtype=np.float64,
        )
    )

    valid = (
        np.isfinite(
            lon
        )
        & np.isfinite(
            lat
        )
    )

    lon = lon[
        valid
    ]

    lat = lat[
        valid
    ]

    if lon.size == 0:
        raise ValueError(
            "No finite daily centroid coordinates."
        )

    x = np.cos(
        lat
    ) * np.cos(
        lon
    )

    y = np.cos(
        lat
    ) * np.sin(
        lon
    )

    z = np.sin(
        lat
    )

    xbar = float(
        np.mean(
            x
        )
    )

    ybar = float(
        np.mean(
            y
        )
    )

    zbar = float(
        np.mean(
            z
        )
    )

    resultant = float(
        np.sqrt(
            xbar**2
            + ybar**2
            + zbar**2
        )
    )

    longitude = float(
        np.rad2deg(
            np.arctan2(
                ybar,
                xbar,
            )
        )
        % 360.0
    )

    latitude = float(
        np.rad2deg(
            np.arctan2(
                zbar,
                np.sqrt(
                    xbar**2
                    + ybar**2
                ),
            )
        )
    )

    return (
        longitude,
        latitude,
        resultant,
    )


def load_mean_daily_centroids() -> pd.DataFrame:
    rows = []

    for threshold_c, source in CENTROID_FILES.items():
        data = pd.read_csv(
            source,
            usecols=[
                "date",
                "lon_360",
                "lat",
                "area_km2",
            ],
        )

        data[
            "date"
        ] = pd.to_datetime(
            data[
                "date"
            ],
            errors="raise",
        ).dt.normalize()

        if len(
            data
        ) != EXPECTED_DAYS:
            raise ValueError(
                f"Unexpected centroid record count in {source}"
            )

        longitude, latitude, resultant = spherical_mean_centroid(
            pd.to_numeric(
                data[
                    "lon_360"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            pd.to_numeric(
                data[
                    "lat"
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
        )

        area = pd.to_numeric(
            data[
                "area_km2"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        rows.append(
            {
                "threshold_c":
                threshold_c,

                "n_days":
                len(
                    data
                ),

                "mean_spherical_lon_360_deg":
                longitude,

                "mean_spherical_lat_deg":
                latitude,

                "mean_spherical_resultant_length":
                resultant,

                "mean_daily_pwp_area_km2":
                float(
                    np.nanmean(
                        area
                    )
                ),

                "median_daily_pwp_area_km2":
                float(
                    np.nanmedian(
                        area
                    )
                ),

                "centroid_source_file":
                str(
                    source
                ),
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            "threshold_c"
        )
        .reset_index(
            drop=True
        )
    )


# =============================================================================
# CLIMATOLOGICAL THRESHOLD SUMMARY
# =============================================================================

def spherical_cell_areas_km2(
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    return P17.spherical_cell_areas_km2(
        latitude,
        longitude,
    )


def build_summary(
    mean_sst: np.ndarray,
    pacific_mask: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    mean_centroids: pd.DataFrame,
) -> pd.DataFrame:
    cell_area = spherical_cell_areas_km2(
        latitude,
        longitude,
    )

    rows = []

    for threshold_c in THRESHOLDS_C:
        mean_field_domain = (
            pacific_mask
            & np.isfinite(
                mean_sst
            )
            & (
                mean_sst
                >= threshold_c
            )
        )

        mean_field_area = float(
            np.sum(
                cell_area[
                    mean_field_domain
                ]
            )
        )

        centroid_row = mean_centroids.loc[
            np.isclose(
                mean_centroids[
                    "threshold_c"
                ].to_numpy(
                    dtype=float
                ),
                threshold_c,
            )
        ].iloc[
            0
        ]

        rows.append(
            {
                "threshold_c":
                threshold_c,

                "target_start_date":
                START_DATE,

                "target_end_date":
                END_DATE,

                "n_target_days":
                EXPECTED_DAYS,

                "mean_sst_threshold_domain_cell_count":
                int(
                    np.count_nonzero(
                        mean_field_domain
                    )
                ),

                "mean_sst_threshold_domain_area_km2":
                mean_field_area,

                "mean_sst_threshold_domain_area_million_km2":
                mean_field_area
                / 1.0e6,

                "mean_daily_spherical_centroid_lon_360_deg":
                float(
                    centroid_row[
                        "mean_spherical_lon_360_deg"
                    ]
                ),

                "mean_daily_spherical_centroid_lat_deg":
                float(
                    centroid_row[
                        "mean_spherical_lat_deg"
                    ]
                ),

                "mean_daily_pwp_area_km2":
                float(
                    centroid_row[
                        "mean_daily_pwp_area_km2"
                    ]
                ),

                "median_daily_pwp_area_km2":
                float(
                    centroid_row[
                        "median_daily_pwp_area_km2"
                    ]
                ),

                "mean_field_domain_interpretation":
                (
                    "threshold applied to long-term mean SST; "
                    "NOT daily occurrence frequency"
                ),

                "centroid_interpretation":
                (
                    "spherical mean of canonical daily Program-05 "
                    "centroid positions; NOT centroid of mean-SST domain"
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# FIGURE UTILITIES
# =============================================================================

def add_panel_label(
    axis: plt.Axes,
    label: str,
) -> None:
    axis.text(
        0.015,
        0.975,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        bbox={
            "facecolor":
            "white",

            "edgecolor":
            "none",

            "alpha":
            0.76,

            "pad":
            2.0,
        },
        zorder=30,
    )


def configure_map_axis(
    axis: plt.Axes,
    title: str,
    show_x_label: bool = True,
    show_y_label: bool = True,
) -> None:
    axis.set_xlim(
        MAP_LONGITUDE_MIN,
        MAP_LONGITUDE_MAX,
    )

    axis.set_ylim(
        MAP_LATITUDE_MIN,
        MAP_LATITUDE_MAX,
    )

    axis.set_xticks(
        np.arange(
            100.0,
            291.0,
            30.0,
        )
    )

    axis.set_yticks(
        np.arange(
            -30.0,
            31.0,
            15.0,
        )
    )

    if show_x_label:
        axis.set_xlabel(
            "Longitude (°E)"
        )
    else:
        axis.set_xlabel(
            ""
        )

    if show_y_label:
        axis.set_ylabel(
            "Latitude (°)"
        )
    else:
        axis.set_ylabel(
            ""
        )

    axis.set_title(
        title,
        fontsize=10,
    )

    axis.grid(
        linestyle=GRID_LINESTYLE,
        linewidth=GRID_LINEWIDTH,
        alpha=GRID_ALPHA,
    )

    axis.set_facecolor(
        "0.92"
    )


def plot_climatological_figure(
    mean_sst: np.ndarray,
    pacific_mask: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    mean_centroids: pd.DataFrame,
) -> tuple[
    Path,
    Path,
]:
    lon2d, lat2d = np.meshgrid(
        longitude,
        latitude,
    )

    mean_pacific = np.where(
        pacific_mask,
        mean_sst,
        np.nan,
    )

    figure = plt.figure(
        figsize=(
            16.0,
            13.2,
        ),
        constrained_layout=True,
    )

    grid = figure.add_gridspec(
        nrows=3,
        ncols=2,
        height_ratios=(
            1.0,
            1.0,
            1.05,
        ),
    )

    axes = [
        figure.add_subplot(
            grid[
                0,
                0
            ]
        ),
        figure.add_subplot(
            grid[
                0,
                1
            ]
        ),
        figure.add_subplot(
            grid[
                1,
                0
            ]
        ),
        figure.add_subplot(
            grid[
                1,
                1
            ]
        ),
        figure.add_subplot(
            grid[
                2,
                :
            ]
        ),
    ]

    # -------------------------------------------------------------------------
    # A) Long-term mean OISST
    # -------------------------------------------------------------------------
    pcm_a = axes[
        0
    ].pcolormesh(
        longitude,
        latitude,
        mean_sst,
        shading="auto",
        vmin=20.0,
        vmax=31.0,
    )

    configure_map_axis(
        axes[
            0
        ],
        (
            "Long-term mean daily OISST\n"
            f"{START_DATE:%Y-%m-%d} to {END_DATE:%Y-%m-%d}"
        ),
        show_x_label=False,
        show_y_label=True,
    )

    add_panel_label(
        axes[
            0
        ],
        "A)",
    )

    cbar_a = figure.colorbar(
        pcm_a,
        ax=axes[
            0
        ],
        fraction=0.035,
        pad=0.02,
    )

    cbar_a.set_label(
        "Mean SST (°C)"
    )

    # -------------------------------------------------------------------------
    # B) Pacific mask
    # -------------------------------------------------------------------------
    mask_display = np.where(
        pacific_mask,
        1.0,
        np.nan,
    )

    axes[
        1
    ].pcolormesh(
        longitude,
        latitude,
        mask_display,
        shading="auto",
        cmap=ListedColormap(
            [
                "0.25"
            ]
        ),
        vmin=0.0,
        vmax=1.0,
    )

    configure_map_axis(
        axes[
            1
        ],
        "Fixed Pacific Ocean mask",
        show_x_label=False,
        show_y_label=False,
    )

    add_panel_label(
        axes[
            1
        ],
        "B)",
    )

    # -------------------------------------------------------------------------
    # C) Mean SST restricted to Pacific mask
    # -------------------------------------------------------------------------
    pcm_c = axes[
        2
    ].pcolormesh(
        longitude,
        latitude,
        mean_pacific,
        shading="auto",
        vmin=20.0,
        vmax=31.0,
    )

    configure_map_axis(
        axes[
            2
        ],
        "Long-term mean SST restricted to the Pacific mask",
        show_x_label=False,
        show_y_label=True,
    )

    add_panel_label(
        axes[
            2
        ],
        "C)",
    )

    cbar_c = figure.colorbar(
        pcm_c,
        ax=axes[
            2
        ],
        fraction=0.035,
        pad=0.02,
    )

    cbar_c.set_label(
        "Mean SST (°C)"
    )

    # -------------------------------------------------------------------------
    # D) Mean Pacific SST + three isotherms
    # -------------------------------------------------------------------------
    pcm_d = axes[
        3
    ].pcolormesh(
        longitude,
        latitude,
        mean_pacific,
        shading="auto",
        vmin=20.0,
        vmax=31.0,
    )

    contour = axes[
        3
    ].contour(
        lon2d,
        lat2d,
        mean_pacific,
        levels=THRESHOLDS_C,
        linewidths=1.5,
    )

    axes[
        3
    ].clabel(
        contour,
        inline=True,
        fmt=lambda value: f"{value:.1f} °C",
        fontsize=8,
    )

    configure_map_axis(
        axes[
            3
        ],
        "Thermal boundaries in the long-term mean SST field",
        show_x_label=False,
        show_y_label=False,
    )

    add_panel_label(
        axes[
            3
        ],
        "D)",
    )

    cbar_d = figure.colorbar(
        pcm_d,
        ax=axes[
            3
        ],
        fraction=0.035,
        pad=0.02,
    )

    cbar_d.set_label(
        "Mean SST (°C)"
    )

    # -------------------------------------------------------------------------
    # E) Nested threshold regions in mean SST + mean DAILY spherical centroids
    # -------------------------------------------------------------------------
    class_field = np.full(
        mean_sst.shape,
        np.nan,
        dtype=float,
    )

    valid_pacific = (
        pacific_mask
        & np.isfinite(
            mean_sst
        )
    )

    # Mutually exclusive classes for visual clarity.
    class_field[
        valid_pacific
        & (
            mean_sst
            >= 28.0
        )
        & (
            mean_sst
            < 28.5
        )
    ] = 1.0

    class_field[
        valid_pacific
        & (
            mean_sst
            >= 28.5
        )
        & (
            mean_sst
            < 29.0
        )
    ] = 2.0

    class_field[
        valid_pacific
        & (
            mean_sst
            >= 29.0
        )
    ] = 3.0

    class_cmap = ListedColormap(
        [
            "#d8e9f2",
            "#8fc5df",
            "#3182bd",
        ]
    )

    class_norm = BoundaryNorm(
        [
            0.5,
            1.5,
            2.5,
            3.5,
        ],
        class_cmap.N,
    )

    axes[
        4
    ].pcolormesh(
        longitude,
        latitude,
        class_field,
        shading="auto",
        cmap=class_cmap,
        norm=class_norm,
    )

    # Draw the exact mean-field boundaries as contours too.
    axes[
        4
    ].contour(
        lon2d,
        lat2d,
        mean_pacific,
        levels=THRESHOLDS_C,
        linewidths=0.85,
    )

    marker_specs = {
        28.0:
        "o",

        28.5:
        "s",

        29.0:
        "*",
    }

    for threshold_c in THRESHOLDS_C:
        row = mean_centroids.loc[
            np.isclose(
                mean_centroids[
                    "threshold_c"
                ].to_numpy(
                    dtype=float
                ),
                threshold_c,
            )
        ].iloc[
            0
        ]

        axes[
            4
        ].scatter(
            [
                float(
                    row[
                        "mean_spherical_lon_360_deg"
                    ]
                )
            ],
            [
                float(
                    row[
                        "mean_spherical_lat_deg"
                    ]
                )
            ],
            marker=marker_specs[
                threshold_c
            ],
            s=(
                190
                if threshold_c
                < 29.0
                else 270
            ),
            facecolors="white",
            edgecolors="black",
            linewidths=1.25,
            zorder=20,
            label=(
                f"{threshold_c:.1f} °C mean daily spherical centroid"
            ),
        )

    configure_map_axis(
        axes[
            4
        ],
        (
            "Mean-SST thermal domains and long-term mean spherical centroids "
            "of the daily PWP populations"
        ),
        show_x_label=True,
        show_y_label=True,
    )

    add_panel_label(
        axes[
            4
        ],
        "E)",
    )

    domain_legend = [
        Line2D(
            [
                0
            ],
            [
                0
            ],
            marker="s",
            linestyle="none",
            markerfacecolor="#d8e9f2",
            markeredgecolor="none",
            markersize=10,
            label="28.0 ≤ mean SST < 28.5 °C",
        ),
        Line2D(
            [
                0
            ],
            [
                0
            ],
            marker="s",
            linestyle="none",
            markerfacecolor="#8fc5df",
            markeredgecolor="none",
            markersize=10,
            label="28.5 ≤ mean SST < 29.0 °C",
        ),
        Line2D(
            [
                0
            ],
            [
                0
            ],
            marker="s",
            linestyle="none",
            markerfacecolor="#3182bd",
            markeredgecolor="none",
            markersize=10,
            label="mean SST ≥ 29.0 °C",
        ),
    ]

    handles, labels = axes[
        4
    ].get_legend_handles_labels()

    axes[
        4
    ].legend(
        domain_legend
        + handles,
        [
            item.get_label()
            for item in domain_legend
        ]
        + labels,
        frameon=False,
        loc="lower left",
        ncol=2,
        fontsize=8.2,
    )

    figure.suptitle(
        (
            "Long-term mean SST background and thermal definitions "
            "of the Pacific Warm Pool"
        ),
        fontsize=14,
    )

    figure.savefig(
        FIGURE_PNG,
        dpi=FIGURE_DPI,
        bbox_inches=SAVE_BBOX,
        pad_inches=SAVE_PAD_INCHES,
        transparent=SAVE_TRANSPARENT,
    )

    figure.savefig(
        FIGURE_PDF,
        bbox_inches=SAVE_BBOX,
        pad_inches=SAVE_PAD_INCHES,
        transparent=SAVE_TRANSPARENT,
    )

    plt.close(
        figure
    )

    return (
        FIGURE_PNG,
        FIGURE_PDF,
    )


# =============================================================================
# EXPORTS / REPORTING
# =============================================================================

def export_processed_arrays(
    mean_sst: np.ndarray,
    valid_count: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    pacific_mask: np.ndarray,
) -> Path:
    np.savez_compressed(
        NPZ_FILE,
        mean_sst_c=mean_sst,
        valid_count=valid_count,
        latitude=latitude,
        longitude_360=longitude,
        pacific_mask=pacific_mask.astype(
            np.uint8
        ),
        start_date=str(
            START_DATE.date()
        ),
        end_date=str(
            END_DATE.date()
        ),
        n_target_days=np.int64(
            EXPECTED_DAYS
        ),
        thresholds_c=np.asarray(
            THRESHOLDS_C,
            dtype=np.float64,
        ),
        interpretation=(
            "long-term mean SST background; thresholding this array is NOT "
            "climatological daily PWP occurrence frequency"
        ),
    )

    return NPZ_FILE


def export_summary(
    summary: pd.DataFrame,
) -> Path:
    output = summary.copy()

    for column in (
        "target_start_date",
        "target_end_date",
    ):
        output[
            column
        ] = pd.to_datetime(
            output[
                column
            ]
        ).dt.strftime(
            "%Y-%m-%d"
        )

    output.to_csv(
        SUMMARY_CSV,
        index=False,
        float_format="%.10g",
    )

    return SUMMARY_CSV


def write_reports(
    target_dates: pd.DatetimeIndex,
    used_files: list[
        Path
    ],
    duplicate_dates_skipped: int,
    mean_sst: np.ndarray,
    valid_count: np.ndarray,
    pacific_mask: np.ndarray,
    summary: pd.DataFrame,
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

    pacific_counts = valid_count[
        pacific_mask
    ]

    payload = {
        "program":
        PROGRAM_NAME,

        "version":
        PROGRAM_VERSION,

        "generated_utc":
        generated,

        "project_root":
        str(
            PROJECT_ROOT
        ),

        "program17_source":
        str(
            PROGRAM17_SOURCE
        ),

        "program17_sha256":
        sha256(
            PROGRAM17_SOURCE
        ),

        "raw_oisst_directory":
        str(
            RAW_OISST_DIR
        ),

        "pacific_mask":
        str(
            PACIFIC_MASK_FILE
        ),

        "latitude_grid":
        str(
            LATITUDE_GRID_FILE
        ),

        "longitude_grid":
        str(
            LONGITUDE_GRID_FILE
        ),

        "target_start_date":
        str(
            target_dates.min().date()
        ),

        "target_end_date":
        str(
            target_dates.max().date()
        ),

        "n_target_dates":
        int(
            len(
                target_dates
            )
        ),

        "thresholds_c":
        list(
            THRESHOLDS_C
        ),

        "n_raw_files_used":
        int(
            len(
                used_files
            )
        ),

        "raw_files_used":
        [
            str(
                path
            )
            for path in used_files
        ],

        "duplicate_target_dates_skipped":
        int(
            duplicate_dates_skipped
        ),

        "minimum_valid_observations_inside_pacific_mask":
        int(
            np.min(
                pacific_counts
            )
        ),

        "median_valid_observations_inside_pacific_mask":
        float(
            np.median(
                pacific_counts
            )
        ),

        "maximum_valid_observations_inside_pacific_mask":
        int(
            np.max(
                pacific_counts
            )
        ),

        "mean_sst_definition":
        (
            "cell-wise temporal mean of finite daily absolute OISST over the "
            "canonical target-date record"
        ),

        "critical_warning":
        (
            "thresholding mean SST is not climatological PWP occurrence. "
            "Occurrence must be calculated from the daily threshold condition."
        ),

        "summary":
        summary.to_dict(
            orient="records"
        ),

        "created_files":
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
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )

    lines = [
        PROGRAM_NAME,
        "=" * 78,
        "",
        f"Version                         : {PROGRAM_VERSION}",
        f"Generated UTC                   : {generated}",
        f"Project root                    : {PROJECT_ROOT}",
        f"Canonical Program-17 source     : {PROGRAM17_SOURCE}",
        f"Raw OISST directory             : {RAW_OISST_DIR}",
        f"Pacific mask                    : {PACIFIC_MASK_FILE}",
        "",
        "CANONICAL TEMPORAL RECORD",
        "-" * 78,
        (
            f"Period                          : "
            f"{target_dates.min():%Y-%m-%d} to {target_dates.max():%Y-%m-%d}"
        ),
        f"Target dates                    : {len(target_dates):,}",
        f"Raw SST files used              : {len(used_files):,}",
        f"Duplicate target dates skipped  : {duplicate_dates_skipped:,}",
        "",
        "SCIENTIFIC DISTINCTION",
        "-" * 78,
        (
            "This figure describes the long-term MEAN SST background. "
            "It does not estimate climatological daily PWP occurrence."
        ),
        "",
        (
            "Mean-SST threshold domain:"
        ),
        (
            "    1[ mean_t(SST_i) >= T ]"
        ),
        "",
        (
            "Daily occurrence/persistence (separate analysis):"
        ),
        (
            "    sum_t 1[SST_i(t) >= T] / sum_t 1[SST_i(t) finite]"
        ),
        "",
        "THRESHOLD SUMMARY",
        "-" * 78,
    ]

    for row in summary.itertuples(
        index=False
    ):
        lines.append(
            (
                f"{row.threshold_c:.1f} °C | "
                f"mean-SST domain area="
                f"{row.mean_sst_threshold_domain_area_million_km2:.3f} "
                f"×10⁶ km² | "
                f"mean daily spherical centroid="
                f"({row.mean_daily_spherical_centroid_lon_360_deg:.3f}°E, "
                f"{row.mean_daily_spherical_centroid_lat_deg:+.3f}°)"
            )
        )

    lines.extend(
        [
            "",
            "FILES CREATED",
            "-" * 78,
        ]
    )

    lines.extend(
        str(
            path
        )
        for path in created_files
    )

    lines.extend(
        [
            "",
            "=" * 78,
            "PROGRAM COMPLETED SUCCESSFULLY.",
            "",
        ]
    )

    REPORT_TXT.write_text(
        "\n".join(
            lines
        ),
        encoding="utf-8",
    )

    return (
        REPORT_TXT,
        REPORT_JSON,
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

    validate_configuration()

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
        "Canonical Program-17 source",
        PROGRAM17_SOURCE,
    )

    item(
        "Raw OISST directory",
        RAW_OISST_DIR,
    )

    item(
        "Pacific mask",
        PACIFIC_MASK_FILE,
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
        "Interpretation",
        "long-term mean SST background; NOT occurrence frequency",
    )

    section(
        "## LOADING CANONICAL TARGET DATES"
    )

    target_dates = load_common_target_dates()

    item(
        "Record",
        (
            f"{target_dates.min():%Y-%m-%d} to "
            f"{target_dates.max():%Y-%m-%d}"
        ),
    )

    item(
        "Canonical target days",
        f"{len(target_dates):,}",
    )

    section(
        "## LOADING FIXED PACIFIC MASK / REFERENCE GRID"
    )

    (
        latitude,
        longitude,
        pacific_mask,
    ) = load_reference_grid_and_mask()

    item(
        "Grid shape",
        pacific_mask.shape,
    )

    item(
        "Pacific mask cells",
        f"{np.count_nonzero(pacific_mask):,}",
    )

    (
        mean_sst,
        valid_count,
        used_files,
        duplicate_dates_skipped,
    ) = calculate_long_term_mean_sst(
        target_dates=target_dates,
        reference_latitude=latitude,
        reference_longitude=longitude,
    )

    section(
        "## LONG-TERM MEAN SPHERICAL CENTROIDS FROM FROZEN DAILY PRODUCTS"
    )

    mean_centroids = load_mean_daily_centroids()

    for row in mean_centroids.itertuples(
        index=False
    ):
        print(
            f"{row.threshold_c:.1f} °C | "
            f"N={row.n_days:,} | "
            f"mean spherical centroid="
            f"({row.mean_spherical_lon_360_deg:.3f}°E, "
            f"{row.mean_spherical_lat_deg:+.3f}°) | "
            f"mean daily area="
            f"{row.mean_daily_pwp_area_km2 / 1.0e6:.3f} ×10⁶ km²"
        )

    section(
        "## BUILDING CLIMATOLOGICAL THRESHOLD SUMMARY"
    )

    summary = build_summary(
        mean_sst=mean_sst,
        pacific_mask=pacific_mask,
        latitude=latitude,
        longitude=longitude,
        mean_centroids=mean_centroids,
    )

    for row in summary.itertuples(
        index=False
    ):
        print(
            f"{row.threshold_c:.1f} °C | "
            f"area of [mean SST >= T] domain="
            f"{row.mean_sst_threshold_domain_area_million_km2:.3f} ×10⁶ km²"
        )

    section(
        "## EXPORTING PROCESSED ARRAYS AND TABLE"
    )

    created_files: list[
        Path
    ] = []

    created_files.append(
        export_processed_arrays(
            mean_sst=mean_sst,
            valid_count=valid_count,
            latitude=latitude,
            longitude=longitude,
            pacific_mask=pacific_mask,
        )
    )

    created_files.append(
        export_summary(
            summary
        )
    )

    for path in created_files:
        print(
            path
        )

    section(
        "## GENERATING PUBLICATION A-E FIGURE"
    )

    figure_files = plot_climatological_figure(
        mean_sst=mean_sst,
        pacific_mask=pacific_mask,
        latitude=latitude,
        longitude=longitude,
        mean_centroids=mean_centroids,
    )

    created_files.extend(
        figure_files
    )

    for path in figure_files:
        print(
            path
        )

    section(
        "## WRITING REPORTS"
    )

    report_files = write_reports(
        target_dates=target_dates,
        used_files=used_files,
        duplicate_dates_skipped=duplicate_dates_skipped,
        mean_sst=mean_sst,
        valid_count=valid_count,
        pacific_mask=pacific_mask,
        summary=summary,
        created_files=created_files,
    )

    created_files.extend(
        report_files
    )

    for path in report_files:
        print(
            path
        )

    print()
    rule()
    print(
        "PROGRAM 17 CLIMATOLOGICAL EXTENSION COMPLETED SUCCESSFULLY."
    )
    print()
    print(
        "A-E figure generated:"
    )
    print(
        "A) long-term mean OISST"
    )
    print(
        "B) fixed Pacific mask"
    )
    print(
        "C) masked long-term mean Pacific SST"
    )
    print(
        "D) 28.0/28.5/29.0 °C isotherms in the mean SST field"
    )
    print(
        "E) mean-SST thermal domains + mean DAILY spherical centroids"
    )
    print()
    print(
        "WARNING: this program does NOT calculate climatological PWP "
        "occurrence/persistence."
    )
    print(
        "The occurrence diagnostic must threshold each DAILY SST field before "
        "temporal averaging."
    )
    rule()


if __name__ == "__main__":
    main()

