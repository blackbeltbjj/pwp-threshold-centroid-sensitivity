# PWP Threshold-Centroid Sensitivity - Paper 1 Reproducibility Package

Computational reproducibility package supporting:

**Machado, F. V.**

*Defining the Pacific Warm Pool: Threshold Dependence of Centroid Geometry and Robustness of Interannual Variance Modulation.*

Journal of Atmospheric and Oceanic Technology manuscript.

## Scientific Release

**Current archival software release:** v1.0.0

[View v1.0.0 release](https://github.com/blackbeltbjj/pwp-threshold-centroid-sensitivity/releases/tag/v1.0.0)

This release preserves the frozen computational workflow supporting the threshold-sensitivity analysis used in the manuscript.

## Frozen Analysis

- Analysis period: 1981-09-01 to 2026-07-29
- Common daily observations: 16,403
- NOAA OISST v2.1
- Spatial resolution: 0.25 x 0.25 degrees
- SST thresholds: 28.0, 28.5, and 29.0 degrees C
- Fixed historical Pacific mask
- Spherical, physical-area-weighted centroid methodology
- Largest-connected-component diagnostics
- Time-frequency and robustness analyses

## Scientific Scope

The repository supports analyses of:

- threshold sensitivity
- spherical centroid geometry
- geodesic displacement
- connectivity
- largest connected components
- Fourier and wavelet variability
- occurrence and persistence
- inter-threshold spatial differences
- robustness of interannual variance
- reproducibility and numerical integrity

## Core Dependency

The authoritative Program 05 snapshot is preserved under:

`src/core_snapshot/05_calculate_pwp_centroid.py`

SHA-256:

`729B80855247DE4F690C790DB3129262BB2E808FC43B9088D9E83B00E3DA3150`

The reusable core methodology is maintained in:

[OSAF-PWP](https://github.com/blackbeltbjj/OSAF-PWP)

## Reproducibility Workflow

The scientific workflow combines the canonical PWP core with threshold-specific analyses:

`NOAA OISST v2.1 -> Pacific mask -> threshold fields -> spherical geometry -> centroid calculations -> connectivity -> sensitivity analyses -> time-frequency diagnostics -> occurrence and persistence -> audited outputs`

## Data

Raw NOAA OISST files are not distributed in this repository.

Data provenance, configuration information, and reproduction instructions are documented within the repository.

## Reproduction

See:

`docs/REPRODUCIBILITY.md`

## Scope Exclusion

Program 33 is intentionally excluded from this release because the manuscript does not depend on the ENSO-association analysis.

## Archival DOI

A Zenodo DOI will be added here only after the archive identifier has been verified against this released software version.

## Author

**Fabio Vieira Machado**

ORCID: [0000-0003-0723-075X](https://orcid.org/0000-0003-0723-075X)

## License

MIT License.
