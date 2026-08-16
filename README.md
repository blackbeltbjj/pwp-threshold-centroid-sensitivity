# PWP Threshold–Centroid Sensitivity — Paper 1 Reproducibility Package

Computational reproducibility package for PWP Paper 1.

## Frozen analysis
- Period: 1981-09-01 to 2026-07-29
- Common daily observations: 16,403
- OISST resolution: 0.25° × 0.25°
- SST thresholds: 28.0, 28.5, and 29.0 °C
- Fixed historical Pacific mask
- Spherical, area-weighted centroid methodology

## Repository status
Pre-publication reproducibility release candidate.

## Core dependency
Program 05 is preserved as an immutable publication snapshot under:

src/core_snapshot/05_calculate_pwp_centroid.py

SHA-256:
729B80855247DE4F690C790DB3129262BB2E808FC43B9088D9E83B00E3DA3150

## Data
Raw NOAA OISST files are not distributed in this repository.
Exact acquisition and provenance documentation will be recorded in docs/DATA_PROVENANCE.md.

## Reproduction
Detailed instructions: docs/REPRODUCIBILITY.md

## Citation
Citation metadata and DOI will be added only after final authorship approval and Zenodo archival.

## Manuscript relationship

This repository is the frozen computational reproducibility package for PWP
Paper 1 targeted to the Journal of Atmospheric and Oceanic Technology (JTECH).

The exact manuscript title will be inserted only after confirmation against
the final JTECH Science-Freeze manuscript.

Program 33 is intentionally excluded from this release because Paper 1 does
not depend on the ENSO-association analysis.

## License

MIT License.
