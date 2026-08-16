from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "outputs" / "reference_checks"

DOMAIN = pd.read_csv(REF / "pwp_domain_size_centroid_sensitivity_summary.csv")
SEP = pd.read_csv(REF / "pwp_interthreshold_centroid_separation_summary.csv")
OCC = pd.read_csv(REF / "pwp_occurrence_persistence_summary.csv")
CONN = pd.read_csv(REF / "pwp_connectivity_threshold_summary.csv")


def row_for(df, threshold):
    return df.loc[df["threshold_c"].astype(float) == float(threshold)].iloc[0]


@pytest.mark.parametrize(
    "threshold,expected",
    [
        (28.0, 37.5006534),
        (28.5, 28.03081125),
        (29.0, 18.54671425),
    ],
)
def test_full_domain_mean_area(threshold, expected):
    row = row_for(DOMAIN, threshold)
    assert row["mean_area_million_km2"] == pytest.approx(expected, abs=1e-8)


@pytest.mark.parametrize(
    "threshold,expected_lon,expected_lat",
    [
        (28.0, 167.8320978, 2.019889073),
        (28.5, 162.9812324, 1.536551344),
        (29.0, 158.7592261, 1.083571183),
    ],
)
def test_full_domain_spherical_centroids(threshold, expected_lon, expected_lat):
    row = row_for(DOMAIN, threshold)
    assert row["reference_lon_360"] == pytest.approx(expected_lon, abs=1e-7)
    assert row["reference_lat"] == pytest.approx(expected_lat, abs=1e-7)


def test_full_domain_28_29_centroid_separation():
    row = SEP[
        (SEP["threshold_1_c"].astype(float) == 28.0)
        & (SEP["threshold_2_c"].astype(float) == 29.0)
    ].iloc[0]

    assert row["n_days"] == 16403
    assert row["median_separation_km"] == pytest.approx(1021.064658, abs=1e-6)
    assert row["separation_p95_km"] == pytest.approx(1869.337271, abs=1e-6)


@pytest.mark.parametrize(
    "threshold,expected",
    [
        (28.0, 0.15857021),
        (28.5, 0.11833991),
        (29.0, 0.07817878),
    ],
)
def test_occurrence_mean_gridcell_fraction(threshold, expected):
    row = row_for(OCC, threshold)
    assert row["n_common_days"] == 16403
    assert row["mean_gridcell_occurrence_fraction"] == pytest.approx(expected, abs=1e-8)


@pytest.mark.parametrize("threshold", [28.0, 28.5, 29.0])
def test_frozen_period_and_record_count(threshold):
    d = row_for(DOMAIN, threshold)
    c = row_for(CONN, threshold)
    o = row_for(OCC, threshold)

    assert d["n_days"] == 16403
    assert c["n_days"] == 16403
    assert o["n_common_days"] == 16403

    assert d["start_date"] == "1981-09-01"
    assert d["end_date"] == "2026-07-29"
    assert c["start_date"] == "1981-09-01"
    assert c["end_date"] == "2026-07-29"
    assert o["start_date"] == "1981-09-01"
    assert o["end_date"] == "2026-07-29"
