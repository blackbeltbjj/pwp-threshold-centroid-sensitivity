from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "outputs"
    / "reference_checks"
    / "PROGRAM35_PWP_LCC_ONLY_FULL_RECORD_ROBUSTNESS.txt"
)

TEXT = REPORT.read_bytes().decode("cp1252")


@pytest.mark.parametrize(
    "threshold,lon,lat",
    [
        ("28.0", "163.675", "1.206"),
        ("28.5", "159.184", "0.636"),
        ("29.0", "156.925", "0.196"),
    ],
)
def test_lcc_spherical_centroids(threshold, lon, lat):
    assert f"{threshold} " in TEXT
    assert f"mean spherical LCC centroid=({lon}" in TEXT
    assert f"+{lat}" in TEXT


def test_lcc_28_29_same_day_separation():
    assert "28.0 vs 29.0" in TEXT
    assert "N=16,403" in TEXT
    assert "median=1015.2 km" in TEXT
    assert "P95=2937.8 km" in TEXT


@pytest.mark.parametrize(
    "metric,start,end,pct",
    [
        ("mean_lcc_area_million_km2", "34.777047", "15.005293", "-56.85%"),
        ("lcc_lon_p95_minus_p05_deg", "30.738599", "47.243073", "+53.69%"),
        ("lcc_lon_sd_deg", "9.434960", "15.312548", "+62.30%"),
        ("lcc_radial_p95_km", "2165.200710", "3503.603550", "+61.81%"),
    ],
)
def test_lcc_28_to_29_changes(metric, start, end, pct):
    assert f"{metric}: {start} -> {end} | {pct}" in TEXT


def test_program35_record_count():
    assert TEXT.count("N=16,403") >= 6
