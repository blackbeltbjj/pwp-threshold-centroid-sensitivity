from pathlib import Path
import hashlib

EXPECTED = "729B80855247DE4F690C790DB3129262BB2E808FC43B9088D9E83B00E3DA3150"

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "src" / "core_snapshot" / "05_calculate_pwp_centroid.py"

def test_program05_snapshot_hash():
    digest = hashlib.sha256(TARGET.read_bytes()).hexdigest().upper()
    assert digest == EXPECTED
