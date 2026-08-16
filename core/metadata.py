# -*- coding: utf-8 -*-
"""Checksums, stable identifiers, JSON serialization, and provenance metadata."""
from __future__ import annotations
import hashlib, json, platform, re, sys
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
import numpy as np
import pandas as pd

def sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()

def slug(text: object) -> str:
    value = re.sub(r"[^a-zA-Z0-9]+", "_", str(text)).strip("_").lower()
    if not value:
        raise ValueError("Cannot create a slug from an empty identifier.")
    return value

def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()

def json_default(value: Any) -> Any:
    if isinstance(value, Path): return str(value)
    if isinstance(value, (datetime, pd.Timestamp)): return value.isoformat()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, (np.integer, np.floating)): return value.item()
    if isinstance(value, np.bool_): return bool(value)
    if is_dataclass(value): return asdict(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable.")

def write_json(path: Path, payload: Any, *, indent: int = 2) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=indent, ensure_ascii=False, sort_keys=True, default=json_default)+"\n", encoding="utf-8")
    return path

@dataclass(frozen=True)
class RuntimeMetadata:
    generated_utc: str
    python_version: str
    python_executable: str
    platform: str

    @classmethod
    def capture(cls) -> "RuntimeMetadata":
        return cls(utc_now_iso(), sys.version.split()[0], sys.executable, platform.platform())

def file_manifest(paths: Iterable[Path]) -> list[dict[str, Any]]:
    result=[]
    for candidate in paths:
        path=Path(candidate)
        result.append({"path":str(path), "exists":path.exists(), "size_bytes":path.stat().st_size if path.is_file() else None, "sha256":sha256(path) if path.is_file() else None})
    return result
