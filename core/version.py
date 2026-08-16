# -*- coding: utf-8 -*-
"""Platform identity and version helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


PLATFORM_NAME: Final[str] = (
    "Pacific Warm Pool Scientific Analysis Platform"
)

PLATFORM_VERSION: Final[str] = "4.0.3"
CORE_VERSION: Final[str] = "1.1.0"


@dataclass(frozen=True)
class VersionInfo:
    """Immutable platform/core version descriptor."""

    platform_name: str = PLATFORM_NAME
    platform_version: str = PLATFORM_VERSION
    core_version: str = CORE_VERSION


def version_string() -> str:
    """Return a human-readable platform version."""

    return (
        f"{PLATFORM_NAME} {PLATFORM_VERSION} "
        f"(core {CORE_VERSION})"
    )
