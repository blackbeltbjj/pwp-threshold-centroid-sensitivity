# -*- coding: utf-8 -*-
"""
===============================================================================
PROJECT
    Pacific Warm Pool Scientific Analysis Platform

MODULE
    core/reporting.py

VERSION
    1.1.0

PURPOSE
    Provide the authoritative, dependency-free console and plain-text reporting
    primitives used throughout the platform.

DESIGN GOALS
    - preserve backward compatibility with rule(), section(), and item();
    - standardize banners, progress, warnings, failures, and success messages;
    - provide an object-oriented ConsoleReport interface;
    - provide a deterministic TextReport builder for technical reports;
    - remain independent of scientific algorithms and third-party packages;
    - support dependency injection of output streams for testing.

PUBLIC API
    DEFAULT_WIDTH
    DEFAULT_KEY_WIDTH
    MessageLevel
    rule
    banner
    section
    item
    progress
    warning
    success
    failure
    render_rule
    render_banner
    render_section
    render_key_value
    render_status
    ConsoleReport
    TextReport

BACKWARD COMPATIBILITY
    Existing calls in Programs 11, 18, 19, and 20 remain valid:

        rule()
        section("CONFIGURATION")
        item("Program version", "4.0.2")

PYTHON
    Python 3.10+

AUTHOR
    Fabio Vieira Machado
===============================================================================
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, TextIO


DEFAULT_WIDTH = 78
DEFAULT_KEY_WIDTH = 31


class MessageLevel(str, Enum):
    """Standardized semantic level for console status messages."""

    PROGRESS = "PROGRESS"
    WARNING = "WARNING"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


def _validate_width(width: int) -> int:
    """Validate and return a positive report width."""

    if not isinstance(width, int):
        raise TypeError("Report width must be an integer.")

    if width < 1:
        raise ValueError("Report width must be positive.")

    return width


def _validate_key_width(key_width: int) -> int:
    """Validate and return a non-negative key-column width."""

    if not isinstance(key_width, int):
        raise TypeError("Key width must be an integer.")

    if key_width < 0:
        raise ValueError("Key width cannot be negative.")

    return key_width


def _rule_character(character: str) -> str:
    """Return the first character used to construct a horizontal rule."""

    if not isinstance(character, str):
        raise TypeError("Rule character must be a string.")

    if not character:
        raise ValueError("Rule character cannot be empty.")

    return character[0]


def render_rule(
    character: str = "=",
    width: int = DEFAULT_WIDTH,
) -> str:
    """Return a deterministic horizontal rule."""

    return _rule_character(character) * _validate_width(width)


def render_banner(
    title: object,
    width: int = DEFAULT_WIDTH,
    character: str = "=",
) -> str:
    """Return a three-line banner."""

    horizontal_rule = render_rule(
        character=character,
        width=width,
    )

    return "\n".join(
        (
            horizontal_rule,
            str(title),
            horizontal_rule,
        )
    )


def render_section(
    title: object,
    width: int = DEFAULT_WIDTH,
) -> str:
    """Return a section heading without a leading blank line."""

    return "\n".join(
        (
            str(title),
            render_rule("-", width),
        )
    )


def render_key_value(
    key: object,
    value: object,
    key_width: int = DEFAULT_KEY_WIDTH,
) -> str:
    """Return one aligned key-value line."""

    validated_width = _validate_key_width(key_width)

    return f"{str(key):<{validated_width}s}: {value}"


def render_status(
    message: object,
    level: MessageLevel | str,
) -> str:
    """Return a standardized status line."""

    try:
        normalized_level = MessageLevel(level)
    except ValueError as error:
        allowed = ", ".join(
            member.value
            for member in MessageLevel
        )
        raise ValueError(
            f"Unsupported message level {level!r}. "
            f"Allowed values: {allowed}."
        ) from error

    return f"[{normalized_level.value}] {message}"


def _write_line(
    text: object = "",
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Write one line to the selected stream."""

    print(
        str(text),
        file=stream,
    )


def rule(
    character: str = "=",
    width: int = DEFAULT_WIDTH,
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Print a horizontal rule."""

    _write_line(
        render_rule(
            character=character,
            width=width,
        ),
        stream=stream,
    )


def banner(
    title: object,
    width: int = DEFAULT_WIDTH,
    character: str = "=",
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Print a standardized program banner."""

    _write_line(
        render_banner(
            title=title,
            width=width,
            character=character,
        ),
        stream=stream,
    )


def section(
    title: object,
    width: int = DEFAULT_WIDTH,
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """
    Print a section heading.

    A leading blank line is intentionally preserved for backward compatibility
    with the established Program 11, 18, 19, and 20 console layout.
    """

    _write_line(
        "",
        stream=stream,
    )

    _write_line(
        render_section(
            title=title,
            width=width,
        ),
        stream=stream,
    )


def item(
    key: object,
    value: object,
    key_width: int = DEFAULT_KEY_WIDTH,
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Print one aligned key-value pair."""

    _write_line(
        render_key_value(
            key=key,
            value=value,
            key_width=key_width,
        ),
        stream=stream,
    )


def progress(
    message: object,
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Print a standardized progress message."""

    _write_line(
        render_status(
            message,
            MessageLevel.PROGRESS,
        ),
        stream=stream,
    )


def warning(
    message: object,
    *,
    stream: TextIO = sys.stderr,
) -> None:
    """Print a standardized warning message."""

    _write_line(
        render_status(
            message,
            MessageLevel.WARNING,
        ),
        stream=stream,
    )


def success(
    message: object,
    *,
    stream: TextIO = sys.stdout,
) -> None:
    """Print a standardized success message."""

    _write_line(
        render_status(
            message,
            MessageLevel.SUCCESS,
        ),
        stream=stream,
    )


def failure(
    message: object,
    *,
    stream: TextIO = sys.stderr,
) -> None:
    """Print a standardized failure message."""

    _write_line(
        render_status(
            message,
            MessageLevel.FAILURE,
        ),
        stream=stream,
    )


@dataclass
class ConsoleReport:
    """
    Stateful console-reporting facade.

    The class does not store scientific results. It only centralizes layout,
    widths, output streams, and semantic status messages.
    """

    width: int = DEFAULT_WIDTH
    key_width: int = DEFAULT_KEY_WIDTH
    output_stream: TextIO = sys.stdout
    error_stream: TextIO = sys.stderr

    def __post_init__(self) -> None:
        self.width = _validate_width(self.width)
        self.key_width = _validate_key_width(self.key_width)

    def rule(
        self,
        character: str = "=",
    ) -> None:
        rule(
            character=character,
            width=self.width,
            stream=self.output_stream,
        )

    def banner(
        self,
        title: object,
        character: str = "=",
    ) -> None:
        banner(
            title=title,
            width=self.width,
            character=character,
            stream=self.output_stream,
        )

    def section(
        self,
        title: object,
    ) -> None:
        section(
            title=title,
            width=self.width,
            stream=self.output_stream,
        )

    def item(
        self,
        key: object,
        value: object,
    ) -> None:
        item(
            key=key,
            value=value,
            key_width=self.key_width,
            stream=self.output_stream,
        )

    def progress(
        self,
        message: object,
    ) -> None:
        progress(
            message,
            stream=self.output_stream,
        )

    def warning(
        self,
        message: object,
    ) -> None:
        warning(
            message,
            stream=self.error_stream,
        )

    def success(
        self,
        message: object,
    ) -> None:
        success(
            message,
            stream=self.output_stream,
        )

    def failure(
        self,
        message: object,
    ) -> None:
        failure(
            message,
            stream=self.error_stream,
        )


@dataclass
class TextReport:
    """
    Deterministic plain-text report builder.

    The builder deliberately avoids automatic timestamps or environment details.
    Those belong to the metadata/provenance layer and must be supplied
    explicitly by the calling program.
    """

    title: str
    width: int = DEFAULT_WIDTH
    lines: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.width = _validate_width(self.width)

        if not isinstance(self.title, str):
            raise TypeError("Report title must be a string.")

        if not self.title.strip():
            raise ValueError("Report title cannot be empty.")

        self.lines.extend(
            (
                self.title,
                render_rule("=", self.width),
            )
        )

    def add(
        self,
        line: object = "",
    ) -> "TextReport":
        """Append one line and return self for fluent composition."""

        self.lines.append(
            str(line)
        )

        return self

    def add_blank(self) -> "TextReport":
        """Append one blank line."""

        self.lines.append("")

        return self

    def add_rule(
        self,
        character: str = "-",
    ) -> "TextReport":
        """Append one horizontal rule."""

        self.lines.append(
            render_rule(
                character=character,
                width=self.width,
            )
        )

        return self

    def add_section(
        self,
        title: object,
    ) -> "TextReport":
        """Append a blank line and a standardized section heading."""

        self.lines.extend(
            (
                "",
                render_section(
                    title=title,
                    width=self.width,
                ),
            )
        )

        return self

    def add_item(
        self,
        key: object,
        value: object,
        key_width: int = DEFAULT_KEY_WIDTH,
    ) -> "TextReport":
        """Append one aligned key-value line."""

        self.lines.append(
            render_key_value(
                key=key,
                value=value,
                key_width=key_width,
            )
        )

        return self

    def add_items(
        self,
        values: Iterable[tuple[object, object]],
        key_width: int = DEFAULT_KEY_WIDTH,
    ) -> "TextReport":
        """Append multiple aligned key-value lines."""

        for key, value in values:
            self.add_item(
                key=key,
                value=value,
                key_width=key_width,
            )

        return self

    def add_status(
        self,
        message: object,
        level: MessageLevel | str,
    ) -> "TextReport":
        """Append one standardized status line."""

        self.lines.append(
            render_status(
                message=message,
                level=level,
            )
        )

        return self

    def extend(
        self,
        lines: Iterable[object],
    ) -> "TextReport":
        """Append an iterable of lines."""

        self.lines.extend(
            str(line)
            for line in lines
        )

        return self

    def text(
        self,
        final_newline: bool = True,
    ) -> str:
        """Render the complete report."""

        content = "\n".join(
            self.lines
        )

        if final_newline:
            return content + "\n"

        return content

    def write(
        self,
        path: Path | str,
        encoding: str = "utf-8",
    ) -> Path:
        """Write the report atomically where possible."""

        output = Path(path)
        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        temporary = output.with_name(
            f".{output.name}.tmp"
        )

        temporary.write_text(
            self.text(),
            encoding=encoding,
        )

        temporary.replace(
            output
        )

        return output
