"""Shared helpers: text normalisation, ancestry codes, subprocess, JSON IO."""

from __future__ import annotations

import contextlib
import json
import math
import os
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}

# Ancestries with a 1000 Genomes LD panel (super-population codes in the panel file).
ANCESTRY_LABELS = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Hispanic or Latin American",
}

ANCESTRY_ALIASES = {
    "eur": "EUR", "european": "EUR",
    "afr": "AFR", "african": "AFR",
    "eas": "EAS", "east asian": "EAS",
    "sas": "SAS", "south asian": "SAS",
    "amr": "AMR", "american": "AMR", "hispanic": "AMR", "latino": "AMR",
    "hispanic or latin american": "AMR",
}


class Gwas2RnaError(RuntimeError):
    """An expected, user-facing failure (printed without a traceback)."""


def banner(text: str, width: int = 88) -> None:
    print()
    print("=" * width)
    print(text)
    print("=" * width, flush=True)


def normalize_text(value: Any) -> str:
    value = "" if value is None else str(value)
    value = re.sub(r"[^a-z0-9]+", " ", value.strip().lower())
    return re.sub(r"\s+", " ", value).strip()


def slugify(value: Any) -> str:
    return normalize_text(value).replace(" ", "_")


def canonical_ancestry(value: str) -> tuple[str, str]:
    """Return (code, label), e.g. ('EUR', 'European')."""
    code = ANCESTRY_ALIASES.get(normalize_text(value))
    if code is None:
        raise Gwas2RnaError(
            f"Unsupported ancestry {value!r}. Use one of: {', '.join(ANCESTRY_LABELS)} "
            "(these are the ancestries with a 1000 Genomes LD reference)."
        )
    return code, ANCESTRY_LABELS[code]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def write_json(path: Path, payload: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Per-process temporary name: parallel jobs may write the same file.
    temporary = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(temporary, path)


def nonempty(path: Path) -> bool:
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def run_command(command, log_file: Path | None = None) -> str:
    """Run a command, echo its output, optionally save it, raise on failure."""
    command = [str(x) for x in command]
    print("$ " + " ".join(shlex.quote(x) for x in command), flush=True)
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    captured = []
    handle = None
    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handle = open(log_file, "w", encoding="utf-8")
    try:
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            captured.append(line)
            if handle:
                handle.write(line)
    finally:
        if handle:
            handle.close()
    code = process.wait()
    output = "".join(captured)
    if code != 0:
        raise RuntimeError(
            f"Command failed with exit code {code}:\n{' '.join(command)}\n\n{output[-5000:]}"
        )
    return output


def find_column(columns, candidates):
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def clean_chr(value: Any) -> int | None:
    text = str(value).strip().lower().replace("chr", "")
    try:
        chrom = int(float(text))
    except (TypeError, ValueError):
        return None
    return chrom if 1 <= chrom <= 22 else None


def clean_allele(value: Any) -> str | None:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    value = str(value).strip().upper()
    if not value or value in {"NA", "NAN", "<NA>", "NONE", "."}:
        return None
    return value


def safe_int(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, (int,)):
        return value if value > 0 else None
    if isinstance(value, float):
        return int(round(value)) if math.isfinite(value) and value > 0 else None
    match = re.search(r"\d+", str(value).strip().replace(",", ""))
    if not match:
        return None
    number = int(match.group())
    return number if number > 0 else None


def parse_numbers(text: Any, minimum: int = 100) -> list[int]:
    """All integers >= minimum in free text such as '10,000 cases, 20,000 controls'."""
    if text is None or (not isinstance(text, str) and pd.isna(text)):
        return []
    values = []
    for token in re.findall(r"\d[\d,]*", str(text)):
        try:
            number = int(token.replace(",", ""))
        except ValueError:
            continue
        if number >= minimum:
            values.append(number)
    return values


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def to_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


@contextlib.contextmanager
def file_lock(path: Path):
    """Exclusive lock so concurrent jobs do not download/build the same file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import fcntl
    except ImportError:  # Windows: no cross-process locking, run setup once first.
        yield
        return
    with open(path, "w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)
