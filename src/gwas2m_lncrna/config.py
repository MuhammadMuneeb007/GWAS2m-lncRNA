"""Configuration: one YAML file, defaults taken from the bundled template."""

from __future__ import annotations

import copy
import os
from importlib import resources as importlib_resources
from pathlib import Path
from typing import Any

import yaml

from gwas2m_lncrna.common import Gwas2RnaError

CONFIG_FILENAME = "gwas2m-lncrna.yaml"
CONFIG_ENV = "GWAS2M_LNCRNA_CONFIG"


def template_text() -> str:
    return (importlib_resources.files("gwas2m_lncrna") / "data" / "config_template.yaml").read_text(
        encoding="utf-8"
    )


def defaults() -> dict:
    return yaml.safe_load(template_text())


def _merge(base: dict, override: dict, where: str = "") -> dict:
    merged = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if key not in base:
            raise Gwas2RnaError(f"Unknown configuration key '{where}{key}'. Check spelling against `gwas2m-lncrna init`.")
        if isinstance(base[key], dict):
            if value is None:
                continue
            if not isinstance(value, dict):
                raise Gwas2RnaError(f"Configuration key '{where}{key}' must be a section (mapping).")
            merged[key] = _merge(base[key], value, f"{where}{key}.")
        else:
            merged[key] = value
    return merged


class Config:
    """Effective settings. Access sections like cfg["finemapping"]["coverage"]."""

    def __init__(self, data: dict, source: Path | None):
        self.data = data
        self.source = source
        self.base_dir = source.parent if source else Path.cwd()

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def path(self, value: str | os.PathLike | None) -> Path | None:
        if value in (None, ""):
            return None
        p = Path(os.path.expandvars(str(value))).expanduser()
        return (p if p.is_absolute() else self.base_dir / p).resolve()

    @property
    def output_dir(self) -> Path:
        return self.path(self.data["output_dir"])

    @property
    def resources_dir(self) -> Path:
        return self.path(self.data["resources_dir"])

    @property
    def threads(self) -> int:
        return max(1, int(self.data["threads"]))

    @property
    def jobs(self) -> int:
        return max(1, int(self.data["jobs"]))

    def describe(self) -> str:
        return str(self.source) if self.source else "built-in defaults (no gwas2m-lncrna.yaml found)"


def find_config_file(explicit: str | None) -> Path | None:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if not path.is_file():
            raise Gwas2RnaError(f"Configuration file not found: {path}")
        return path
    if os.environ.get(CONFIG_ENV):
        path = Path(os.environ[CONFIG_ENV]).expanduser().resolve()
        if not path.is_file():
            raise Gwas2RnaError(f"{CONFIG_ENV} points to a missing file: {path}")
        return path
    local = Path.cwd() / CONFIG_FILENAME
    return local.resolve() if local.is_file() else None


def load_config(explicit: str | None = None, overrides: dict | None = None) -> Config:
    source = find_config_file(explicit)
    user = {}
    if source is not None:
        try:
            user = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError as exc:
            raise Gwas2RnaError(f"Cannot parse {source}: {exc}") from exc
        if not isinstance(user, dict):
            raise Gwas2RnaError(f"{source} must contain a YAML mapping.")
    data = _merge(defaults(), user)
    for key, value in (overrides or {}).items():
        if value is None:
            continue
        section, _, name = key.partition(".")
        if name:
            data[section][name] = value
        else:
            data[section] = value
    return Config(data, source)


def write_template(path: Path, force: bool = False) -> Path:
    path = Path(path).resolve()
    if path.exists() and not force:
        raise Gwas2RnaError(f"{path} already exists (use --force to overwrite).")
    path.write_text(template_text(), encoding="utf-8")
    return path
