"""Repository-level checks for the rdpower multi-language package repo."""

from __future__ import annotations

import configparser
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def require_path(relative_path: str) -> Path:
    path = ROOT / relative_path
    if not path.exists():
        fail(f"Missing required path: {relative_path}")
    return path


def read_text(relative_path: str) -> str:
    return require_path(relative_path).read_text(encoding="utf-8")


def check_required_layout() -> None:
    required_paths = [
        "README.md",
        "LICENSE.md",
        "R/rdpower/DESCRIPTION",
        "Python/rdpower/pyproject.toml",
        "Python/rdpower/setup.cfg",
        "stata/stata.toc",
        "stata/rdpower.pkg",
    ]

    for relative_path in required_paths:
        require_path(relative_path)


def check_python_metadata() -> str:
    config = configparser.ConfigParser()
    config.read(ROOT / "Python/rdpower/setup.cfg", encoding="utf-8")

    if "metadata" not in config:
        fail("Python/rdpower/setup.cfg is missing [metadata].")

    name = config["metadata"].get("name", "").strip()
    version = config["metadata"].get("version", "").strip()

    if name != "rdpower":
        fail(f"Unexpected Python package name: {name!r}")
    if not version:
        fail("Python package version is empty.")

    return version


def check_r_metadata() -> str:
    description = read_text("R/rdpower/DESCRIPTION")
    package_match = re.search(r"^Package:\s*(\S+)\s*$", description, re.MULTILINE)
    version_match = re.search(r"^Version:\s*(\S+)\s*$", description, re.MULTILINE)

    if not package_match:
        fail("R/rdpower/DESCRIPTION is missing Package.")
    if package_match.group(1) != "rdpower":
        fail(f"Unexpected R package name: {package_match.group(1)!r}")
    if not version_match:
        fail("R/rdpower/DESCRIPTION is missing Version.")

    return version_match.group(1)


def check_stata_manifest() -> str:
    package_text = read_text("stata/rdpower.pkg")
    distribution_match = re.search(r"^d Distribution-Date:\s*(\d{8})\s*$", package_text, re.MULTILINE)

    if not distribution_match:
        fail("stata/rdpower.pkg is missing Distribution-Date.")

    missing_files = []
    for line in package_text.splitlines():
        if not line.startswith("f "):
            continue
        relative_file = line[2:].strip()
        if not (ROOT / "stata" / relative_file).is_file():
            missing_files.append(relative_file)

    if missing_files:
        missing = ", ".join(missing_files)
        fail(f"Files listed in stata/rdpower.pkg are missing: {missing}")

    return distribution_match.group(1)


def main() -> None:
    check_required_layout()
    r_version = check_r_metadata()
    python_version = check_python_metadata()
    stata_distribution_date = check_stata_manifest()

    print("Repository checks passed.")
    print(f"R package version: {r_version}")
    print(f"Python package version: {python_version}")
    print(f"Stata distribution date: {stata_distribution_date}")


if __name__ == "__main__":
    main()
