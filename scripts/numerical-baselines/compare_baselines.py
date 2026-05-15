#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def flatten(prefix: str, value: Any, out: dict[str, Any]) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            flatten(f"{prefix}.{key}" if prefix else key, child, out)
    else:
        out[prefix] = value


def as_number(value: Any) -> float | None:
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        return float(value)
    if isinstance(value, str):
        try:
            parsed = float(value.strip())
        except ValueError:
            return None
        return parsed if math.isfinite(parsed) else None
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two rdpower numerical baseline JSON files.")
    parser.add_argument("expected")
    parser.add_argument("actual")
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument(
        "--common-only",
        action="store_true",
        help="Compare only fields present in both files.",
    )
    args = parser.parse_args()

    expected = json.loads(Path(args.expected).read_text(encoding="utf-8"))
    actual = json.loads(Path(args.actual).read_text(encoding="utf-8"))

    left: dict[str, Any] = {}
    right: dict[str, Any] = {}
    flatten("", expected.get("cases", {}), left)
    flatten("", actual.get("cases", {}), right)

    failures: list[str] = []
    all_keys = sorted(set(left) | set(right))
    for key in all_keys:
        if key not in left:
            if args.common_only:
                continue
            failures.append(f"missing in expected: {key}")
            continue
        if key not in right:
            if args.common_only:
                continue
            failures.append(f"missing in actual: {key}")
            continue
        a = left[key]
        b = right[key]
        a_num = as_number(a)
        b_num = as_number(b)
        if a_num is not None and b_num is not None:
            if not math.isclose(a_num, b_num, rel_tol=args.rtol, abs_tol=args.atol):
                failures.append(f"{key}: expected {a_num:.17g}, actual {b_num:.17g}, diff {b_num - a_num:.3g}")
        elif a != b:
            failures.append(f"{key}: expected {a!r}, actual {b!r}")

    if failures:
        for failure in failures:
            print(failure)
        print(f"Compared {len(all_keys)} fields; {len(failures)} differences outside tolerance.")
        return 1

    print(f"Compared {len(all_keys)} fields; all matched within tolerance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
