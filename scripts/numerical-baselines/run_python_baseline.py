#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import configparser
import io
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[2]
PY_SRC = ROOT / "Python" / "rdpower" / "src"
RDROBUST_PY_SRC = ROOT.parent / "rdrobust" / "Python" / "rdrobust" / "src"
if RDROBUST_PY_SRC.exists():
    sys.path.insert(0, str(RDROBUST_PY_SRC))
sys.path.insert(0, str(PY_SRC))

import numpy as np
import pandas as pd

from rdpower import rdmde, rdpower, rdsampsi


def rdrobust_version() -> str | None:
    setup_cfg = ROOT.parent / "rdrobust" / "Python" / "rdrobust" / "setup.cfg"
    if setup_cfg.exists():
        config = configparser.ConfigParser()
        config.read(setup_cfg, encoding="utf-8")
        return config["metadata"].get("version")
    try:
        from importlib import metadata

        return metadata.version("rdrobust")
    except Exception:
        return None


def finite_or_none(value: Any) -> float | int | None:
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, np.ndarray) and value.size == 1:
        value = value.item()
    if isinstance(value, (int, float)):
        return value if math.isfinite(value) else None
    return value


def quiet_call(fn: Callable[..., dict[str, Any]], *args: Any, **kwargs: Any) -> dict[str, Any]:
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        return fn(*args, **kwargs)


def keep(result: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    return {field: finite_or_none(result.get(field)) for field in fields}


def summarize_rdpower(result: dict[str, Any]) -> dict[str, Any]:
    fields = [
        "power_rbc",
        "se_rbc",
        "power_conv",
        "se_conv",
        "sampsi_l",
        "sampsi_r",
        "samph_l",
        "samph_r",
        "N_l",
        "N_r",
        "Nh_l",
        "Nh_r",
        "tau",
        "bias_l",
        "bias_r",
        "Vl_rb",
        "Vr_rb",
        "alpha",
    ]
    return keep(result, fields)


def summarize_rdsampsi(result: dict[str, Any]) -> dict[str, Any]:
    fields = [
        "sampsi.h.tot",
        "sampsi.h.l",
        "sampsi.h.r",
        "sampsi.tot",
        "sampsi.h.tot.cl",
        "sampsi.h.l.cl",
        "sampsi.h.r.cl",
        "sampsi.tot.cl",
        "N.l",
        "N.r",
        "Nh.l",
        "Nh.r",
        "bias.l",
        "bias.r",
        "var.l",
        "var.r",
        "samph.l",
        "samph.r",
        "tau",
        "beta",
        "alpha",
        "init.cond",
        "no.iter",
    ]
    return keep(result, fields)


def summarize_rdmde(result: dict[str, Any]) -> dict[str, Any]:
    fields = [
        "mde",
        "mde_conv",
        "se_rbc",
        "se_conv",
        "sampsi_l",
        "sampsi_r",
        "samph_l",
        "samph_r",
        "N_l",
        "N_r",
        "Nh_l",
        "Nh_r",
        "bias_l",
        "bias_r",
        "Vl_rb",
        "Vr_rb",
        "alpha",
        "beta",
    ]
    return keep(result, fields)


def run_cases() -> dict[str, Any]:
    data = pd.read_csv(ROOT / "Python" / "rdpower_senate.csv")
    z = data[["demvoteshfor2", "demmv"]]
    covs = data[["population", "dopen", "dmidterm"]]
    cluster = data["state"]

    cases: dict[str, Any] = {}
    cases["senate_default"] = {
        "rdpower": summarize_rdpower(quiet_call(rdpower, data=z, tau=5)),
        "rdsampsi": summarize_rdsampsi(quiet_call(rdsampsi, data=z, tau=5)),
        "rdmde": summarize_rdmde(quiet_call(rdmde, data=z)),
    }
    cases["senate_covs"] = {
        "rdpower": summarize_rdpower(quiet_call(rdpower, data=z, tau=5, covs=covs)),
        "rdsampsi": summarize_rdsampsi(quiet_call(rdsampsi, data=z, tau=5, covs=covs)),
        "rdmde": summarize_rdmde(quiet_call(rdmde, data=z, covs=covs)),
    }
    cases["senate_fixed_bandwidth"] = {
        "rdpower": summarize_rdpower(quiet_call(rdpower, data=z, tau=5, h=[16, 18], b=[18, 20], all=True)),
        "rdsampsi": summarize_rdsampsi(quiet_call(rdsampsi, data=z, tau=5, beta=0.9, samph=[18, 19], nratio=0.5, all=True)),
        "rdmde": summarize_rdmde(quiet_call(rdmde, data=z, beta=0.75, samph=[12, 13])),
    }
    cases["senate_cluster"] = {
        "rdpower": summarize_rdpower(quiet_call(rdpower, data=z, tau=5, cluster=cluster)),
        "rdsampsi": summarize_rdsampsi(quiet_call(rdsampsi, data=z, tau=5, cluster=cluster)),
        "rdmde": summarize_rdmde(quiet_call(rdmde, data=z, cluster=cluster)),
    }
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Python rdpower numerical baselines.")
    parser.add_argument(
        "--output",
        default=str(ROOT / "docs" / "audit" / "baselines" / "python-current.json"),
        help="Output JSON path.",
    )
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "schema_version": 1,
        "package": "rdpower",
        "language": "python",
        "source": "working-tree",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "rdrobust": rdrobust_version(),
            "rdrobust_source": str(RDROBUST_PY_SRC) if RDROBUST_PY_SRC.exists() else "installed",
        },
        "cases": run_cases(),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
