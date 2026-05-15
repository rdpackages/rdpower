# rdpower Numerical Baselines

Local baseline harness for checking numerical stability across the R, Python,
and Stata implementations.

The scripts write JSON files under `docs/audit/baselines/` by default.
Reference `*-modernized.json` files are committed; local `*-current.json`
outputs are ignored.

## Cases

- `senate_default`: package Senate data with default `rdpower`, `rdsampsi`, and
  `rdmde`.
- `senate_covs`: package Senate data with covariates.
- `senate_fixed_bandwidth`: package Senate data with fixed bandwidth or sample
  design options.
- `senate_cluster`: package Senate data with cluster-robust calculations.

## Commands

From the repository root:

```powershell
python scripts/numerical-baselines/run_python_baseline.py
Rscript scripts/numerical-baselines/run_r_baseline.R
```

If Stata is available:

```powershell
& "C:\Program Files\Stata18\StataMP-64.exe" /e do scripts/numerical-baselines/run_stata_baseline.do
```

Then compare two JSON files:

```powershell
python scripts/numerical-baselines/compare_baselines.py docs/audit/baselines/r-current.json docs/audit/baselines/python-current.json --common-only
```

Use tight tolerances for baseline-to-baseline checks within the same language.
Use slightly looser tolerances for cross-language checks because package data
storage and dependency versions can differ across languages.
