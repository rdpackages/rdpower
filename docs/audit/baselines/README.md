# Numerical Baselines

Reference JSON outputs for the release-readiness numerical audit.

- `*-modernized.json` files are checked into the repository as the current reference baseline.
- `*-current.json` files are generated locally by the scripts in `scripts/numerical-baselines/` and are ignored by Git.

Regenerate local outputs from the repository root with:

```powershell
python scripts/numerical-baselines/run_python_baseline.py
Rscript scripts/numerical-baselines/run_r_baseline.R
& "C:\Program Files\Stata18\StataMP-64.exe" /e do scripts/numerical-baselines/run_stata_baseline.do
```

Then compare with:

```powershell
python scripts/numerical-baselines/compare_baselines.py docs/audit/baselines/python-modernized.json docs/audit/baselines/python-current.json
```
