# Changelog

## 1.3.0

### Added
- `--help` and `--version` on `sk-run`, `sk-prepare` and `sk-synth`; `sk-run --quiet`; `sk-synth --seed N --outdir DIR`.
- Exit codes 0/1/2 with a single `<prog>: error: ...` line on stderr instead of tracebacks.
- Validation of `.dat` files and config keys (`exchange`, `model`, `init.mode`).
- `influential` key in the jackknife JSON.
- Relative dataset paths fall back to the config file's directory.
- Golden-value tests (`pytest`) and this changelog.
- A warning on stderr when the optimizer doesn't converge.

### Changed
- Derived uncertainties propagate the fit covariance, so `phiStdDev` changes slightly.
- `sk-prepare` writes diagnostics to stderr, so stdout is always valid JSON.
- Version comes from a single source (`cpmg.__version__`); the banner shows the real citation.
- Requires Python >= 3.9.
- Code split into four modules (`models`, `fitting`, `analysis`, `report`) behind `cpmg.model_2state.CPMG_model`; matplotlib is imported lazily.
- kex/pB/phi uncertainties are now consistent with the fit covariance.
- The final stdout block is strict JSON (non-finite values become `null`).

### Fixed
- Crashes and empty results for slow+Meiboom and `exchange: "undefined"` configs.
- NaN/inf in edge cases (repeated eigenvalues, kex=0, dw=0, London coalescence).
- Numerical fallbacks for singular parameters instead of failures.
- Docs: `field` is documented as the observed-nucleus Larmor frequency; stdout layout described exactly.
- Demo evaluators tolerate unavailable values.

## 1.2.1

- Docs: English reference `USAGE.en.md`.
- Docs: `exchange: "undefined"` documented.
- Release workflow added; `dist/` is no longer tracked.
