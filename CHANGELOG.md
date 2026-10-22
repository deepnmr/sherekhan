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
- CI runs the tests on Python 3.9 through 3.14 and adds a build job; the release workflow runs only after the tests pass.
- `sk-run` checks that the output directory exists (and is a directory) before fitting.

### Changed
- Derived uncertainties propagate the fit covariance, so `phiStdDev` changes slightly.
- `sk-prepare` writes diagnostics to stderr, so stdout is always valid JSON.
- Version comes from a single source (`cpmg.__version__`); the banner shows the real citation.
- Packaging: `setuptools>=77`, version read dynamically from `cpmg.__version__`, `requires-python >= 3.9`, and a `[test]` extra.
- One echo-count rule (round tcp*nu, tolerance 1e-3, at least one cycle) shared by the scalar, eigenvector and batch Matrix paths, raising ValueError instead of exiting.
- Library hygiene: the `cpmg` package no longer calls `exit`, reads `sys.argv`, uses `np.matrix` or calls `os.uname`.
- Usage errors are a single line on stderr.
- `.dat` validation: finite positive field, tcp and nu_cpmg, finite R2, sigma > 0, label lines `# <one token>`; residues with a label but no data are dropped.
- Code split into four modules (`models`, `fitting`, `analysis`, `report`) behind `cpmg.model_2state.CPMG_model`; matplotlib is imported lazily.
- kex/pB/phi uncertainties are now consistent with the fit covariance.
- The final stdout block is strict JSON (non-finite values become `null`).

### Fixed
- Crashes and empty results for slow+Meiboom and `exchange: "undefined"` configs.
- NaN/inf in edge cases (repeated eigenvalues, kex=0, dw=0, London coalescence).
- Numerical fallbacks for singular parameters instead of failures.
- Docs: `field` is documented as the observed-nucleus Larmor frequency; stdout layout described exactly.
- The alpha (field-dependence) value is computed from the field matching each dataset.
- Demo evaluators tolerate unavailable values.

## 1.2.1

- Docs: English reference `USAGE.en.md`.
- Docs: `exchange: "undefined"` documented.
- Release workflow added; `dist/` is no longer tracked.
