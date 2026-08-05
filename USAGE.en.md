# ShereKhan 1.2.1 Reference

> 🌐 [한국어](USAGE.md) · **English**

NMR CPMG relaxation dispersion analysis tool — reference for scripts, config, and data formats.  
Author: Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)  
Original code: Adam Mazur, Bjoern Hammesfahr, Christian Griesinger, Donghan Lee, Martin Kollmar (Max-Planck-Institute for Biophysical Chemistry, 2012)  
Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)

> For the **decision and validation procedure** (individual vs global fit, AIC comparison,
> jackknife, heterogeneity diagnosis), see the step-by-step manual
> [`MANUAL.en.md`](MANUAL.en.md) · [한국어](MANUAL.ko.md).
> This document covers only the **format reference** for scripts, config, and data.

---

## Quick Start

```bash
# 1. Generate synthetic data (for testing)
python sk_createSyntheticDataset.py

# 2. Generate a config file
python sk_prepare.py synth-60.dat synth-90.dat > run.conf

# 3. Edit the config (exchange, model, residue flags, ...)
#    → edit run.conf directly

# 4. Run the analysis
python sk_run.py run.conf
```

---

## Scripts

### `sk_createSyntheticDataset.py` — Synthetic data generation

Generates a CPMG dataset for testing.

```bash
python sk_createSyntheticDataset.py
```

**Output files:**
- `synth-60.dat` — 60.12 MHz field data
- `synth-90.dat` — 90.23 MHz field data

**Parameters (edit inside the script):**

| Variable | Default | Description |
|------|--------|------|
| `fields` | `[60.12, 90.23]` | 1H field strengths (MHz) |
| `tcp` | `0.04` | CPMG pulse delay (s) |
| `kAB` | `20.0` | A→B exchange rate (s⁻¹) |
| `kBA` | `1000.0` | B→A exchange rate (s⁻¹) |
| `noiseRatio` | `0.02` | Noise level (2%) |
| `res` | (residue list) | Residues and their chemical shift differences (csd, ppm) |

---

### `sk_prepare.py` — Config file generation

Reads experimental data files (.dat) and generates an analysis config (JSON).

```bash
python sk_prepare.py <field1.dat> <field2.dat> [field3.dat ...] > myproject.conf
```

**Example:**
```bash
python sk_prepare.py exp-600MHz.dat exp-900MHz.dat > project.conf
```

**Output:** config in JSON format (stdout). Redirect to a file, then edit.

**Automatic detection:** computes each residue's alpha value to estimate the fast/slow exchange regime and suggests a suitable model.

---

### `sk_run.py` — Run the analysis

Reads a config file, fits the CPMG dispersion curves, and writes the results.

```bash
python sk_run.py <config.conf>
```

**Example:**
```bash
python sk_run.py project.conf
```

**Output files (named from `Project Name`):**
- `{name}.log` — text log of the fit results
- `{name}.pdf` — dispersion curve plots

**stdout:** fitted values for all residues, in JSON format

---

## Config file format (`.conf`)

```json
{
    "Project Name": "myproject",
    "comments": "",
    "datasets": [
        "exp-600MHz.dat",
        "exp-900MHz.dat"
    ],
    "exchange": "fast",
    "model": "Meiboom",
    "init": {
        "mode": "guess"
    },
    "residues": [
        {"name": "A10", "alpha": 1.44, "flag": "on"},
        {"name": "G15", "alpha": 0.83, "flag": "off"}
    ]
}
```

### Main fields

| Field | Values | Description |
|------|-----|------|
| `exchange` | `"fast"` / `"slow"` / `"undefined"` | Exchange regime. `"undefined"` is always paired with model `Matrix`. Ambiguous multi-field data (fast/slow tie) → every residue is automatically set to `flag: "off"`. Single field → also `"undefined"`/`Matrix`, but every residue gets `flag: "on"` (alpha needs ≥ 2 fields) |
| `model` | `"Meiboom"` / `"Matrix"` / `"London"` | Fitting model |
| `init.mode` | `"guess"` / `"values"` | How initial values are set |
| `residues[].flag` | `"on"` / `"off"` | Include this residue in the analysis |
| `compare_aic` | `true` / `false` | Compare global vs individual fit by AIC, plus a per-residue verdict (`sk_prepare` defaults to `true`; `false` when the key is absent) |
| `jackknife` | `true` / `false` | Validate the global fit by leave-one-residue-out jackknife (`sk_prepare` defaults to `true`; `false` when the key is absent) |

### Explicit initial values (`init.mode = "values"`)

**Matrix / London models:**
```json
"init": {
    "mode": "values",
    "kex": 1000.0,
    "pB": 0.02,
    "csd": 2.5
}
```

**Meiboom model:**
```json
"init": {
    "mode": "values",
    "kex": 1000.0,
    "phi": 0.5
}
```

---

## Model selection guide

| Exchange regime | Recommended model | Description |
|-----------|-----------|------|
| Fast | `Meiboom` | Fast-exchange approximation, few parameters (φ, kex) |
| Fast | `Matrix` | Exact matrix calculation (kAB, kBA, Δω) |
| Slow | `London` | Slow-exchange approximation |
| Slow | `Matrix` | Exact matrix calculation |

> `alpha < 1` → slow exchange likely  
> `alpha > 1` → fast exchange likely

---

## Decision and validation (AIC comparison · jackknife)

With `compare_aic: true` / `jackknife: true`, the global vs individual verdict and the
jackknife validation results are written to the `.log` file and to stdout. **For how to read
them and for the heterogeneity diagnosis procedure, see the step-by-step manual**
[`MANUAL.en.md`](MANUAL.en.md) §5–8.

Machine-readable stdout JSON — labelled blocks:

```
##### model_comparison
{ ... compareModelsAIC result (including per_residue) ... }
##### jackknife
{ ... jackknifeGlobal result ... }
#####
[ ... final per-residue fit values ... ]
```

- The `##### model_comparison` / `##### jackknife` blocks are strict JSON (non-finite values become `null`).
- **The lines after the last `#####` are the residue result list** (compatible with existing parsers, `split('#####')[-1]`).
- When `compare_aic`/`jackknife` are off, the corresponding block is not written.

---

## Data file format (`.dat`)

```
 60.12          ← field strength (MHz, 1H)
0.040000        ← TCP (s)
#nu_cpmg(Hz)        R2(1/s)      Esd(R2)
# A10           ← residue name (# comment)
  50.000   22.415    0.448
 100.000   21.893    0.438
 200.000   20.912    0.418
 ...
# G15
  50.000   18.234    0.365
 ...
```

---

## Full validation tests

```bash
bash test/run_tests.sh        # fast/slow exchange test cases
```

---

## Typical workflow

```
experimental data (.dat files)
        ↓
sk_prepare.py → draft config (.conf)
        ↓
edit config (exchange, model, flags, initial values)
        ↓
sk_run.py → run the fit
        ↓
inspect results (.log, .pdf)
        ↓
adjust config and re-run if needed
```
