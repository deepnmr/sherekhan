# ShereKhan 1.2.0

NMR CPMG relaxation-dispersion curve analysis — two-state chemical-exchange fitting.

- **Author:** Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
- **Original code:** Adam Mazur, Bjoern Hammesfahr, Christian Griesinger,
  Donghan Lee, Martin Kollmar (Max-Planck-Institute for Biophysical Chemistry, 2012)
- **License:** MIT (see `LICENSE`)

## Install

```bash
pip install .
```

This installs the `cpmg` package and three console scripts:

| Command      | Module                        | Purpose                          |
|--------------|-------------------------------|----------------------------------|
| `sk-synth`   | `sk_createSyntheticDataset`   | Generate synthetic test data     |
| `sk-prepare` | `sk_prepare`                  | Build a JSON config from `.dat`  |
| `sk-run`     | `sk_run`                      | Run the fit, write `.log`/`.pdf` |

## Quick start

```bash
sk-synth                                      # → synth-60.dat, synth-90.dat
sk-prepare synth-60.dat synth-90.dat > run.conf
sk-run run.conf                               # → run.log, run.pdf, JSON on stdout
```

Run without installing:

```bash
python sk_createSyntheticDataset.py
python sk_prepare.py synth-60.dat synth-90.dat > run.conf
python sk_run.py run.conf
```

See [`USAGE.md`](USAGE.md) for full documentation (config format, models, AIC comparison).

## Models

| Model     | Regime      | Description                                        |
|-----------|-------------|----------------------------------------------------|
| `Matrix`  | any         | Exact Liouvillian matrix-exponential solution      |
| `London`  | slow        | Analytical slow-exchange approximation             |
| `Meiboom` | fast        | Analytical fast-exchange approximation (fits φ)    |

## Requirements

Python ≥ 3.8, `numpy>=1.20`, `scipy>=1.6`, `matplotlib>=3.3`.

## Tests

```bash
bash test/run_tests.sh
```
