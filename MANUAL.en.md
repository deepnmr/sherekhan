# ShereKhan Step-by-Step Manual — Individual vs Global fit + Jackknife validation

> 🌐 [한국어](MANUAL.ko.md) · **English**

A practical procedure for deciding, from CPMG relaxation-dispersion data, whether
**all residues share a single exchange process (global fit)** or **each residue has
its own kex (individual fit)**, and for validating the robustness of the global fit
with the **jackknife**.

- Audience: users of the `sk-prepare` → `sk-run` pipeline
- Prerequisites: `README.md` (install), `USAGE.md` (config/model details)
- This document covers the **order of steps** from raw `.dat` to verdict and validation.

---

## 0. Install (once)

```bash
pip install .
```

Installs three console scripts: `sk-synth`, `sk-prepare`, `sk-run`.
To run without installing, use `python sk_prepare.py ...`, `python sk_run.py ...`.

Requirements: Python ≥ 3.8, `numpy`, `scipy`, `matplotlib`.

---

## 1. Prepare input data (`.dat`)

One file per magnetic field. Format:

```
 60.12                        # field (MHz, 1H Larmor)
0.040000                      # tcp (s, total CPMG echo time)
#nu_cpmg(Hz)  R2(1/s)  Esd(R2)
# R1                          # residue label
  50.000  35.80  0.36
 100.000  35.14  0.35
 ...
# R2
 ...
```

- The same residue must have the **same label in every file** to be linked across fields.
- To determine kex well you need **two or more fields** (jackknife/comparison are only
  meaningful at field ≥ 2).
- No test data? Generate synthetic data with `sk-synth`.

---

## 2. Step 1 — build a config (`sk-prepare`)

```bash
sk-prepare f60.dat f90.dat > run.conf
```

What `sk-prepare` does:
- Computes per-residue alpha → estimates fast/slow exchange regime, suggests a model
  (fast→Meiboom, slow→London, ambiguous→Matrix)
- **Sets `compare_aic: true` and `jackknife: true` by default** (per-residue comparison
  + jackknife enabled)

Output is JSON on stdout. Save to a file and review it in the next step.

---

## 3. Step 2 — review/edit the config (`run.conf`)

Key fields:

| Field | Value | Meaning |
|-------|-------|---------|
| `exchange` | `"fast"` / `"slow"` / `"undefined"` | exchange regime, always paired with model `Matrix` when `"undefined"`. **Ambiguous** multi-field data (equal fast/slow tally) → every residue auto-set `flag: "off"` (turn the ones you want on manually). **Single-field** data → also `"undefined"`/`Matrix` but every residue is left `flag: "on"` (alpha needs field ≥ 2, so nothing is disabled) |
| `model` | `"Meiboom"` / `"London"` / `"Matrix"` | fitting model (Matrix = exact, slow) |
| `init.mode` | `"guess"` / `"values"` | initial values: grid search / explicit |
| `residues[].flag` | `"on"` / `"off"` | include the residue or not |
| `compare_aic` | `true` / `false` | per-residue individual-vs-global comparison |
| `jackknife` | `true` / `false` | leave-one-residue-out validation |

Review checklist:
- Confirm residues to analyze are `flag: "on"`; set obvious noise residues to `"off"`.
- Confirm `compare_aic`/`jackknife` are on (the key verdict/validation output).
- If you know initial values, use `init.mode: "values"` to skip the grid search. Example:
  ```json
  "init": {"mode": "values", "kex": 1000.0, "pB": 0.02, "csd": 2.5}
  ```
  (For Meiboom specify `"phi"` instead of `pB`/`csd`.)

---

## 4. Step 3 — run (`sk-run`)

```bash
sk-run run.conf
```

Outputs:
- `<Project Name>.log` — human-readable full report (global fit stats + comparison table + jackknife)
- `<Project Name>.pdf` — per-residue dispersion-curve plots (exp points + calc curves)
- **stdout** — JSON blocks (for programmatic use, §9)

---

## 5. Step 4 — read the global fit result (`.log`)

The statistics block in `.log`:

```
npar=13 nvar=66 ndof=53 chi2=  63.629 chi2/dof=   1.201
kex: 1051.054  +-    94.327
 pB:    0.020  +-     0.002
resId:  csd [ppm]   R2_0 (60.1 MHz)   R2_0 (90.2 MHz)
  R1:   7.669 ...
```

- `chi2/dof ≈ 1` means the global fit is consistent with the data and errors. `≫ 1` means
  model misfit or underestimated errors.
- `kex ± std`, `pB ± std` are the shared exchange parameters. Per-residue `csd` (or phi)
  and `R2_0` are the local parameters.

This block is the **global fit** result. The next step compares it to the individual fit.

---

## 6. Step 5 — read the individual-vs-global verdict (core)

At the bottom of `.log` / on stdout, the comparison block:

```
Model               chi2      k      n          AIC         AICc
global            86.139     16    110      118.139      123.989
individual        83.084     20    110      123.084      132.523
delta AICc (individual - global):        8.534
Preferred model (whole dataset): global
---
Per-residue preference (individual vs global fit):
 resId    chi2_glob   chi2_indiv        dAICc       z(kex)     better
   R1        23.31        21.12          4.23         1.55     global
   R5      1164.47        11.38      -1146.68        26.31 individual
---
Residues preferring global: 4   preferring individual: 1
```

**Whole-dataset verdict** — `Preferred model`:
- `delta AICc > 0` → **global** (the shared kex is better per parameter)
- `delta AICc < 0` → **individual** (each residue needs its own kex)

**Per-residue verdict** — the `better` column:
- `chi2_glob` = that residue's chi² contribution under the shared kex
- `chi2_indiv` = its chi² when fitted alone
- `dAICc` = AICc(individual) − AICc(global): **> 0 → global**, **< 0 → individual**
  - Convention: individual pays for its own exchange parameter(s); global borrows the
    shared kex for free → "is this residue's own kex worth its cost?"
- `z(kex)` = (kex_individual − kex_leave-one-out) / σ(kex_individual)
  - Referenced to the **leave-one-out** kex (fit from the other residues), so it is unbiased.
  - **|z| > 2** → the residue statistically disagrees with the shared kex (outlier candidate).

Reading order: ① the whole-dataset verdict for the big picture → ② the `better` column for
which residues split off → ③ the size of `z(kex)`/`chi2_glob` for the split residues to
gauge severity.

---

## 7. Step 6 — validate the global fit with the jackknife

```
Jackknife validation of the global fit (leave-one-residue-out)
kex (full global fit):            1134.895
kex (jackknife mean):             1133.920
jackknife std error (kex):          67.180
relative std error:                  5.92%
---
 resId      kex(-res)      delta_kex        |z|
   R1        1166.713        -31.818       0.47
   R4        1179.261        -44.366       0.66
---
(* = influential: dropping this residue shifts kex by > 2 jackknife SE)
```

Drop each residue in turn and refit globally → track the change in the shared kex.

- **Small relative SE** (e.g. < a few %) + no `*`-flagged residue → **robust global fit**; no
  single residue dominates kex.
- **Large relative SE** (e.g. > 15–20%) → unstable kex; a sign of kex disagreement among residues.
- `delta_kex` = kex_full − kex(that residue removed). `|delta_kex| > 2 SE` is flagged `*` →
  influential residue.

> Note: a single strong outlier inflates the SE itself, which can **mask** the `*` flag
> (single-deletion jackknife breakdown point = 1). In that case `z(kex)`/`chi2_glob` in §6 are
> more sensitive. Relative SE remains a valid overall-stability indicator.

---

## 8. Step 7 — heterogeneity diagnosis: individual residues (N) & group structure (K)

A whole-dataset verdict of individual is not the end. Dig into **why** in two ways.

> **Key principle**: least-squares global fitting is **not robust** to outliers. Even one
> divergent residue can dominate the total chi², drag the shared kex, and then even the
> normal residues no longer match the contaminated global kex — so **they all flip to
> `better = individual`**. The binary `preferred` flag alone hides the cause. The **magnitude
> (ranking) of `chi2_glob`** is the signal that survives contamination.

### 8-A. Diagnose individual residues (N) → exclude → re-verify

**① Pinpoint (chi2_glob ranking)** — sort the comparison table by `chi2_glob` descending: the
true individual residues cluster at the top. Even with several individuals (N of them), the top
N are those residues, and the ranking holds even under heavy contamination (e.g. 4 of 6 are
individuals). `z(kex)` is a secondary signal.

```
ranked chi2_global:  R6=1241*  R3=659*  R4=298*  R5=288*  R2=82  R1=56   (* = true individual)
```
→ the top 4 (R6, R3, R4, R5) are individuals, clearly separated from the normal residues (R1, R2).

**② Exclude and refit** — set the pinpointed residues' `flag` to `"off"`:
```bash
# in run.conf, set the individual residues' "flag": "on" -> "off"
sk-run run.conf
```

**③ Re-verify** — if the remaining residues truly shared one kex:
- the whole-dataset verdict returns to **global**, all remaining residues `better = global`
- jackknife relative SE drops sharply (e.g. 12.6% → 0.6%; ≤ ~3% even after removing N residues)
- `kex` converges to the true consensus value

### 8-B. Detect multiple global groups (K)

If exclusion does not restore global and **the fitted kex values cluster into a few groups**,
the data may consist not of one global fit but of **K distinct global fits** (e.g. residues
split into a kex≈600 group and a kex≈1800 group). Procedure:

**① Group hypothesis** — from the per-residue individual fit kex (the `fit_kex`/individual
`kex` in the comparison table), split residues into K groups by clustering value.

**② Per-group global fit** — build a config with only one group `flag: on` and run each separately:
```bash
# groupA.conf : only group A residues on   /   groupB.conf : only group B residues on
sk-run groupA.conf
sk-run groupB.conf
```
If each group has **all `better = global`** + **small jackknife relative SE**, that group is an
internally consistent global fit.

**③ Structure selection (AICc)** — compare the total AICc of three hypotheses (same number of
data points n):

| Structure | chi² | parameters k |
|-----------|------|--------------|
| 1-global | whole global fit chi² (large) | 1 kex + locals |
| **K-global** | Σ(per-group global chi²) | K kex + locals |
| individual | Σ(per-residue individual chi²) | N kex + locals |

`AICc = chi² + 2k + 2k(k+1)/(n−k−1)`, **lowest wins**. Here `n` is the **total number of data
points** (summed over all residues, not per-group n) — the three hypotheses describe the same
data, so n is the same. Reading each group config's `global` `chi2`/`k` from its `Model
comparison` table and summing gives the K-global chi²/k. For K-group data, K-global beats both
1-global (underfit) and individual (overfit) — individual's chi² is slightly lower but pays for
N exchange-parameter values.

```
1-global     chi2=  2955.88  AICc=  3003.45
2-global     chi2=    91.69  AICc=   144.97   <-- WINNER
individual   chi2=    87.57  AICc=   165.99
```
→ "this data contains 2 global fits" is statistically confirmed.

> The 8-A / 8-B procedures are reproduced automatically in `demo/` (single outlier) and
> `demo/multi/` (N individuals, K groups) — see §10.

---

## 9. stdout JSON integration (programmatic processing)

Labelled JSON blocks on `sk-run` stdout:

```
##### model_comparison
{ ... compareModelsAIC result (includes per_residue) ... }
##### jackknife
{ ... jackknifeGlobal result ... }
#####
[ ... per-residue final fit values (reportAllValues) ... ]
```

- The `##### model_comparison` / `##### jackknife` blocks are **strict JSON** (non-finite → `null`).
- **The line after the last `#####` is the residue-results list** — always last, for backward
  compatibility. Read "the line after the last `#####`" to be safe (`split('#####')[-1]`).
- If `compare_aic`/`jackknife` are off, those blocks are not emitted.

---

## 10. Quick validation (optional)

Confirm the feature judges known ground truth correctly:

```bash
bash demo/run_demo.sh          # global / individual / mixed (single outlier)
bash demo/multi/run_multi.sh   # N individuals (m1..m4), K global groups (g2, g3)
```

- `demo/run_demo.sh` — generates, runs, and evaluates the global/individual/mixed scenarios
  (exit 0 = all pass).
- `demo/multi/run_multi.sh` — reproduces 8-A (identify/recover 1–4 individuals) and 8-B (AICc
  structure selection for 2/3 global groups).
- See `demo/README.md`, `demo/multi/README.md` for the scenarios and expected results.

---

## 11. Decision-rule summary

```
Whole-dataset verdict:
  delta AICc (individual - global) > 0  ->  GLOBAL
                                    < 0  ->  INDIVIDUAL

Per-residue (better column):
  dAICc > 0  ->  global      |  dAICc < 0  ->  individual
  |z(kex)| > 2               ->  disagrees with the shared kex (outlier candidate)

Jackknife:
  small relative SE and no *  ->  robust global fit
  *-flagged residue           ->  influential residue

Heterogeneity diagnosis (when the whole dataset is individual):
  top N by chi2_glob ranking  ->  individual residues (robust to contamination; §8-A)
    -> flag off -> refit -> confirm global is restored
  fitted kex clusters into K  ->  per-group fit + AICc(1-global vs K-global vs individual)
    -> if the lowest AICc is K-global, "K global fits" structure is confirmed (§8-B)
```

---

## 12. Troubleshooting

| Symptom | Cause / action |
|---------|----------------|
| `Two or more fields strengths are necessary for alpha calculation` | Only 1 field → no comparison/alpha. Need field ≥ 2 data. |
| `Jackknife skipped: need at least 2 active residues.` | Fewer than 2 active residues. Jackknife needs ≥ 2. |
| `covariance matrix is None` warning | Singular / under-constrained fit. Stds shown as 0. Check initial values, model, data. |
| `n is not an integer` | tcp × nu_CPMG is not an integer. Check the tcp/frequencies in the `.dat`. |
| Matrix model is slow | Jackknife does N refits. For big data use Meiboom/London or `init.mode: values` to speed up. |
| Whole dataset is individual, cause unclear | Use the §8 procedure: suspect residues with large `z(kex)`/`chi2_glob` as outliers → exclude and re-verify. |
