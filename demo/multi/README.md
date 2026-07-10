# Structured heterogeneity demo — N individuals & K global groups

Extends the single-outlier `demo/` to two systematic families (6 residues,
2 fields, exact Matrix model, `pB=0.1`, 1% noise, known ground truth).

## Family A — one shared kex (1000) + N individual residues

| Scenario | # individual | Question |
|----------|--------------|----------|
| `m1`..`m4` | 1, 2, 3, 4 | Does the tool *identify* which residues are individual, and is a clean global fit restored once they are removed? |

Result: the single-`kex` verdict flips to **individual** for every N (least-squares
is not robust). But ranking residues by **`chi2_global`** identifies the true
individuals with 100% precision at every N (even N=4, where 4 of 6 residues are
individuals and heavily contaminate the shared kex). Removing them restores a
clean global fit (jackknife relative SE ≤ ~3%, kex ≈ 1000). Takeaway: the binary
`preferred` flag saturates under contamination; the **`chi2_global` magnitude**
is the robust individual detector.

## Family B — data best described by K separate global fits

| Scenario | Groups (kex) | Question |
|----------|--------------|----------|
| `g2` | {600} \| {1800} | Is a **2-global** structure selected over 1-global and fully-individual? |
| `g3` | {400} \| {1000} \| {2500} | Is a **3-global** structure selected? |

Result: one global fit is strongly rejected (chi² explodes). AICc model selection
across **1-global vs K-global vs individual** picks the **K-global** structure —
it beats 1-global (underfit) and fully-individual (overfit: individual's chi² is
marginally lower but pays for N exchange params). Each group's internal jackknife
relative SE is small (≤ ~2%), confirming each is a robust global fit. I.e. "the
data contains K global fits" is statistically confirmed.

## Run

```bash
bash demo/multi/run_multi.sh          # generate -> fit -> evaluate (exit 0 = all pass)
PYTHON=/path/to/python bash demo/multi/run_multi.sh
```

Generated `<scenario>/` outputs are reproducible and git-ignored.
