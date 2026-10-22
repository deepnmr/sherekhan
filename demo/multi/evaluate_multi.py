#!/usr/bin/env python3
"""Evaluate the N-individual (Family A) and K-global-group (Family B) scenarios.

Family A: does the tool identify the N individual residues (rank by chi2_global),
          and is a clean global fit restored once they are removed?
Family B: across three structural hypotheses -- 1 global fit, K global fits (one
          per ground-truth group), or fully individual -- which does AICc pick?
          The K-group model should win for K-group data; each group's own
          jackknife SE should be small (internally consistent global fit).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
from cpmg.model_2state import CPMG_model

FAMILY_A = ['m1', 'm2', 'm3', 'm4']
FAMILY_B = ['g2', 'g3']


def fmt(v, spec):
    """Format v with a format spec such as '.0f'; None -> 'n/a'."""
    return 'n/a' if v is None else format(v, spec)


def pct(v):
    return 'n/a' if v is None else '%.1f%%' % v


def relse(jk):
    """Jackknife relative SE in percent, or None when any input is null."""
    if not jk or jk['se'] is None or jk['kex_full'] is None:
        return None
    return jk['se'] / abs(jk['kex_full']) * 100.0


def build_model(name):
    d = os.path.join(HERE, name)
    truth = json.load(open(os.path.join(d, 'truth.json')))
    m = CPMG_model()
    m.verbose = False
    m.dataset.addField(os.path.join(d, 'f60.dat'))
    m.dataset.addField(os.path.join(d, 'f90.dat'))
    m.dataset.exchange = 'slow'
    m.model = 'Matrix'
    m.init_dd(2.5)
    m.kAB, m.kBA = 0.1 * 1000.0, 0.9 * 1000.0
    m.init_R2_0()
    return m, truth


def set_active(m, labels):
    for r in m.dataset.rsds:
        r.active = r.label in labels
    m._cache_arrays()


def global_score(m):
    """Fit the currently-active set as one global fit; return (chi2, k, n)."""
    m.fit()
    p = m._build_p0_current()
    r = m.errFunc(p)
    return float(np.dot(r, r)), len(p), len(r)


def eval_family_A(name):
    m, truth = build_model(name)
    individuals = set(truth['individuals'])
    N = len(individuals)
    alllabels = [r.label for r in m.dataset.rsds]

    set_active(m, alllabels)
    m.fit()
    cmp = m.compareModelsAIC()

    # Identification: rank residues by their chi2 contribution under global.
    # The top-N needs the complete ranking: any null chi2_global skips it.
    na = [d['label'] for d in cmp['per_residue'] if d['chi2_global'] is None]
    ranked = sorted((d for d in cmp['per_residue'] if d['chi2_global'] is not None),
                    reverse=True, key=lambda d: d['chi2_global'])
    topN = {d['label'] for d in ranked[:N]}
    hit = len(topN & individuals)

    # Recovery: refit with the true individuals removed.
    shared = [l for l in alllabels if l not in individuals]
    set_active(m, shared)
    m.fit()
    rec = m.compareModelsAIC()
    rjk = m.jackknifeGlobal()
    rel = relse(rjk)

    print('-' * 84)
    print('%s  (Family A: %d individual%s + %d shared)   true individuals: %s' % (
        name.upper(), N, '' if N == 1 else 's', len(shared), sorted(individuals)))
    print('  whole-dataset preferred : %s   (prefer global=%d, individual=%d)' % (
        cmp['preferred'], cmp['n_pref_global'], cmp['n_pref_individual']))
    # Aggregation rule: a check is SKIPPED when any value it needs is None.
    # Scenario = FAIL if any check FAILs; else PASS if any check PASSes; else SKIPPED.
    checks = []
    if na:
        checks.append('SKIPPED')
        print('  identification (top-%d by chi2_global): skipped: n/a (chi2_global unavailable for %s)' % (
            N, ', '.join(na)))
    else:
        checks.append('PASS' if hit == N else 'FAIL')
        print('  identification (top-%d by chi2_global): %s   -> %d/%d individuals found' % (
            N, sorted(topN), hit, N))
    print('  ranked chi2_global: ' + '  '.join(
        '%s=%s%s' % (d['label'], fmt(d['chi2_global'], '.0f'), '*' if d['label'] in individuals else '')
        for d in ranked))
    print('  recovery (individuals removed): whole=%s  jackknife relSE=%s  kex=%s' % (
        fmt(rec['preferred'], 's'), pct(rel), fmt(rjk['kex_full'] if rjk else None, '.0f')))
    # Recovery is a separate check; it reads rec['preferred'].
    if rec['preferred'] is None:
        checks.append('SKIPPED')
        print('  recovery RESULT: skipped: n/a (preferred unavailable for whole dataset)')
    else:
        checks.append('PASS' if rec['preferred'] == 'global' else 'FAIL')
        print('  recovery RESULT: %s  (global restored after removal)' % checks[-1])
    if 'FAIL' in checks:
        label = 'FAIL'
    elif 'PASS' in checks:
        label = 'PASS'
    else:
        label = 'SKIPPED'
    print('  RESULT: %s  (%d check(s) skipped)' % (label, checks.count('SKIPPED')))
    return name, label, label != 'FAIL'


def eval_family_B(name):
    m, truth = build_model(name)
    groups = {}
    for label, gid in truth['group'].items():
        groups.setdefault(gid, []).append(label)
    K = len(groups)
    alllabels = [r.label for r in m.dataset.rsds]

    # Hypothesis 1: one global fit for everyone.
    set_active(m, alllabels)
    chi2_1, k_1, n = global_score(m)
    cmp = m.compareModelsAIC()          # also gives the individual aggregate
    chi2_ind, k_ind = cmp['individual']['chi2'], cmp['individual']['k']

    # Hypothesis K: one global fit per ground-truth group.
    chi2_K, k_K = 0.0, 0
    group_relse = {}
    for gid, labels in groups.items():
        set_active(m, labels)
        c, k, _ = global_score(m)
        chi2_K += c
        k_K += k
        jk = m.jackknifeGlobal()
        group_relse[gid] = relse(jk)

    set_active(m, alllabels)            # restore

    kg_label = '%d-global' % K
    # Selection needs the individual chi2 and k; a null there skips the check.
    have_ind = chi2_ind is not None and k_ind is not None
    aicc = {'1-global': m._aic(chi2_1, k_1, n)[1], kg_label: m._aic(chi2_K, k_K, n)[1]}
    if have_ind:
        aicc['individual'] = m._aic(chi2_ind, k_ind, n)[1]
    winner = min(aicc, key=aicc.get)

    print('-' * 84)
    print('%s  (Family B: %d global groups)   groups: %s' % (
        name.upper(), K, {g: groups[g] for g in groups}))
    print('  single global fit on all -> whole-dataset preferred: %s' % cmp['preferred'])
    print('  structure model selection (lower AICc wins):')
    amin = min(aicc.values())
    for label in aicc:
        print('    %-12s chi2=%9.2f  AICc=%10.2f  dAICc=%+9.2f %s' % (
            label,
            {'1-global': chi2_1, kg_label: chi2_K, 'individual': chi2_ind}[label],
            aicc[label], aicc[label] - amin, '  <-- WINNER' if label == winner else ''))
    print('  per-group internal jackknife relSE: ' + '  '.join(
        '%s=%s' % (g, pct(group_relse[g])) for g in groups))
    if not have_ind:
        label = 'SKIPPED'
        print('  RESULT: skipped: n/a (individual chi2/k unavailable)')
    else:
        label = 'PASS' if winner == kg_label else 'FAIL'
        print('  RESULT: %s  (%d-global model wins as expected)' % (label, K))
    return name, label, label != 'FAIL'


if __name__ == '__main__':
    print('=' * 84)
    print('FAMILY A -- majority global + N individuals (does the tool find them?)')
    print('=' * 84)
    ra = [eval_family_A(n) for n in FAMILY_A]
    print()
    print('=' * 84)
    print('FAMILY B -- data best described by K separate global fits')
    print('=' * 84)
    rb = [eval_family_B(n) for n in FAMILY_B]

    print()
    print('=' * 84)
    print('SUMMARY')
    for name, label, _ in ra + rb:
        print('  %-4s %s' % (name, label))
    raise SystemExit(1 if any(label == 'FAIL' for _, label, _ in ra + rb) else 0)
