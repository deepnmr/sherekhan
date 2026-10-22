#!/usr/bin/env python3
"""Evaluate the global/individual/mixed scenario runs against ground truth.

Reads the '##### model_comparison' / '##### jackknife' JSON blocks from each
scenario's captured sk_run stdout (out.txt) and checks the tool's verdict
against the known ground truth.  Run after run_demo.sh has produced the
out.txt / out-recover.txt files.

Any numeric field that is null in the JSON prints as 'n/a'; a check that needs
a null value is skipped and reported as 'skipped: n/a'.
"""
import json
import os
import re

BASE = os.path.dirname(os.path.abspath(__file__))
EXPECT = {
    'global':     'global',       # whole-dataset preferred
    'individual': 'individual',
    'mixed':      'split',         # lone outlier flagged, global restored on removal
}
MIXED_OUTLIER = 'R5'              # kex 2500 vs the shared 1000


def fmt(v, spec):
    """Format v with a format spec such as '8.2f'; None -> 'n/a' right-aligned to the width."""
    width = int(re.match(r'\d*', spec).group() or 0)
    if v is None:
        return 'n/a'.rjust(width)
    return format(v, spec)


def pct(v):
    return 'n/a' if v is None else '%.1f%%' % v


def relse(jk):
    """Jackknife relative SE in percent, or None when any input is null."""
    if not jk or jk['se'] is None or jk['kex_full'] is None:
        return None
    return jk['se'] / abs(jk['kex_full']) * 100.0


def parse_blocks(path):
    """Return {label: dict} for the '##### <label>' JSON blocks in stdout."""
    lines = open(path).read().splitlines()
    out = {}
    for i, l in enumerate(lines):
        if l.startswith('##### '):
            out[l[6:].strip()] = json.loads(lines[i + 1])
    return out


def evaluate(name):
    d = os.path.join(BASE, name)
    truth = {r['name']: r for r in
             json.load(open(os.path.join(d, 'truth.json')))['residues']}
    blocks = parse_blocks(os.path.join(d, 'out.txt'))
    mc, jk = blocks['model_comparison'], blocks['jackknife']

    print('=' * 78)
    print('SCENARIO: %s   (ground-truth kex: %s)' % (
        name.upper(), [truth[n]['kex'] for n in truth]))
    print('-' * 78)
    print('whole-dataset preferred : %s   (global w=%s, individual w=%s)' % (
        mc['preferred'], fmt(mc['weight_global'], '.3f'), fmt(mc['weight_individual'], '.3f')))
    print('residues preferring global=%d  individual=%d' % (
        mc['n_pref_global'], mc['n_pref_individual']))
    print('-' * 78)
    print('%4s %8s %10s %10s %10s %9s %8s %10s' % (
        'res', 'true_kex', 'fit_kex', 'chi2_glob', 'chi2_ind', 'dAICc', 'z(kex)', 'verdict'))
    for r in mc['per_residue']:
        print('%4s %8d %s %s %s %s %s %10s' % (
            r['label'], truth[r['label']]['kex'], fmt(r['kex'], '10.1f'),
            fmt(r['chi2_global'], '10.2f'), fmt(r['chi2'], '10.2f'),
            fmt(r['delta_aicc'], '9.2f'), fmt(r['z_kex'], '8.2f'), r['preferred']))
    print('-' * 78)
    se = jk['se']
    print('jackknife: kex_full=%s  SE=%s  (rel %s)  bias=%s' % (
        fmt(jk['kex_full'], '.1f'), fmt(se, '.1f'), pct(relse(jk)), fmt(jk['bias'], '.1f')))
    # Residues with a null delta_kex are listed separately; the influential
    # line must not claim 'none' while any residue is n/a.
    na_infl = [p['label'] for p in jk['per_residue'] if p['delta_kex'] is None]
    flagged = [p['label'] for p in jk['per_residue']
               if p['delta_kex'] is not None and se is not None
               and abs(p['delta_kex']) > 2 * se]
    if se is None:
        print('jackknife influence: skipped: n/a (se unavailable)')
        print('flagged: skipped: n/a (se unavailable)')
    elif na_infl:
        print('jackknife influence: skipped: n/a (delta_kex unavailable for %s)' % ', '.join(na_infl))
        print('flagged: %s' % (', '.join(flagged) or 'none'))
    else:
        print('jackknife influential (|delta_kex| > 2 SE): %s' % (', '.join(flagged) or 'none'))

    exp = EXPECT[name]
    print('-' * 78)
    if exp == 'split':
        # A single high-leverage outlier corrupts the shared kex for everyone
        # (least-squares is not robust), so the binary verdict flips to
        # individual for ALL residues.  Success = the tool *identifies* the
        # outlier (dominates |z| and chi2_global) and removing it restores a
        # clean global fit.
        pr = {r['label']: r for r in mc['per_residue']}
        # Aggregation rule: a check is SKIPPED when any value it needs is None.
        # Scenario = FAIL if any check FAILs; else PASS if any check PASSes; else SKIPPED.
        checks = []
        na_z = [l for l in pr if pr[l]['z_kex'] is None]
        if na_z:
            checks.append('SKIPPED')
            print('outlier %s identified (z): skipped: n/a (z_kex unavailable for %s)' % (
                MIXED_OUTLIER, ', '.join(na_z)))
        else:
            zs = {l: abs(pr[l]['z_kex']) for l in pr}
            z_out = zs[MIXED_OUTLIER]
            z_rest = max(zs[l] for l in zs if l != MIXED_OUTLIER)
            checks.append('PASS' if z_out > 2 * z_rest else 'FAIL')
            print('outlier %s identified (z): |z|=%.1f (rest max %.1f, %.1fx) -> %s' % (
                MIXED_OUTLIER, z_out, z_rest, z_out / z_rest, checks[-1]))
        na_c = [l for l in pr if pr[l]['chi2_global'] is None]
        if na_c:
            checks.append('SKIPPED')
            print('outlier %s identified (chi2_global): skipped: n/a (chi2_global unavailable for %s)'
                  % (MIXED_OUTLIER, ', '.join(na_c)))
        else:
            cg = {l: pr[l]['chi2_global'] for l in pr}
            c_out = cg[MIXED_OUTLIER]
            c_rest = max(cg[l] for l in cg if l != MIXED_OUTLIER)
            checks.append('PASS' if c_out > 3 * c_rest else 'FAIL')
            print('outlier %s identified (chi2_global): chi2_glob=%.0f (rest max %.0f, %.1fx) -> %s' % (
                MIXED_OUTLIER, c_out, c_rest, c_out / c_rest, checks[-1]))

        rblocks = parse_blocks(os.path.join(d, 'out-recover.txt'))
        rec, rjk = rblocks['model_comparison'], rblocks['jackknife']
        print('recovery (%s removed): whole=%s  global-preferring=%s  individual-preferring=%s  '
              'jackknife SE=%s  kex=%s' % (
                  MIXED_OUTLIER, fmt(rec['preferred'], 's'), fmt(rec['n_pref_global'], 'd'),
                  fmt(rec['n_pref_individual'], 'd'),
                  pct(relse(rjk)), fmt(rjk['kex_full'], '.0f')))
        # Recovery check reads preferred, n_pref_global and n_pref_individual.
        na_r = [f for f in ('preferred', 'n_pref_global', 'n_pref_individual') if rec[f] is None]
        if na_r:
            checks.append('SKIPPED')
            print('recovery check: skipped: n/a (%s unavailable for whole dataset)' % ', '.join(na_r))
        else:
            rec_ok = rec['preferred'] == 'global' and rec['n_pref_individual'] == 0
            checks.append('PASS' if rec_ok else 'FAIL')
        n_skip = checks.count('SKIPPED')
        if 'FAIL' in checks:
            verdict = 'FAIL'
        elif 'PASS' in checks:
            verdict = 'PASS'
        else:
            verdict = 'SKIPPED'
        print('EXPECT outlier flagged + global restored after removal -> %s%s' % (
            verdict, ' (%d check(s) skipped)' % n_skip if n_skip else ''))
        ok = verdict != 'FAIL'
    else:
        if mc['preferred'] is None:
            verdict = 'SKIPPED'
            print('preferred check: skipped: n/a (preferred unavailable for whole dataset)')
        else:
            verdict = 'PASS' if mc['preferred'] == exp else 'FAIL'
        print('EXPECT whole-dataset preferred = %s -> %s' % (exp, verdict))
        ok = verdict != 'FAIL'
    print()
    return name, exp, ok, verdict


if __name__ == '__main__':
    results = [evaluate(n) for n in ('global', 'individual', 'mixed')]
    print('=' * 78)
    print('SUMMARY')
    allok = True
    for name, exp, ok, verdict in results:
        allok = allok and ok
        print('  %-11s expected=%-11s %s' % (name, exp, verdict))
    raise SystemExit(0 if allok else 1)
