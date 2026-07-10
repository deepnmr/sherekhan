#!/usr/bin/env python3
"""Evaluate the global/individual/mixed scenario runs against ground truth.

Reads the '##### model_comparison' / '##### jackknife' JSON blocks from each
scenario's captured sk_run stdout (out.txt) and checks the tool's verdict
against the known ground truth.  Run after run_demo.sh has produced the
out.txt / out-recover.txt files.
"""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))
EXPECT = {
    'global':     'global',       # whole-dataset preferred
    'individual': 'individual',
    'mixed':      'split',         # lone outlier flagged, global restored on removal
}
MIXED_OUTLIER = 'R5'              # kex 2500 vs the shared 1000


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
    print('whole-dataset preferred : %s   (global w=%.3f, individual w=%.3f)' % (
        mc['preferred'], mc['weight_global'], mc['weight_individual']))
    print('residues preferring global=%d  individual=%d' % (
        mc['n_pref_global'], mc['n_pref_individual']))
    print('-' * 78)
    print('%4s %8s %10s %10s %10s %9s %8s %10s' % (
        'res', 'true_kex', 'fit_kex', 'chi2_glob', 'chi2_ind', 'dAICc', 'z(kex)', 'verdict'))
    for r in mc['per_residue']:
        z = r['z_kex']
        zs = '%8.2f' % z if z is not None else '     n/a'
        print('%4s %8d %10.1f %10.2f %10.2f %9.2f %s %10s' % (
            r['label'], truth[r['label']]['kex'], r['kex'], r['chi2_global'],
            r['chi2'], r['delta_aicc'], zs, r['preferred']))
    print('-' * 78)
    rel = jk['se'] / abs(jk['kex_full']) * 100.0
    print('jackknife: kex_full=%.1f  SE=%.1f  (rel %.1f%%)  bias=%.1f' % (
        jk['kex_full'], jk['se'], rel, jk['bias']))
    flagged = [p['label'] for p in jk['per_residue'] if abs(p['delta_kex']) > 2 * jk['se']]
    print('jackknife influential (|delta_kex| > 2 SE): %s' % (flagged or 'none'))

    exp = EXPECT[name]
    print('-' * 78)
    if exp == 'split':
        # A single high-leverage outlier corrupts the shared kex for everyone
        # (least-squares is not robust), so the binary verdict flips to
        # individual for ALL residues.  Success = the tool *identifies* the
        # outlier (dominates |z| and chi2_global) and removing it restores a
        # clean global fit.
        pr = {r['label']: r for r in mc['per_residue']}
        z_out = abs(pr[MIXED_OUTLIER]['z_kex'])
        z_rest = max(abs(pr[l]['z_kex']) for l in pr if l != MIXED_OUTLIER)
        c_out = pr[MIXED_OUTLIER]['chi2_global']
        c_rest = max(pr[l]['chi2_global'] for l in pr if l != MIXED_OUTLIER)
        id_ok = z_out > 2 * z_rest and c_out > 3 * c_rest
        print('outlier %s identified: |z|=%.1f (rest max %.1f, %.1fx)  '
              'chi2_glob=%.0f (rest max %.0f, %.1fx)' % (
                  MIXED_OUTLIER, z_out, z_rest, z_out / z_rest,
                  c_out, c_rest, c_out / c_rest))

        rblocks = parse_blocks(os.path.join(d, 'out-recover.txt'))
        rec, rjk = rblocks['model_comparison'], rblocks['jackknife']
        rel_r = rjk['se'] / abs(rjk['kex_full']) * 100.0
        rec_ok = rec['preferred'] == 'global' and rec['n_pref_individual'] == 0
        print('recovery (%s removed): whole=%s  global-preferring=%d/%d  '
              'jackknife SE=%.1f%%  kex=%.0f' % (
                  MIXED_OUTLIER, rec['preferred'], rec['n_pref_global'],
                  rec['n_pref_global'] + rec['n_pref_individual'], rel_r, rjk['kex_full']))
        ok = id_ok and rec_ok
        print('EXPECT outlier flagged + global restored after removal -> %s' % (
            'PASS' if ok else 'FAIL'))
    else:
        ok = mc['preferred'] == exp
        print('EXPECT whole-dataset preferred = %s -> %s' % (
            exp, 'PASS' if ok else 'FAIL'))
    print()
    return name, exp, ok


if __name__ == '__main__':
    results = [evaluate(n) for n in ('global', 'individual', 'mixed')]
    print('=' * 78)
    print('SUMMARY')
    allok = True
    for name, exp, ok in results:
        allok = allok and ok
        print('  %-11s expected=%-11s %s' % (name, exp, 'PASS' if ok else 'FAIL'))
    raise SystemExit(0 if allok else 1)
