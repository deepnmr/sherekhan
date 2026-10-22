"""Golden regression tests: sk-run output vs values captured at the baseline run.

Each test/golden/<name>.json holds the full per-record output of one bundled
config (every key, lists included) plus the pinned model-comparison and
jackknife values for fast-globalfit. Comparison is strict on structure:
record key sets must match exactly, lists are compared element-wise, and
numbers use rel_tol=REL_TOL / abs_tol=ABS_TOL.
Run with pytest, or as a script: python test/test_golden.py
"""
import glob
import json
import math
import os
import re
import subprocess
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(TEST_DIR)
SK_RUN = os.path.join(ROOT_DIR, 'sk_run.py')
GOLDEN_DIR = os.path.join(TEST_DIR, 'golden')
REL_TOL = 1e-4
ABS_TOL = 1e-9
CHI2_LINE = re.compile(r'npar=(\d+) nvar=(\d+) ndof=(\d+) chi2=\s*([-\d.eE+]+)')


def _conf_dir(name):
    return os.path.join(TEST_DIR, 'fast' if name.startswith('fast') else 'slow')


def _run(name):
    conf = name + '.conf'
    env = dict(os.environ, MPLBACKEND='Agg')
    proc = subprocess.run([sys.executable, SK_RUN, conf], cwd=_conf_dir(name),
                          env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    assert proc.returncode == 0, 'sk_run %s failed:\n%s' % (conf, proc.stderr)
    return proc.stdout


def _blocks(stdout):
    """Map each '#####' marker to the JSON line that follows it."""
    lines = stdout.split('\n')
    markers = [i for i, line in enumerate(lines) if line.startswith('#####')]
    return [lines[i] for i in markers], [json.loads(lines[i + 1]) for i in markers]


def _chi2_fields(stdout):
    lines = [l for l in stdout.split('\n') if l.startswith('npar=')]
    assert lines, 'no npar=... chi2 line in stdout'
    m = CHI2_LINE.search(lines[-1])
    assert m, 'unparsable chi2 line: %r' % lines[-1]
    return int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4))


def _is_num(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _compare(actual, expected, where):
    """Recursive strict comparison: exact key sets, element-wise lists, tolerant numbers."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict), '%s: expected object, got %r' % (where, type(actual).__name__)
        assert set(actual) == set(expected), \
            '%s: key set changed: added %s, removed %s' % (
                where, sorted(set(actual) - set(expected)), sorted(set(expected) - set(actual)))
        for key in expected:
            _compare(actual[key], expected[key], '%s.%s' % (where, key))
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), \
            '%s: list length %s != golden %d' % (where, len(actual) if isinstance(actual, list) else type(actual).__name__, len(expected))
        for i, (got, want) in enumerate(zip(actual, expected)):
            _compare(got, want, '%s[%d]' % (where, i))
    elif _is_num(expected):
        assert _is_num(actual), '%s: expected number, got %r' % (where, actual)
        assert math.isclose(actual, expected, rel_tol=REL_TOL, abs_tol=ABS_TOL), \
            '%s: %r != golden %r' % (where, actual, expected)
    else:
        assert actual == expected, '%s: %r != golden %r' % (where, actual, expected)


def _project(actual, expected, where):
    """Keep only the keys pinned by the golden (for the partial model_comparison / jackknife blocks)."""
    if isinstance(expected, dict):
        missing = sorted(set(expected) - set(actual))
        assert not missing, '%s: missing keys %s' % (where, missing)
        return {k: _project(actual[k], expected[k], '%s.%s' % (where, k)) for k in expected}
    if isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), '%s: list length differs' % where
        return [_project(a, e, '%s[%d]' % (where, i)) for i, (a, e) in enumerate(zip(actual, expected))]
    return actual


def _check_stdout(name, stdout):
    with open(os.path.join(GOLDEN_DIR, name + '.json')) as f:
        golden = json.load(f)
    npar, nvar, ndof, chi2 = _chi2_fields(stdout)
    assert (npar, nvar, ndof) == (golden['npar'], golden['nvar'], golden['ndof']), \
        '%s: npar/nvar/ndof %s != golden %s' % (name, (npar, nvar, ndof),
                                              (golden['npar'], golden['nvar'], golden['ndof']))
    _compare(chi2, golden['chi2'], name + ' chi2')
    markers, blocks = _blocks(stdout)
    assert markers[-1] == '#####', 'final marker missing: %r' % markers
    actual = [r for res in blocks[-1] for r in res]
    records = golden['records']
    assert len(actual) == len(records), '%s: %d records, golden %d' % (name, len(actual), len(records))
    for i, (got, want) in enumerate(zip(actual, records)):
        _compare(got, want, '%s record %d (%s)' % (name, i, want.get('label')))
    if 'model_comparison' in golden:
        mc = blocks[markers.index('##### model_comparison')]
        _compare(_project(mc, golden['model_comparison'], name + ' model_comparison'),
                 golden['model_comparison'], name + ' model_comparison')
        jk = blocks[markers.index('##### jackknife')]
        for entry in jk['per_residue']:
            entry.pop('influential', None)  # the only key jackknife per_residue may add
        _compare(_project(jk, golden['jackknife'], name + ' jackknife'),
                 golden['jackknife'], name + ' jackknife')


def _check_golden(name):
    _check_stdout(name, _run(name))


def _names():
    return sorted(os.path.splitext(os.path.basename(p))[0]
                  for p in glob.glob(os.path.join(GOLDEN_DIR, '*.json')))


def _replace_final_record(stdout, mutate):
    lines = stdout.rstrip('\n').split('\n')
    record = json.loads(lines[-1])
    mutate(record)
    lines[-1] = json.dumps(record)
    return '\n'.join(lines) + '\n'


def _mutate_r2_clc_1pct(record):
    record[0][0]['r2_clc'][0] *= 1.01


def _mutate_renamed_key(record):
    record[0][0]['kex_renamed'] = record[0][0].pop('kex')


def test_golden_all():
    names = _names()
    assert len(names) == 5, names
    for name in names:
        _check_golden(name)


def test_compare_rejects_drift():
    """The compare must fail on a +1% r2_clc element and on a renamed record key."""
    stdout = _run('fast')
    for mutate, label in ((_mutate_r2_clc_1pct, 'r2_clc +1%'), (_mutate_renamed_key, 'renamed key')):
        try:
            _check_stdout('fast', _replace_final_record(stdout, mutate))
        except AssertionError as err:
            print('expected failure (%s): %s' % (label, err))
        else:
            raise AssertionError('compare accepted drift: %s' % label)


if __name__ == '__main__':
    for n in _names():
        _check_golden(n)
        print('golden ok: %s' % n)
    test_compare_rejects_drift()
    print('drift detection ok')
    print('test_golden: all passed')
