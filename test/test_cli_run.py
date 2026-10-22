#!/usr/bin/env python3
"""CLI contract tests for sk_run.py (help/version, error matrix, paths, --quiet).

Run directly: `python3 test/test_cli_run.py` (or `python -m pytest test/test_cli_run.py`).
Every case runs sk_run.py as a subprocess in a fresh temporary directory.
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SK_RUN = os.path.join(ROOT, 'sk_run.py')
DATA = os.path.join(ROOT, 'test', 'fast', 'in', 'fast-synth-60.dat')
SLOW_CONF = os.path.join(ROOT, 'test', 'slow', 'slow.conf')
ENV = dict(os.environ, MPLBACKEND='Agg')


def _run(args, cwd):
    return subprocess.run([sys.executable, SK_RUN] + args, cwd=cwd, env=ENV,
                          capture_output=True, text=True)


def _base_conf(**over):
    conf = {
        'Project Name': 'cli_case',
        'datasets': [DATA],
        'exchange': 'fast',
        'model': 'Matrix',
        'init': {'mode': 'guess'},
        'residues': [{'name': 'K1f', 'flag': 'on'}],
    }
    conf.update(over)
    return conf


def _write_conf(d, conf, name='run.conf'):
    path = os.path.join(d, name)
    with open(path, 'w') as f:
        if isinstance(conf, str):
            f.write(conf)
        else:
            json.dump(conf, f)
    return path


def test_help_and_version_write_nothing():
    with tempfile.TemporaryDirectory() as d:
        r = _run(['--help'], d)
        assert r.returncode == 0, r.stderr
        assert '--quiet' in r.stdout
        r = _run(['--version'], d)
        assert r.returncode == 0, r.stderr
        assert r.stdout.startswith('sk-run '), r.stdout
        assert os.listdir(d) == [], os.listdir(d)


def test_no_args_is_usage_error():
    with tempfile.TemporaryDirectory() as d:
        r = _run([], d)
        assert r.returncode == 2, r.returncode
        assert 'Traceback' not in r.stderr


def test_invalid_configs_give_one_line_error():
    cases = {
        'missing file': None,
        'malformed json': '{"Project Name": ',
        'missing dataset': _base_conf(datasets=[os.path.join(ROOT, 'nope.dat')]),
        'unknown model': _base_conf(model='Foo'),
        'exchange case': _base_conf(exchange='Fast'),
        'meiboom with slow': _base_conf(model='Meiboom', exchange='slow',
                                        init={'mode': 'guess'}),
        'bad init mode': _base_conf(init={'mode': 'foo'}),
        'bad flag': _base_conf(residues=[{'name': 'K1f', 'flag': 'maybe'}]),
        'pB out of range': _base_conf(init={'mode': 'values', 'kex': 100,
                                            'pB': 1.5, 'csd': 1.0}),
        'non-bool jackknife': _base_conf(jackknife='yes'),
    }
    with tempfile.TemporaryDirectory() as d:
        for label, conf in cases.items():
            if conf is None:
                path = os.path.join(d, 'nonexistent.conf')
            else:
                path = _write_conf(d, conf)
            r = _run([path], d)
            combined = r.stdout + r.stderr
            assert r.returncode == 1, (label, r.returncode, combined)
            assert r.stderr.startswith('sk-run: error:'), (label, r.stderr)
            assert 'Traceback' not in combined, (label, combined)


def test_overflowing_number_is_invalid_config():
    huge = int('1' + '0' * 400)
    conf = _base_conf(init={'mode': 'values', 'kex': huge, 'pB': 0.05, 'csd': 1.0})
    with tempfile.TemporaryDirectory() as d:
        path = _write_conf(d, conf)
        r = _run([path], d)
        assert r.returncode == 1, (r.returncode, r.stderr)
        assert r.stderr.count('\n') == 1, r.stderr
        assert '"init.kex"' in r.stderr, r.stderr
        assert 'Traceback' not in r.stdout + r.stderr


def test_missing_output_directory_fails_before_fitting():
    conf = _base_conf(**{'Project Name': 'no_such_dir/run'})
    with tempfile.TemporaryDirectory() as d:
        path = _write_conf(d, conf)
        r = _run([path], d)
        assert r.returncode == 1, (r.returncode, r.stderr)
        assert r.stderr.count('\n') == 1, r.stderr
        assert 'output directory does not exist' in r.stderr, r.stderr
        assert 'chi2' not in r.stdout and 'red_chi2=' not in r.stdout
        assert os.listdir(d) == ['run.conf'], os.listdir(d)


def test_usage_error_is_one_stderr_line():
    with tempfile.TemporaryDirectory() as d:
        r = _run(['--bogus'], d)
        assert r.returncode == 2, r.returncode
        lines = [ln for ln in r.stderr.splitlines() if ln]
        assert len(lines) == 1, r.stderr
        assert lines[0].startswith('sk-run: error: '), lines[0]
        assert lines[0].endswith('(use -h for help)'), lines[0]
        assert 'Traceback' not in r.stderr


def test_relative_datasets_resolve_from_config_dir():
    with tempfile.TemporaryDirectory() as cfg_dir, tempfile.TemporaryDirectory() as cwd:
        rel = os.path.relpath(DATA, os.path.realpath(cfg_dir))
        path = _write_conf(cfg_dir, _base_conf(datasets=[rel]))
        r = _run([path], cwd)
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(cwd, 'cli_case.log'))


def test_quiet_keeps_final_json_and_drops_chi2_lines():
    with tempfile.TemporaryDirectory() as d:
        loud = _run([SLOW_CONF], d)
        assert loud.returncode == 0, loud.stderr
        quiet = _run(['--quiet', SLOW_CONF], d)
        assert quiet.returncode == 0, quiet.stderr
        assert 'red_chi2=' in loud.stdout
        assert 'red_chi2=' not in quiet.stdout
        last_loud = loud.stdout.rstrip('\n').split('\n')[-1]
        last_quiet = quiet.stdout.rstrip('\n').split('\n')[-1]
        assert last_loud == last_quiet
        # Strict JSON: NaN/Infinity tokens would make this raise.
        json.loads(last_quiet, parse_constant=lambda tok: 1 / 0)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for t in tests:
        t()
        print('ok  ' + t.__name__)
    print('all %d cli tests passed' % len(tests))


if __name__ == '__main__':
    main()
