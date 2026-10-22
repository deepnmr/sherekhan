#!/usr/bin/env python3
"""
sk_run.py — ShereKhan fitting runner.

Usage
-----
    sk-run [--quiet] <config.json>

Reads a JSON configuration file that specifies the input .dat files, the
exchange regime, the fitting model, initial parameter values (or 'guess'),
and which residues to include.  Runs the CPMG_model fit, writes a .log file
and a PDF plot, and prints the detailed JSON results to stdout.

Config file keys
----------------
  Project Name : str   — base name for .log and .pdf output files
  datasets     : list  — paths to input .dat files (one per field); relative
                         paths resolve from the CWD, then from the config's dir
  exchange     : str   — 'fast', 'slow' or 'undefined' (runs as Matrix/London)
  model        : str   — 'Matrix', 'London', or 'Meiboom' (Meiboom needs 'fast')
  init
    mode       : str   — 'guess' (grid search) or 'values' (explicit)
    kex        : float — total exchange rate (s^-1), used when mode='values'
    pB         : float — minor-state population, used when mode='values'
                         and model != 'Meiboom'
    csd        : float — chemical shift difference (ppm), mode='values',
                         model != 'Meiboom'
    phi        : float — effective phi parameter (ppm^2), mode='values',
                         model == 'Meiboom'
  residues     : list  — dicts with 'name' and 'flag' ('on'/'off') keys
  compare_aic  : bool  — optional, global vs individual AIC comparison
  jackknife    : bool  — optional, leave-one-residue-out validation

Exit status: 0 on success, 1 on an input or runtime error (one
`sk-run: error: ...` line on stderr), 2 on a usage error.

Authors
-------
Prof. Dr. Donghan Lee
Korea Basic Science Institute (KBSI)

Original code:
Adam Mazur, Bjoern Hammesfahr, Christian Griesinger,
Donghan Lee, Martin Kollmar
Max-Planck-Institute for Biophysical Chemistry, Goettingen, 2012

Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
"""

###
# Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
# Original code: Adam Mazur et al., MPI Goettingen, 2012
###

from sys import stdout, stderr

import argparse
import json
import math
import os

from cpmg import __version__
from cpmg.model_2state import CPMG_model

EXCHANGES = ('fast', 'slow', 'undefined')
MODELS    = ('Matrix', 'London', 'Meiboom')
INIT_MODES = ('guess', 'values')
FLAGS     = ('on', 'off')


def _finite_json(obj):
    """Serialise to strict JSON, mapping inf/nan to null.

    The analysis blocks can hold non-finite floats (AICc is +inf when a residue
    is data-starved; z_kex is nan when the individual kex has no error estimate).
    Python's json emits Infinity/NaN tokens that strict parsers (e.g. JS
    JSON.parse) reject, so replace them with null for portability.
    """
    def clean(x):
        if isinstance(x, float):
            return x if math.isfinite(x) else None
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items()}
        if isinstance(x, list):
            return [clean(v) for v in x]
        return x
    return json.dumps(clean(obj))


def _is_number(x):
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return False
    try:
        return math.isfinite(float(x))
    except OverflowError:
        return False


def _check_config(conf, path):
    """Raise ValueError naming the offending key if the config is not runnable."""
    def fail(msg):
        raise ValueError('config %s: %s' % (path, msg))

    if not isinstance(conf, dict):
        fail('top level must be a JSON object')

    name = conf.get('Project Name')
    if not isinstance(name, str) or not name.strip():
        fail('"Project Name" must be a non-empty string')

    datasets = conf.get('datasets')
    if (not isinstance(datasets, list) or not datasets
            or not all(isinstance(d, str) for d in datasets)):
        fail('"datasets" must be a non-empty list of file paths')

    exchange = conf.get('exchange')
    if exchange not in EXCHANGES:
        fail('"exchange" must be one of %s (got %s)'
             % (', '.join(EXCHANGES), json.dumps(exchange)))

    model = conf.get('model')
    if model not in MODELS:
        fail('"model" must be one of %s (got %s)'
             % (', '.join(MODELS), json.dumps(model)))
    if model == 'Meiboom' and exchange != 'fast':
        fail('"model" Meiboom requires "exchange" to be "fast" (got %s)'
             % json.dumps(exchange))

    init = conf.get('init')
    if not isinstance(init, dict):
        fail('"init" must be an object')
    mode = init.get('mode')
    if mode not in INIT_MODES:
        fail('"init.mode" must be one of %s (got %s)'
             % (', '.join(INIT_MODES), json.dumps(mode)))
    if mode == 'values':
        if not _is_number(init.get('kex')) or init['kex'] <= 0:
            fail('"init.kex" must be a number > 0')
        if model == 'Meiboom':
            if not _is_number(init.get('phi')):
                fail('"init.phi" must be a number')
        else:
            pB = init.get('pB')
            if not _is_number(pB) or not 0 < pB < 1:
                fail('"init.pB" must be a number between 0 and 1')
            if not _is_number(init.get('csd')):
                fail('"init.csd" must be a number')

    residues = conf.get('residues')
    if not isinstance(residues, list):
        fail('"residues" must be a list')
    for i, rsd in enumerate(residues):
        if not isinstance(rsd, dict) or not isinstance(rsd.get('name'), str):
            fail('"residues[%d]" must be an object with a string "name"' % i)
        if rsd.get('flag') not in FLAGS:
            fail('"residues[%d].flag" must be "on" or "off" (got %s)'
                 % (i, json.dumps(rsd.get('flag'))))

    for key in ('compare_aic', 'jackknife'):
        if key in conf and not isinstance(conf[key], bool):
            fail('"%s" must be true or false' % key)


def _resolve_dataset(dataset, conf_dir):
    """Return the first existing path: as given (CWD), then relative to the config."""
    for candidate in (dataset, os.path.join(conf_dir, dataset)):
        if os.path.exists(candidate):
            return candidate
    raise ValueError('dataset not found: %s' % dataset)


def run(config_path, quiet):
    """Run the fit described by config_path. Raises on any input error."""
    with open(config_path) as configFile:
        try:
            conf = json.load(configFile)
        except json.JSONDecodeError as exc:
            raise ValueError('config %s: invalid JSON: %s' % (config_path, exc))
    _check_config(conf, config_path)
    conf_dir = os.path.dirname(os.path.abspath(config_path))

    logDir = os.path.dirname(conf['Project Name'] + '.log') or '.'
    if not os.path.isdir(logDir):
        raise ValueError('output directory does not exist: %s' % logDir)

    # -----------------------------------------------------------------------
    # Initialise the fitting model and print the program banner
    # -----------------------------------------------------------------------
    m2 = CPMG_model()
    m2.verbose = not quiet  # print chi² at each function evaluation

    print('********')
    stdout.write(m2.programName)
    print('********')

    projectName   = conf['Project Name']  # base name for output files
    datasetsNames = conf['datasets']      # list of .dat file paths

    # Add each field's dispersion data to the dataset
    for dataset in datasetsNames:
        m2.dataset.addField(_resolve_dataset(dataset, conf_dir))

    exchange = conf['exchange']  # 'fast', 'slow' or 'undefined'
    model    = conf['model']     # 'Matrix', 'London', or 'Meiboom'

    print('Exchange is: ' + exchange)
    print('Model is: '    + model)

    # -----------------------------------------------------------------------
    # Parse initial parameter settings
    # -----------------------------------------------------------------------
    initMode = conf['init']['mode']  # 'guess' or 'values'

    kex_0 = pB_0 = dd_0 = None  # initialised only when mode == 'values'

    if initMode == 'values':
        kex_0 = conf['init']['kex']
        if model != 'Meiboom':
            pB_0  = conf['init']['pB']
            dd_0  = conf['init']['csd']   # chemical shift difference (ppm)
        else:
            dd_0  = conf['init']['phi']   # phi = pB*(1-pB)*csd^2 (ppm^2)

    if initMode == 'guess':
        print('Initial parameters will be guessed')
    else:
        if model != 'Meiboom':
            print('Initial parameters are:\nkex: %8.3f\npB: %8.3f\ndd: %8.3f'
                  % (kex_0, pB_0, dd_0))
        else:
            print('Initial parameters are:\nkex: %8.3f\ndd: %8.3f' % (kex_0, dd_0))

    # -----------------------------------------------------------------------
    # Apply per-residue active/inactive flags from the config
    # -----------------------------------------------------------------------
    residues = conf['residues']  # list of {name, flag} dicts

    # Build a set of all residue labels present in the loaded dataset for fast lookup
    dataset_labels = {rsd.label for rsd in m2.dataset.rsds}

    for residue_conf in residues:
        resId  = residue_conf['name']
        active = residue_conf['flag']  # 'on' or 'off'
        print('Residue: ' + resId + ' ' + active)

        # Skip residues listed in the config but absent from the data (with a warning)
        if resId not in dataset_labels:
            print('Warning: residue ' + resId + ' not found in dataset, skipping.')
            continue

        for rsd in m2.dataset.rsds:
            if rsd.label == resId:
                rsd.active = (active == 'on')

    # -----------------------------------------------------------------------
    # Set the exchange regime and fitting model
    # -----------------------------------------------------------------------
    m2.dataset.exchange = exchange
    m2.model            = model

    # -----------------------------------------------------------------------
    # Initialise parameters (grid search or explicit values)
    # -----------------------------------------------------------------------
    if initMode == 'guess':
        # Grid-search over kex / pB / csd to find the best starting point
        m2.initGuessAll()
    elif initMode == 'values':
        m2.init_dd(dd_0)  # set chemical shift difference (or phi) for all residues
        if m2.model == 'Meiboom':
            m2.kex = kex_0
        elif m2.model in ('Matrix', 'London'):
            # Convert kex + pB into individual rate constants
            kAB    = pB_0 * kex_0
            kBA    = (1.0 - pB_0) * kex_0
            m2.kAB = kAB
            m2.kBA = kBA

    # Initialise R2_0 for every residue from the experimental data
    m2.init_R2_0()

    # -----------------------------------------------------------------------
    # Run the fitting
    # -----------------------------------------------------------------------
    m2.fit()

    # -----------------------------------------------------------------------
    # Optional: AIC comparison of the global fit against per-residue individual fits
    # -----------------------------------------------------------------------
    # Enabled with "compare_aic": true in the config.  Must run before toMeiboom()
    # so both models are compared in the same (csd) parameter domain.
    aicReportText = None
    aicCmp        = None
    if conf.get('compare_aic', False):
        print('Comparing global vs individual fits (AIC)...')
        aicCmp        = m2.compareModelsAIC()
        aicReportText = m2.aicReport(aicCmp)
        print(aicReportText)

    # -----------------------------------------------------------------------
    # Optional: jackknife validation of the global fit
    # -----------------------------------------------------------------------
    # Enabled with "jackknife": true in the config.  Leave-one-residue-out
    # refits check how strongly the shared global kex depends on any single
    # residue.  Like compare_aic, it must run before toMeiboom().
    jkReportText = None
    jk           = None
    if conf.get('jackknife', False):
        print('Jackknife validation of the global fit (leave-one-residue-out)...')
        jk = m2.jackknifeGlobal()
        if jk is not None:
            jkReportText = m2.jackknifeReport(jk)
            print(jkReportText)
        else:
            print('Jackknife skipped: need at least 2 active residues.')

    # If fast exchange was fitted with Matrix or London, convert results to
    # the Meiboom phi representation for a more interpretable parameter set
    if m2.dataset.exchange == 'fast' and m2.model != 'Meiboom':
        m2.toMeiboom()

    # -----------------------------------------------------------------------
    # Write outputs: log file, PDF, and JSON to stdout
    # -----------------------------------------------------------------------
    logBuf = m2.getLogBuffer(config_path)

    # Append the AIC model-comparison block to the log when it was computed
    if aicReportText is not None:
        logBuf += aicReportText

    # Append the jackknife validation block to the log when it was computed
    if jkReportText is not None:
        logBuf += jkReportText

    with open(projectName + '.log', 'w') as file1:
        print(logBuf)
        file1.write(logBuf)

    m2.pdf(projectName + '.pdf')  # multi-page PDF with one residue per page

    out = m2.reportAllValues()    # nested list of per-residue result dicts

    print('Calculations finished successfully.')

    # Optional analysis blocks: labelled JSON, printed before the residues so
    # the final line stays the residues list for backward-compatible parsing.
    if aicCmp is not None:
        print("##### model_comparison")
        print(_finite_json(aicCmp))
    if jk is not None:
        print("##### jackknife")
        print(_finite_json(jk))

    # JSON results are printed to stdout so they can be piped or captured
    print("#####")
    print(_finite_json(out))


class _UsageParser(argparse.ArgumentParser):
    """ArgumentParser that reports usage errors as one stderr line."""

    def error(self, message):
        self.exit(2, '%s: error: %s (use -h for help)\n' % (self.prog, message))


def main(argv=None):
    """Command-line entry point: parse options, run, and map input errors to exit 1."""
    parser = _UsageParser(
        prog='sk-run',
        description='Run a ShereKhan CPMG relaxation-dispersion fit from a JSON config.')
    parser.add_argument('config', help='JSON run configuration file (e.g. run.conf)')
    parser.add_argument('--quiet', action='store_true',
                        help='do not print chi2 progress lines during fitting')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)

    try:
        run(args.config, args.quiet)
    except KeyError as exc:
        stderr.write('sk-run: error: missing key %s\n' % exc)
        raise SystemExit(1)
    except (OSError, ValueError, TypeError) as exc:
        stderr.write('sk-run: error: %s\n' % exc)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
