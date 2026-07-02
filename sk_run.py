#!/usr/bin/env python3
"""
sk_run.py — ShereKhan fitting runner.

Usage
-----
    python sk_run.py <config.json>

Reads a JSON configuration file that specifies the input .dat files, the
exchange regime, the fitting model, initial parameter values (or 'guess'),
and which residues to include.  Runs the CPMG_model fit, writes a .log file
and a PDF plot, and prints the detailed JSON results to stdout.

Config file keys
----------------
  Project Name : str   — base name for .log and .pdf output files
  datasets     : list  — paths to input .dat files (one per field)
  exchange     : str   — 'fast' or 'slow'
  model        : str   — 'Matrix', 'London', or 'Meiboom'
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

from sys import argv, stdout

import json

from cpmg.model_2state import CPMG_model


def main():
    """Run a ShereKhan fit from a JSON config file passed as argv[1]."""
    if len(argv) < 2:
        print('Usage: sk_run.py <config.json>')
        raise SystemExit(1)

    # -----------------------------------------------------------------------
    # Initialise the fitting model and print the program banner
    # -----------------------------------------------------------------------
    m2 = CPMG_model()
    m2.verbose = True  # print chi² at each function evaluation

    print('********')
    stdout.write(m2.programName)
    print('********')

    # -----------------------------------------------------------------------
    # Load the JSON configuration file
    # -----------------------------------------------------------------------
    with open(argv[1]) as configFile:
        conf = json.load(configFile)

    projectName   = conf['Project Name']  # base name for output files
    datasetsNames = conf['datasets']      # list of .dat file paths

    # Add each field's dispersion data to the dataset
    for dataset in datasetsNames:
        m2.dataset.addField(dataset)

    exchange = conf['exchange']  # 'fast' or 'slow'
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
                if active == 'on':
                    rsd.active = True
                elif active == 'off':
                    rsd.active = False
                else:
                    print('Error: wrong flag for the residue' + resId + ' ' + active)
                    raise SystemExit(1)

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
        elif m2.model == 'Matrix' or m2.model == 'London':
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
    if conf.get('compare_aic', False):
        print('Comparing global vs individual fits (AIC)...')
        aicCmp        = m2.compareModelsAIC()
        aicReportText = m2.aicReport(aicCmp)
        print(aicReportText)

    # If fast exchange was fitted with Matrix or London, convert results to
    # the Meiboom phi representation for a more interpretable parameter set
    if m2.dataset.exchange == 'fast' and m2.model != 'Meiboom':
        m2.toMeiboom()

    # -----------------------------------------------------------------------
    # Write outputs: log file, PDF, and JSON to stdout
    # -----------------------------------------------------------------------
    logBuf = m2.getLogBuffer()

    # Append the AIC model-comparison block to the log when it was computed
    if aicReportText is not None:
        logBuf += aicReportText

    with open(projectName + '.log', 'w') as file1:
        print(logBuf)
        file1.write(logBuf)

    m2.pdf(projectName + '.pdf')  # multi-page PDF with one residue per page

    out = m2.reportAllValues()    # nested list of per-residue result dicts

    print('Calculations finished successfully.')

    # JSON results are printed to stdout so they can be piped or captured
    print("#####")
    print(json.dumps(out))


if __name__ == '__main__':
    main()
