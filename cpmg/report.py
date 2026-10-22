"""
report.py - Log, PDF and JSON reporting for the two-state CPMG fit.

Functions take the CPMG_model instance as first argument `m` and only read
its state. matplotlib is imported inside write_pdf(), so importing this module
does not load matplotlib.
"""

from getpass import getuser
from time import ctime
import platform
import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, exp, zeros, sort, tanh, isfinite
from . import models



def calc_r2(m, r, dspCurve, R2_0, v_array):
    """Compute R2 calculated values for a residue at given CPMG frequencies.

    Dispatches to the appropriate model (Matrix / London / Meiboom) and
    returns an array of R2eff values at all requested frequencies.

    For Matrix and London models, uses r.dd_copy (original csd) when
    available, falling back to r.dd.  r.dd_copy is set by toMeiboom() for
    active residues in fast-exchange runs; inactive residues retain r.dd
    as the original csd.

    Args:
        r        (Residue)         : Residue object with dd, R2_0, dd_copy.
        dspCurve (DispersionCurve) : Curve supplying field and tcp values.
        R2_0     (float)           : Intrinsic R2 for this field.
        v_array  (np.ndarray)      : CPMG frequencies at which to evaluate.

    Returns:
        np.ndarray: Calculated R2eff values (same length as v_array).
    """
    ty = zeros(v_array.shape)

    if m.model == 'Matrix':
        # Use the original csd (dd_copy) when toMeiboom() has already
        # replaced r.dd with phi; otherwise fall back to r.dd
        dd = r.dd_copy if r.dd_copy is not None else r.dd
        domega = dd * 2 * pi * dspCurve.field  # rad/s
        for ii in range(len(v_array)):
            ty[ii] = models.r2_matrix(m.kAB, m.kBA, domega, R2_0, v_array[ii], dspCurve.tcp)

    elif m.model == 'London':
        # Same dd_copy guard as Matrix model above
        dd = r.dd_copy if r.dd_copy is not None else r.dd
        domega = dd * 2 * pi * dspCurve.field  # rad/s
        for ii in range(len(v_array)):
            ty[ii] = models.r2_london(m.kAB, m.kBA, domega, R2_0, v_array[ii])

    elif m.model == 'Meiboom':
        # r.dd holds phi = pB*(1-pB)*dd^2 for the Meiboom model
        for ii in range(len(v_array)):
            ty[ii] = models.r2_meiboom(m.kex, r.dd, R2_0, v_array[ii], dspCurve.field)

    return ty


def log_buffer(m, config_path=None):
    """Build a formatted text report of the fitting results.

    Includes the program header, run metadata (user, hostname, timestamp),
    the contents of the input config file, all datasets with experimental
    and calculated R2 values, fit statistics, and per-residue parameters
    with standard deviations.

    Returns:
        str: Multi-line text suitable for writing directly to a .log file.
    """
    buf = '********\n'
    buf += m.programName
    buf += '********\n\n'

    # Run metadata: user, hostname, and current time
    buf += '********\n'
    user     = getuser()  # getlogin() needs a tty; getuser() reads env, works headless
    hostname = platform.node()
    buf += 'User: %s@%s\n' % (user, hostname)
    buf += '%s\n' % ctime()
    buf += '********\n\n'

    # Echo the input config file for full reproducibility
    buf += '********\n\n'
    if config_path is None:
        buf += 'Input file: (not provided)\n'
    else:
        buf += 'Input file:\n'
        with open(config_path, 'r') as file1:
            buf += file1.read()
    buf += '********\n\n'

    # Dataset summary: for each field, list experimental and calculated R2
    buf += 'Dataset(s):\n'
    for fi in range(0, len(m.dataset.fields)):
        buf += 'Field, tcp [MHz, s]:\n'
        buf += '%5.2f\n' % m.dataset.fields[fi]
        buf += '%f\n'    % m.dataset.tcps[fi]
        buf += '#v_cpmg(Hz)   R2(1/s)      Esd(R2)      R2_calc\n'

        for r in m.dataset.rsds:
            if r.active:
                buf += '# %s\n' % r.label
                for i in range(0, len(r.dspCurves)):
                    dspCurve = r.dspCurves[i]
                    R2_0 = r.R2_0[i]

                    va     = array(dspCurve.v)
                    R2_calc = calc_r2(m, r, dspCurve, R2_0, va)

                    # Only print the row for the field being reported
                    if dspCurve.field == m.dataset.fields[fi]:
                        for j in range(0, len(dspCurve.v)):
                            buf += '%8.3f %12.3f %12.3f %12.3f\n' % (
                                dspCurve.v[j], dspCurve.R2exp[j],
                                dspCurve.R2stddev[j], R2_calc[j])

    # Fit statistics
    buf += '\n\n\n---\n'
    buf += 'Exchange regime: %s\n' % m.dataset.exchange
    buf += 'Model: %s\n' % m.model
    buf += '---\n'
    # Guard against dof ≤ 0 to prevent division by zero
    red_chi2 = m.chi2 / m.dof if m.dof > 0 else float('nan')
    buf += 'npar=%d nvar=%d ndof=%d chi2=%8.3f chi2/dof=%8.3f\n' % (
        m.npar, m.nvar, m.dof, m.chi2, red_chi2)
    buf += '---\n'

    # Per-exchange-regime parameter summary
    if m.dataset.exchange in ('slow', 'undefined'):
        # Derive kex and pB from the individual rate constants
        kex, kex_std, pB, pB_std = m._kex_pB()

        buf += 'kAB: %8.3f  +-  %8.3f\n' % (m.kAB, m.kAB_std)
        buf += 'kBA: %8.3f  +-  %8.3f\n' % (m.kBA, m.kBA_std)
        buf += '---\n'
        buf += 'kex: %8.3f  +-  %8.3f\n' % (kex, kex_std)
        buf += ' pB: %8.3f  +-  %8.3f\n' % (pB, pB_std)
        buf += '---\n'
        buf += 'resId:         csd [ppm]   '
        for field in m.dataset.fields:
            buf += '       R2_0 (%5.1f MHz) ' % field
        buf += '\n---\n'

    if m.dataset.exchange == 'fast':
        buf += 'kex: %8.3f  +-  %8.3f\n' % (m.kex, m.kex_std)
        buf += '---\n'
        buf += 'resId:         phi [ppm^2] '
        for field in m.dataset.fields:
            buf += '       R2_0 (%5.1f MHz) ' % field
        buf += '\n---\n'

    # Per-residue parameters (dd/phi and R2_0 at each field)
    for r in m.dataset.rsds:
        if r.active:
            buf += '%5s: %8.3f +- %8.3f' % (r.label, r.dd, r.dd_std)
            for cf in m.dataset.fields:
                flag = False
                for i in range(0, len(r.dspCurves)):
                    if r.dspCurves[i].field == cf:
                        buf += '    %8.3f' % r.R2_0[i]
                        buf += ' +- %8.3f' % r.R2_0_std[i]
                        flag = True
                if not flag:
                    buf += '    --------'
                    buf += ' +- --------'
            buf += '\n'
        else:
            buf += '%5s:    excluded from the fit\n' % r.label
    buf += '---\nEnd.\n'

    return buf


def report_all_values(m):
    """Export all per-residue fit results as a nested list of dictionaries.

    For each residue and each dispersion curve (one per field), builds a
    dictionary containing the calculated R2 curve, exchange parameters,
    and original experimental data.  The structure is::

        [
          [  # residue 0
            { field 0 data },
            { field 1 data },
            ...
          ],
          ...  # residue 1, 2, ...
        ]

    Returns:
        list: Nested list suitable for JSON serialisation via json.dumps().
    """
    all_data = []

    for r in m.dataset.rsds:
        all_data.append([])
        for i, dspCurve in enumerate(r.dspCurves):
            tx   = array(dspCurve.v)
            R2_0 = r.R2_0[i]
            ty   = calc_r2(m, r, dspCurve, R2_0, tx)  # calculated R2

            if m.dataset.exchange in ('slow', 'undefined'):
                kAB        = m.kAB
                kBA        = m.kBA
                kex        = kAB + kBA
                pB         = kAB / kex
                chemShift  = r.dd       # chemical shift difference (ppm)
                chemShift_std = r.dd_std

            elif m.dataset.exchange == 'fast':
                kex     = m.kex
                phi     = r.dd          # effective Meiboom phi (ppm^2)
                phi_std = r.dd_std

            # Append the per-field result dictionary
            if m.dataset.exchange in ('slow', 'undefined'):
                all_data[-1].append({
                    'tx':                 list(tx),
                    'r2_clc':             list(ty),
                    'kAB':                kAB,
                    'kBA':                kBA,
                    'kex':                kex,
                    'pB':                 pB,
                    'chemShiftDiff':      chemShift,
                    'chemShiftDiffStdDev': chemShift_std,
                    'r20':                R2_0,
                    'x1':                 dspCurve.v,
                    'y_val':              dspCurve.R2exp,
                    'y_val_err':          dspCurve.R2stddev,
                    'field':              dspCurve.field,
                    'label':              r.label
                })
            elif m.dataset.exchange == 'fast':
                all_data[-1].append({
                    'tx':        list(tx),
                    'r2_clc':    list(ty),
                    'kex':       kex,
                    'phi':       phi,
                    'phiStdDev': phi_std,
                    'r20':       R2_0,
                    'x1':        dspCurve.v,
                    'y_val':     dspCurve.R2exp,
                    'y_val_err': dspCurve.R2stddev,
                    'field':     dspCurve.field,
                    'label':     r.label
                })

    return all_data


def aic_report(model_name, cmp):
    """Format the AIC model-comparison result as a text block.

    Args:
        model_name (str): 'Matrix', 'London' or 'Meiboom' (selects the column layout).
        cmp (dict): the dict returned by compareModelsAIC().

    Returns:
        str: multi-line report suitable for the log file / stdout.
    """
    g = cmp['global']
    i = cmp['individual']
    buf  = '\n========\n'
    buf += 'Model comparison: global vs individual fit (AIC)\n'
    buf += '========\n'
    buf += 'Model         %10s %6s %6s %12s %12s\n' % (
        'chi2', 'k', 'n', 'AIC', 'AICc')
    buf += 'global      %12.3f %6d %6d %12.3f %12.3f\n' % (
        g['chi2'], g['k'], g['n'], g['aic'], g['aicc'])
    buf += 'individual  %12.3f %6d %6d %12.3f %12.3f\n' % (
        i['chi2'], i['k'], i['n'], i['aic'], i['aicc'])
    buf += '---\n'
    buf += 'delta AIC  (individual - global): %12.3f\n' % cmp['delta_aic']
    buf += 'delta AICc (individual - global): %12.3f\n' % cmp['delta_aicc']
    buf += 'Akaike weights (AICc): global=%.4f  individual=%.4f\n' % (
        cmp['weight_global'], cmp['weight_individual'])
    buf += 'Preferred model (whole dataset): %s\n' % cmp['preferred']
    buf += '---\n'
    buf += 'Per-residue individual fit:\n'
    if model_name == 'Meiboom':
        buf += '%6s %10s %12s %10s\n' % ('resId', 'kex', 'phi', 'chi2')
        for d in cmp['per_residue']:
            buf += '%6s %10.3f %12.5f %10.3f\n' % (
                d['label'], d['kex'], d['phi'], d['chi2'])
    else:
        buf += '%6s %10s %8s %10s %10s\n' % (
            'resId', 'kex', 'pB', 'csd', 'chi2')
        for d in cmp['per_residue']:
            buf += '%6s %10.3f %8.4f %10.4f %10.3f\n' % (
                d['label'], d['kex'], d['pB'], d['dd'], d['chi2'])

    # Per-residue verdict: individual vs global for each residue
    buf += '---\n'
    buf += 'Per-residue preference (individual vs global fit):\n'
    buf += '%6s %12s %12s %12s %12s %10s\n' % (
        'resId', 'chi2_glob', 'chi2_indiv', 'dAICc', 'z(kex)', 'better')
    for d in cmp['per_residue']:
        z = d.get('z_kex', float('nan'))
        buf += '%6s %12.3f %12.3f %12.3f %12.2f %10s\n' % (
            d['label'], d['chi2_global'], d['chi2'],
            d['delta_aicc'], z, d['preferred'])
    buf += '---\n'
    buf += ('dAICc = AICc(individual) - AICc(global): '
            '> 0 favours global, < 0 favours individual.\n')
    buf += ('z(kex) = (kex_individual - kex_leave_one_out) / sigma(kex_individual), '
            'i.e. vs the kex fitted from the other residues (unbiased): '
            '|z| > 2 means the residue disagrees with the shared kex.\n')
    buf += 'Residues preferring global: %d   preferring individual: %d\n' % (
        cmp['n_pref_global'], cmp['n_pref_individual'])
    buf += '========\n'
    return buf


def jackknife_report(jk):
    """Format the jackknife validation result as a text block.

    Args:
        jk (dict): the dict returned by jackknifeGlobal().

    Returns:
        str: multi-line report suitable for the log file / stdout.
    """
    buf  = '\n========\n'
    buf += 'Jackknife validation of the global fit (leave-one-residue-out)\n'
    buf += '========\n'
    buf += 'active residues resampled:    %12d\n'   % jk['n']
    buf += 'kex (full global fit):        %12.3f\n' % jk['kex_full']
    buf += 'kex (jackknife mean):         %12.3f\n' % jk['kex_jack_mean']
    buf += 'kex (bias-corrected):         %12.3f\n' % jk['kex_jack_est']
    buf += 'jackknife bias:               %12.3f\n' % jk['bias']
    buf += 'jackknife std error (kex):    %12.3f\n' % jk['se']
    rel = jk['se'] / abs(jk['kex_full']) * 100.0 if jk['kex_full'] else float('nan')
    buf += 'relative std error:           %11.2f%%\n' % rel
    buf += '---\n'
    buf += 'Per-residue influence (delta_kex = kex_full - kex_without_residue):\n'
    buf += '%6s %14s %14s %10s\n' % ('resId', 'kex(-res)', 'delta_kex', '|z|')
    for d in jk['per_residue']:
        z = abs(d['delta_kex']) / jk['se'] if jk['se'] > 0 else float('nan')
        flag = '  *' if (jk['se'] > 0 and abs(d['delta_kex']) > 2.0 * jk['se']) else ''
        buf += '%6s %14.3f %14.3f %10.2f%s\n' % (
            d['label'], d['kex_drop'], d['delta_kex'], z, flag)
    buf += '---\n'
    buf += '(* = influential: dropping this residue shifts kex by > 2 jackknife SE)\n'
    buf += 'The global fit is robust when no residue is flagged and the '
    buf += 'relative SE is small.\n'
    buf += '========\n'
    return buf


def write_pdf(m, pdfFileName):
    """Generate a multi-page PDF with one dispersion-curve plot per residue.

    Each page shows experimental R2 values (with error bars) and the
    corresponding calculated curve (solid line) for each magnetic field,
    colour-coded consistently.

    Args:
        pdfFileName (str): Path of the output PDF file to create.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(pdfFileName)

    for r in m.dataset.rsds:
        if r.active:
            fig = Figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            # Fixed colour set, indexed by curve position
            colorsSet = ['b', 'g', 'r', 'c', 'm', 'y']

            max_x = []
            min_y = []
            max_y = []

            for i in range(0, len(r.dspCurves)):
                dspCurve = r.dspCurves[i]
                max_x.append(max(dspCurve.v))

                # Assign colour by position (fall back to black when exhausted)
                c = colorsSet[i] if i < len(colorsSet) else 'k'

                # Plot experimental data as dots with error bars
                ax.errorbar(dspCurve.v, dspCurve.R2exp,
                            yerr=dspCurve.R2stddev,
                            fmt='%so' % c,
                            label=' exp %5.1f' % dspCurve.field)

                # Evaluate calculated curve on a sorted frequency array
                tx = sort(array(dspCurve.v))
                R2_0 = r.R2_0[i]
                ty = calc_r2(m, r, dspCurve, R2_0, tx)
                ax.plot(tx, ty, '%s-' % c, label='calc %5.1f' % dspCurve.field)

                # Collect y-range for axis scaling
                min_yv = [v - e for v, e in zip(dspCurve.R2exp, dspCurve.R2stddev)]
                max_yv = [v + e for v, e in zip(dspCurve.R2exp, dspCurve.R2stddev)]
                min_y.append(min(min_yv))
                max_y.append(max(max_yv))

            ax.set_title(r.label)
            ax.set_xlabel('v [Hz]')
            ax.set_ylabel('R2 [1/s]')

            # Shrink plot area by 20% on the right to make room for the legend
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            ax.grid(True)
            ax.set_ylim(0.95 * min(min_y), 1.05 * max(max_y))
            ax.set_xlim(0.0, 1.05 * max(max_x))
            pdf.savefig(fig)

    pdf.close()
