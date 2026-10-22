"""
models.py - Pure two-state CPMG R2 dispersion formulas.

Every function here is stateless: it takes rates, shifts and CPMG frequencies
and returns R2eff. The CPMG_model facade in model_2state.py keeps the public
method names and delegates to these functions.
"""

import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, exp, zeros, sort, tanh, isfinite
from numpy.linalg import matrix_power
from numpy import diag
from scipy.linalg import expm


def _echo_count(tcp, vCPMG):
    """Return the integer number of CPMG echo repeats n = tcp * vCPMG.

    Raises:
        ValueError: when tcp * vCPMG is not a positive integer within 1e-3.
    """
    fn = tcp * vCPMG
    n = int(round(fn))
    if abs(fn - n) > 1e-3 or n < 1:
        raise ValueError('tcp*nu_cpmg = %g is not a positive integer number of CPMG '
                         'cycles (tcp=%g s, nu=%g Hz)' % (fn, tcp, vCPMG))
    return n

def _london_cr72(kAB, kBA, domega, R2_0, nu):
    """Infinite-train (Carver-Richards 1972) R2eff for a London point whose closed form is non-finite."""
    if domega == 0 or kAB == 0:
        return R2_0
    d = 1.0 / (4.0 * nu)
    A = array([[-kAB, kBA], [kAB, -kBA + 1.0j * domega]])
    B = dot(dot(expm(A * d), expm(A.conj() * 2.0 * d)), expm(A * d))
    rho = max(abs(np.linalg.eigvals(B)))
    return R2_0 - nu * log(rho)

def r2_matrix(kAB, kBA, domega, R2_0, vCPMG, tcp):
    """Compute R2eff using the exact Liouvillian matrix-exponential method.

    Constructs the 2×2 Liouvillian superoperator for two-site chemical
    exchange and propagates the magnetisation vector through one complete
    CPMG echo train of n repeats, each consisting of four delta intervals
    (delta = 1 / (4 * nu_CPMG))::

        B = exp(A*delta) · exp(A*·delta) · exp(A*·delta) · exp(A*delta)
        C = B^n
        R2eff = -1/(4*n*delta) * ln(C[0,0] / M0[0])

    Args:
        kAB    (float): A→B exchange rate (s^-1).
        kBA    (float): B→A exchange rate (s^-1).
        domega (float): Chemical shift difference in rad/s.
        R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
        vCPMG  (float): CPMG pulsing frequency (Hz).
        tcp    (float): Total CPMG pulse delay (s); used to derive n.

    Returns:
        float: Calculated effective R2 relaxation rate (s^-1).
    """
    delta = 1.0 / (4.0 * vCPMG)   # time for one quarter-echo interval
    n = _echo_count(tcp, vCPMG)    # number of CPMG repeats (integer)
    if kAB == 0 or kAB + kBA == 0 or domega == 0:
        return R2_0                # no exchange or no chemical shift: no R2 contribution

    # Build the Liouvillian matrix A for state A (real) and state B (complex)
    Alist = [[-R2_0 - kAB, kBA],
             [kAB, -R2_0 - kBA + 1.0j * domega]]

    A = array(Alist)
    eA = expm(A * delta)           # propagator for forward interval

    Astar = A.conj()               # complex conjugate (180° pulse effect)
    eAstar = expm(Astar * delta)   # propagator for refocused interval

    # One full echo: delta → delta* → delta* → delta (CPMG refocusing)
    B = dot(dot(dot(eA, eAstar), eAstar), eA)
    # Repeat n times for the full CPMG train
    C = matrix_power(B, n)

    # Equilibrium populations
    pB = kAB / (kAB + kBA)
    M0 = array([1.0 - pB, pB])  # initial magnetisation vector

    MT = dot(C, M0)  # magnetisation after the full CPMG sequence

    # Recover R2eff from the decay of the ground-state magnetisation
    R2_clc = (-1.0 / (4.0 * n * delta)) * log(MT[0] / M0[0])

    return real(R2_clc)

def r2_matrix_eig(kAB, kBA, domega, R2_0, vCPMG, tcp):
    """Compute R2eff via eigenvalue decomposition of the Liouvillian.

    Alternative matrix formulation used inside errFunc during optimisation.
    The propagator is built from the analytic eigenvalues (l1, l2) and
    eigenvectors of A, which can be more numerically stable for certain
    parameter regimes.

    Args:
        kAB    (float): A→B exchange rate (s^-1).
        kBA    (float): B→A exchange rate (s^-1).
        domega (float): Chemical shift difference in rad/s.
        R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
        vCPMG  (float): CPMG pulsing frequency (Hz).
        tcp    (float): Total CPMG pulse delay (s); used to derive n.

    Returns:
        float: Calculated effective R2 relaxation rate (s^-1).
    """
    delta = 1.0 / (4.0 * vCPMG)
    n = _echo_count(tcp, vCPMG)

    if kAB == 0 or kAB + kBA == 0 or domega == 0:
        return R2_0

    # Analytic square-root factor appearing in both eigenvalues
    AA = sqrt(kBA ** 2 + (2.0 * kAB - 2.0 * 1.0j * domega) * kBA
              + kAB ** 2 + 2.0 * 1.0j * domega * kAB - domega ** 2)

    # Eigenvalues of A * delta
    l1 = -(delta * AA + delta * kBA + delta * kAB + (2.0 * R2_0 - 1.0j * domega) * delta) / 2.0
    l2 =  (delta * AA - delta * kBA - delta * kAB + (1.0j * domega - 2.0 * R2_0) * delta) / 2.0

    L = array([l1, l2])

    # Eigenvectors of A and their inverse
    v1 = (2 * kAB) / (-AA + kBA - kAB - 1.0j * domega)
    v2 = (2 * kAB) / ( AA + kBA - kAB - 1.0j * domega)
    X  = array([[1.0, 1.0], [v1, v2]])    # columns are eigenvectors
    Xi = array([
            [ v2 / (v2 - v1), -1.0 / (v2 - v1)],
            [-v1 / (v2 - v1),  1.0 / (v2 - v1)]
           ])                              # inverse of X

    # Complex conjugates for the refocused (Astar) propagator
    Lstar  = L.conj()
    Xstar  = X.conj()
    Xistar = Xi.conj()

    # Forward propagator:  exp(A*delta) = X · diag(exp(L)) · X^-1
    eA = dot(dot(X, diag(exp(L))), Xi)

    # Refocused propagator squared:  exp(Astar*2*delta)
    eAstar2 = dot(dot(Xstar, diag(exp(Lstar) ** 2)), Xistar)

    # One full echo and n repeats
    B = dot(dot(eA, eAstar2), eA)
    C = matrix_power(B, n)

    pB = kAB / (kAB + kBA)
    M0 = array([1.0 - pB, pB])
    MT = dot(C, M0)

    R2_clc = (-1.0 / (4.0 * n * delta)) * log(MT[0] / M0[0])
    if not np.isfinite(R2_clc):
        return r2_matrix(kAB, kBA, domega, R2_0, vCPMG, tcp)

    return real(R2_clc)

def r2_london(kAB, kBA, domega, R2_0, vCPMG):
    """Compute R2eff using the London slow-exchange approximation.

    Closed-form analytical expression valid when kex << domega::

        Rex = (kex/2) - (1/tcp) * ln(sqrt(D+*cosh^2(eta+) - D-*cos^2(eta-))
                                     + sqrt(D+*sinh^2(eta+) + D-*sin^2(eta-)))

    where tcp = 1/(2*nu_CPMG).  Adds Rex to the intrinsic R2_0.

    Args:
        kAB    (float): A→B exchange rate (s^-1).
        kBA    (float): B→A exchange rate (s^-1).
        domega (float): Chemical shift difference in rad/s.
        R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
        vCPMG  (float): CPMG pulsing frequency (Hz).

    Returns:
        float: Calculated effective R2 relaxation rate (s^-1).
    """
    # Auxiliary variables following London et al. notation
    zet = 2.0 * domega * (kAB - kBA)
    psi = (kAB - kBA) ** 2.0 - domega ** 2 + 4.0 * kAB * kBA
    ss  = sqrt(psi ** 2.0 + zet ** 2.0)

    tcp = 1.0 / (2.0 * vCPMG)  # effective inter-pulse delay for this nu_CPMG

    # eta± and D± factors of the London formula
    np = tcp / sqrt(8.0) * sqrt( psi + ss)
    nm = tcp / sqrt(8.0) * sqrt(-psi + ss)
    Dp = 1.0 / 2.0 * (1.0 + (psi + 2.0 * domega ** 2) / ss)
    Dm = 1.0 / 2.0 * (-1.0 + (psi + 2.0 * domega ** 2) / ss)

    Rex = ((kAB + kBA) / 2.0
           - (1.0 / tcp) * log(sqrt(Dp * cosh(np) ** 2 - Dm * cos(nm) ** 2)
                               + sqrt(Dp * sinh(np) ** 2 + Dm * sin(nm) ** 2)))

    R2_clc = R2_0 + Rex
    if not isfinite(R2_clc):
        return _london_cr72(kAB, kBA, domega, R2_0, vCPMG)
    return R2_clc

def r2_meiboom(kex, phi, R2_0, vCPMG, field):
    """Compute R2eff using the Meiboom fast-exchange approximation.

    Closed-form formula valid when kex >> domega::

        Rex = (phi_omega / kex) * [1 - (4*nu/kex) * tanh(kex / (4*nu))]

    where phi_omega = phi * (2*pi*field)^2 converts phi (in ppm^2) to
    rad^2/s^2 at the given field strength.

    Args:
        kex   (float): Total exchange rate kAB + kBA (s^-1).
        phi   (float): Effective exchange parameter
                       phi = pB*(1-pB)*delta_omega_ppm^2 (ppm^2).
        R2_0  (float): Intrinsic transverse relaxation rate (s^-1).
        vCPMG (float): CPMG pulsing frequency (Hz).
        field (float): Static field strength in MHz (proton Larmor frequency).

    Returns:
        float: Calculated effective R2 relaxation rate (s^-1).
    """
    if kex == 0:
        return R2_0
    # Convert phi from ppm^2 to rad^2/s^2 at this magnetic field
    phi_omega = phi * (2.0 * pi * field) ** 2
    x = kex / (4.0 * vCPMG)
    if abs(x) < 1e-3:
        # Series of the same expression; avoids cancellation at small kex/(4 nu)
        Rex = phi_omega * kex / (48.0 * vCPMG ** 2) * (1.0 - 2.0 * x ** 2 / 5.0 + 17.0 * x ** 4 / 105.0
                                                       - 62.0 * x ** 6 / 945.0)
    else:
        Rex = (phi_omega / kex) * (1.0 - (4.0 * vCPMG / kex) * tanh(kex / (4.0 * vCPMG)))
    R2_clc = R2_0 + Rex
    return R2_clc

def r2_matrix_eig_batch(kAB, kBA, domega_arr, R2_0_arr,
                         v_array, delta_arr, n_arr, tcp):
    """Vectorised R2eff over many residues *and* all CPMG frequencies.

    Evaluates a batch of curves that share the same CPMG frequency grid
    (and tcp).  Per-residue quantities (domega, R2_0)
    vary along the curve axis C; the eigen-decomposition and propagator are
    evaluated for the full (C, V) grid in one set of array operations.

    Falls back to the per-curve implementation for any curve whose
    eigenvectors are degenerate (kAB ≈ 0 or repeated eigenvalues), so the
    numerical result is identical to the scalar path in every regime.

    Args:
        kAB, kBA   (float)     : exchange rate constants (s^-1).
        domega_arr (np.ndarray): chemical shift differences, shape (C,) rad/s.
        R2_0_arr   (np.ndarray): intrinsic R2 per curve, shape (C,) s^-1.
        v_array    (np.ndarray): shared CPMG frequencies, shape (V,) Hz.
        delta_arr  (np.ndarray): 1/(4 v), shape (V,).
        n_arr      (np.ndarray): integer echo counts round(tcp v), shape (V,).
        tcp        (float)     : CPMG pulse delay (s), for the scalar fallback.

    Returns:
        np.ndarray: real R2eff values, shape (C, V).
    """
    domega = np.asarray(domega_arr, dtype=float)      # (C,)
    R2_0   = np.asarray(R2_0_arr,   dtype=float)       # (C,)
    C = domega.shape[0]
    V = v_array.shape[0]

    # No exchange, or no exchange-free denominator: no R2 contribution.
    if kAB == 0 or kAB + kBA == 0:
        return np.repeat(R2_0[:, None], V, axis=1)

    # Curves without chemical shift difference: R2_0, no eigenvectors needed.
    zero_dw = domega == 0
    if np.any(zero_dw):
        out = np.empty((C, V))
        out[zero_dw] = R2_0[zero_dw][:, None]
        nz = ~zero_dw
        if np.any(nz):
            out[nz] = r2_matrix_eig_batch(
                kAB, kBA, domega[nz], R2_0[nz], v_array, delta_arr, n_arr, tcp)
        return out

    # --- Per-curve scalar quantities (independent of vCPMG) ---
    AA = np.sqrt(kBA ** 2 + (2.0 * kAB - 2.0j * domega) * kBA
                 + kAB ** 2 + 2.0j * domega * kAB - domega ** 2)   # (C,)

    L1_base = -(AA + kBA + kAB + (2.0 * R2_0 - 1.0j * domega)) / 2.0  # (C,)
    L2_base =  (AA - kBA - kAB + (1.0j * domega - 2.0 * R2_0)) / 2.0  # (C,)

    denom_v1 = -AA + kBA - kAB - 1.0j * domega        # (C,)
    denom_v2 =  AA + kBA - kAB - 1.0j * domega
    v1 = (2.0 * kAB) / denom_v1
    v2 = (2.0 * kAB) / denom_v2
    diff_v = v2 - v1

    # Degenerate eigenvectors: hand those curves to the scalar path.
    bad = ((np.abs(denom_v1) < 1e-10) | (np.abs(denom_v2) < 1e-10)
           | (np.abs(diff_v) < 1e-10))
    if np.any(bad):
        out = np.empty((C, v_array.shape[0]))
        good = ~bad
        if np.any(good):
            out[good] = r2_matrix_eig_batch(
                kAB, kBA, domega[good], R2_0[good],
                v_array, delta_arr, n_arr, tcp)
        for ci in np.nonzero(bad)[0]:
            out[ci] = [r2_matrix(kAB, kBA, domega[ci], R2_0[ci],
                                          v, tcp) for v in v_array]
        return out

    # Eigenvector matrices X and their inverse Xi, per curve — shape (C,2,2)
    X = np.empty((C, 2, 2), dtype=complex)
    X[:, 0, 0] = 1.0
    X[:, 0, 1] = 1.0
    X[:, 1, 0] = v1
    X[:, 1, 1] = v2
    Xi = np.empty((C, 2, 2), dtype=complex)
    Xi[:, 0, 0] =  v2 / diff_v
    Xi[:, 0, 1] = -1.0 / diff_v
    Xi[:, 1, 0] = -v1 / diff_v
    Xi[:, 1, 1] =  1.0 / diff_v
    Xstar  = X.conj()
    Xistar = Xi.conj()

    # --- (C, V) diagonal exponentials ---
    # Forward propagator eigenvalues: exp(L_base * delta)
    el1 = np.exp(L1_base[:, None] * delta_arr[None, :])   # (C, V)
    el2 = np.exp(L2_base[:, None] * delta_arr[None, :])
    # Conjugate propagator, doubled interval: exp(L_base* * 2 delta)
    el1s2 = np.exp(L1_base.conj()[:, None] * 2.0 * delta_arr[None, :])
    el2s2 = np.exp(L2_base.conj()[:, None] * 2.0 * delta_arr[None, :])

    # eA[c,v] = X[c] diag(el1,el2) Xi[c]; build by scaling X's columns.
    M_el  = np.stack([el1, el2], axis=2)                  # (C, V, 2)
    eA    = (X[:, None, :, :] * M_el[:, :, None, :]) @ Xi[:, None, :, :]
    M_es2 = np.stack([el1s2, el2s2], axis=2)              # (C, V, 2)
    eAstar2 = (Xstar[:, None, :, :] * M_es2[:, :, None, :]) @ Xistar[:, None, :, :]

    # One full CPMG echo, shape (C, V, 2, 2)
    B = eA @ eAstar2 @ eA

    # B^n via the analytic 2x2 matrix-function formula (diagonalisable 2x2)
    a = B[:, :, 0, 0]; b = B[:, :, 0, 1]
    c = B[:, :, 1, 0]; d = B[:, :, 1, 1]
    half_tr = (a + d) * 0.5
    disc    = np.sqrt(((a - d) * 0.5) ** 2 + b * c)
    lb1 = half_tr + disc
    lb2 = half_tr - disc
    denom = lb1 - lb2
    if np.any(np.abs(denom) < 1e-14):
        # Repeated eigenvalue somewhere — fall back per curve for safety
        return np.array([[r2_matrix(kAB, kBA, domega[ci], R2_0[ci],
                                             v, tcp) for v in v_array]
                         for ci in range(C)])
    lb1_n = lb1 ** n_arr[None, :]
    lb2_n = lb2 ** n_arr[None, :]
    f1 = (lb1_n - lb2_n) / denom
    f2 = (lb1 * lb2_n - lb2 * lb1_n) / denom
    Cm = f1[:, :, None, None] * B
    Cm[:, :, 0, 0] += f2
    Cm[:, :, 1, 1] += f2

    pB   = kAB / (kAB + kBA)
    M0_0 = 1.0 - pB
    M0   = np.array([M0_0, pB], dtype=complex)
    MT0  = Cm[:, :, 0, :] @ M0                            # (C, V)

    R2 = (-1.0 / (4.0 * n_arr[None, :] * delta_arr[None, :])) * np.log(MT0 / M0_0)
    out = np.real(R2)
    for ci, vi in zip(*np.nonzero(~np.isfinite(out))):
        out[ci, vi] = r2_matrix(kAB, kBA, domega[ci], R2_0[ci], v_array[vi], tcp)
    return out

def r2_london_vec(kAB, kBA, domega, R2_0, v_array):
    """Vectorised R2eff using the London slow-exchange approximation.

    Args:
        kAB     (float)     : A→B exchange rate (s^-1).
        kBA     (float)     : B→A exchange rate (s^-1).
        domega  (float)     : Chemical shift difference in rad/s.
        R2_0    (float)     : Intrinsic transverse relaxation rate (s^-1).
        v_array (np.ndarray): CPMG pulsing frequencies in Hz, shape (n_v,).

    Returns:
        np.ndarray: Real R2eff values, shape (n_v,).
    """
    zet = 2.0 * domega * (kAB - kBA)
    psi = (kAB - kBA) ** 2.0 - domega ** 2 + 4.0 * kAB * kBA
    ss  = sqrt(psi ** 2.0 + zet ** 2.0)

    tcp = 1.0 / (2.0 * v_array)   # effective inter-pulse delay, shape (n_v,)

    # Guard against negative argument under sqrt due to floating-point error
    eta_p = tcp / sqrt(8.0) * np.sqrt(np.maximum(0.0,  psi + ss))
    eta_m = tcp / sqrt(8.0) * np.sqrt(np.maximum(0.0, -psi + ss))

    Dp = 0.5 * (1.0 + (psi + 2.0 * domega ** 2) / ss)
    Dm = 0.5 * (-1.0 + (psi + 2.0 * domega ** 2) / ss)

    Rex = ((kAB + kBA) / 2.0
           - (1.0 / tcp) * np.log(
               np.sqrt(Dp * np.cosh(eta_p) ** 2 - Dm * np.cos(eta_m) ** 2)
               + np.sqrt(Dp * np.sinh(eta_p) ** 2 + Dm * np.sin(eta_m) ** 2)))
    out = R2_0 + Rex
    for i in np.nonzero(~np.isfinite(out))[0]:
        out[i] = _london_cr72(kAB, kBA, domega, R2_0, v_array[i])
    return out

def r2_meiboom_vec(kex, phi, R2_0, v_array, field):
    """Vectorised R2eff using the Meiboom fast-exchange approximation.

    Args:
        kex     (float)     : Total exchange rate kAB + kBA (s^-1).
        phi     (float)     : phi = pB*(1-pB)*delta_omega_ppm^2 (ppm^2).
        R2_0    (float)     : Intrinsic transverse relaxation rate (s^-1).
        v_array (np.ndarray): CPMG pulsing frequencies in Hz, shape (n_v,).
        field   (float)     : Static field strength in MHz.

    Returns:
        np.ndarray: Real R2eff values, shape (n_v,).
    """
    if kex == 0:
        return R2_0 + np.zeros_like(v_array, dtype=float)
    phi_omega = phi * (2.0 * pi * field) ** 2
    x = kex / (4.0 * v_array)
    Rex = (phi_omega / kex) * (1.0 - (4.0 * v_array / kex) * np.tanh(x))
    small = np.abs(x) < 1e-3
    if np.any(small):
        series = (phi_omega * kex / (48.0 * v_array ** 2)
                  * (1.0 - 2.0 * x ** 2 / 5.0 + 17.0 * x ** 4 / 105.0 - 62.0 * x ** 6 / 945.0))
        Rex = np.where(small, series, Rex)
    return R2_0 + Rex
