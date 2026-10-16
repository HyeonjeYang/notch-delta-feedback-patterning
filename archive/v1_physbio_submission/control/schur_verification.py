"""Exact algebraic verification that the Schur-complement characteristic
polynomial equals det(sI - J_mu) directly (the 4x4 Jacobian's own
eigenvalues). Logic copied unmodified from the source repository's
``preprint_final_analysis/scripts/verify_schur.py``.

Requires sympy (symbolic adjugate/determinant expansion) -- used only here,
nowhere else in this repository, and not required for any figure script.
"""

from __future__ import annotations

import numpy as np


def full_char_poly_coeffs(J_mu: np.ndarray) -> np.ndarray:
    """Direct characteristic-polynomial coefficients of the 4x4 Jacobian."""
    return np.poly(J_mu)


def exact_quartic_coeffs_via_schur(A_mu: np.ndarray, b: np.ndarray, c: np.ndarray, a_S: float):
    """Build the exact quartic char-poly coefficients from the
    Schur-complement construction via sympy, using the adjugate (not the
    inverse) of (sI-A_mu) so the expansion stays valid even near A_mu's own
    eigenvalues."""
    import sympy as sp

    s = sp.symbols("s")
    A_sym = sp.Matrix(A_mu.tolist())
    b_sym = sp.Matrix(b.tolist())
    c_sym = sp.Matrix([c.tolist()])
    M = s * sp.eye(3) - A_sym
    adj = M.adjugate()
    det_A = M.det()
    quartic = sp.expand((s - a_S) * det_A - (c_sym * adj * b_sym)[0, 0])
    poly = sp.Poly(quartic, s)
    return [complex(coef) for coef in poly.all_coeffs()]


def verify(J_mu: np.ndarray) -> dict:
    """Compare the direct and Schur-complement characteristic polynomials
    for one Jacobian block; returns max coefficient/root mismatch."""
    from control.state_space import partition

    A_mu, b, c, a_S = partition(J_mu)
    direct = np.array(full_char_poly_coeffs(J_mu), dtype=complex)
    schur = np.array(exact_quartic_coeffs_via_schur(A_mu, b, c, a_S), dtype=complex)
    direct = direct / direct[0]
    schur = schur / schur[0]
    max_coeff_diff = float(np.max(np.abs(direct - schur)))

    direct_roots = np.roots(direct)
    schur_roots = list(np.roots(schur))
    root_diffs = []
    for r in direct_roots:
        idx = int(np.argmin([abs(r - u) for u in schur_roots]))
        root_diffs.append(abs(r - schur_roots.pop(idx)))
    max_root_mismatch = float(max(root_diffs))

    return {"max_coeff_diff": max_coeff_diff, "max_root_mismatch": max_root_mismatch,
           "direct_eigenvalues": direct_roots}
