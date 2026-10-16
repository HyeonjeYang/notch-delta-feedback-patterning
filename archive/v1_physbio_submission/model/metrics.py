"""Spatial order-parameter diagnostics, copied unmodified from the source
repository's ``core.py::morans_i`` and
``amplitude_metrics.py::amplitude_order_metrics``.

``q_amp = -Moran(D) * CV(D)^2`` is amplitude-sensitive (degree-2
homogeneous in the deviation field); Moran's I alone is amplitude-blind
(scale-invariant) -- see docs/CONTROL_THEORY_DERIVATION.md's note on Figure
5C (Moran's I pinned at exactly -0.5 throughout growth and saturation).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from geometry.lattice import hex_adj


def morans_i(arr: np.ndarray) -> float:
    """Moran's I spatial autocorrelation on the six-neighbour contact
    graph. Negative values indicate nearest-neighbour alternation; -1 is
    not attainable on this non-bipartite (triangular) graph."""
    arr = np.asarray(arr, dtype=np.float64)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        raise ValueError("arr must be a non-empty finite array")
    z = arr - arr.mean()
    den = (z ** 2).sum()
    if den == 0:
        return 0.0
    return float((arr.size / (arr.size * 6)) * (z * hex_adj(z)).sum() / den)


def order_parameters(field: Any) -> dict[str, float]:
    """Amplitude and amplitude-weighted order for one positive field
    (typically the Delta field D)."""
    values = np.asarray(field, dtype=float)
    if values.ndim != 2 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("field must be a non-empty finite 2-D array")
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if mean == 0.0:
        cv = q_amp = contrast = float("nan")
    else:
        cv = float(sd / abs(mean))
        q_amp = float(-morans_i(values) * cv ** 2)
        high = values[values > mean]
        low = values[values <= mean]
        contrast = (float((high.mean() - low.mean()) / abs(mean))
                   if high.size and low.size else 0.0)
    return {
        "mean": mean, "std": sd, "cv": cv,
        "relative_contrast": contrast,
        "moran": float(morans_i(values)),
        "q_amp": q_amp,
    }
