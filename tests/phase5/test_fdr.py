"""Benjamini-Hochberg FDR tests."""

from __future__ import annotations

import numpy as np

from signflow.phase5.nulls import benjamini_hochberg


def test_bh_fdr_monotone_and_bounded():
    p = np.array([0.001, 0.02, 0.03, 0.5])
    q = benjamini_hochberg(p)
    assert np.all(q >= 0) and np.all(q <= 1)
    assert q[0] <= q[-1] or np.allclose(sorted(q), q)
