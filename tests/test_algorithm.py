from __future__ import annotations

import numpy as np

from echo.algorithm import BayesianChangePointDetector


def test_bpsd_change_point_detector_identifies_depth_ratio_shift() -> None:
    """BPSD should localize a sustained structural shift in ordered PDS signal."""

    positions = np.arange(42128780, 42128820)
    signal = np.array([1.0] * 20 + [1.55] * 20)
    detector = BayesianChangePointDetector(min_segment_size=5, log_bayes_factor_threshold=3.0)

    calls = detector.detect(positions, signal)

    assert calls
    assert abs(calls[0].coordinate - 42128800) <= 1
