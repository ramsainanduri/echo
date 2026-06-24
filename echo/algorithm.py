"""Bayesian Paralog Signal Deconvolution implementation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import gammaln


@dataclass(frozen=True)
class ChangePoint:
    """A Bayesian change-point call."""

    left_index: int
    right_index: int
    coordinate: int
    log_bayes_factor: float
    left_mean: float
    right_mean: float


class BayesianChangePointDetector:
    """Conjugate Bayesian change-point detector for Gaussian depth ratios.

    For an ordered signal ``x`` this detector compares a one-segment model
    against every two-segment split. Each segment mean and variance are
    integrated out under a Normal-Inverse-Gamma prior:

    ``sigma^2 ~ InvGamma(alpha0, beta0)``
    ``mu | sigma^2 ~ Normal(mu0, sigma^2 / kappa0)``

    The log marginal likelihood for a segment with ``n`` observations is:

    ``lgamma(alpha_n) - lgamma(alpha0) + 0.5(log(kappa0)-log(kappa_n))``
    ``+ alpha0 log(beta0) - alpha_n log(beta_n) - n/2 log(pi)``
    """

    def __init__(
        self,
        min_segment_size: int = 8,
        log_bayes_factor_threshold: float = 8.0,
        mu0: float = 1.0,
        kappa0: float = 0.05,
        alpha0: float = 2.0,
        beta0: float = 0.05,
    ) -> None:
        self.min_segment_size = min_segment_size
        self.log_bayes_factor_threshold = log_bayes_factor_threshold
        self.mu0 = mu0
        self.kappa0 = kappa0
        self.alpha0 = alpha0
        self.beta0 = beta0

    def detect(self, positions: np.ndarray, signal: np.ndarray) -> list[ChangePoint]:
        """Detect one or more change points in an ordered signal.

        Parameters
        ----------
        positions
            Genomic positions corresponding to the signal.
        signal
            Normalized PDS depth ratios.

        Returns
        -------
        list of ChangePoint
            Recursive split calls passing the Bayes-factor threshold.
        """

        clean_mask = np.isfinite(signal)
        clean_positions = positions[clean_mask]
        clean_signal = signal[clean_mask]
        if clean_signal.size < self.min_segment_size * 2:
            return []
        return self._detect_recursive(clean_positions, clean_signal, offset=0)

    def _detect_recursive(
        self, positions: np.ndarray, signal: np.ndarray, offset: int
    ) -> list[ChangePoint]:
        best = self._best_split(positions, signal, offset)
        if best is None or best.log_bayes_factor < self.log_bayes_factor_threshold:
            return []
        split = best.right_index - offset
        left = self._detect_recursive(positions[:split], signal[:split], offset)
        right = self._detect_recursive(positions[split:], signal[split:], best.right_index)
        return [*left, best, *right]

    def _best_split(
        self, positions: np.ndarray, signal: np.ndarray, offset: int
    ) -> ChangePoint | None:
        if signal.size < self.min_segment_size * 2:
            return None
        null_logp = self._log_marginal_likelihood(signal)
        best_cp: ChangePoint | None = None
        for split in range(self.min_segment_size, signal.size - self.min_segment_size + 1):
            left = signal[:split]
            right = signal[split:]
            alt_logp = self._log_marginal_likelihood(left) + self._log_marginal_likelihood(right)
            log_bf = alt_logp - null_logp
            if best_cp is None or log_bf > best_cp.log_bayes_factor:
                best_cp = ChangePoint(
                    left_index=offset + split - 1,
                    right_index=offset + split,
                    coordinate=int(positions[split]),
                    log_bayes_factor=float(log_bf),
                    left_mean=float(np.mean(left)),
                    right_mean=float(np.mean(right)),
                )
        return best_cp

    def _log_marginal_likelihood(self, values: np.ndarray) -> float:
        values = values[np.isfinite(values)]
        n = values.size
        if n == 0:
            return float("-inf")
        mean = float(np.mean(values))
        sumsq = float(np.sum((values - mean) ** 2))
        kappa_n = self.kappa0 + n
        alpha_n = self.alpha0 + n / 2.0
        beta_n = self.beta0 + 0.5 * sumsq
        beta_n += (self.kappa0 * n * (mean - self.mu0) ** 2) / (2.0 * kappa_n)
        return float(
            gammaln(alpha_n)
            - gammaln(self.alpha0)
            + 0.5 * (np.log(self.kappa0) - np.log(kappa_n))
            + self.alpha0 * np.log(self.beta0)
            - alpha_n * np.log(beta_n)
            - (n / 2.0) * np.log(np.pi)
        )


class BPSDAlgorithm:
    """Bayesian Paralog Signal Deconvolution.

    This class normalizes test-sample PDS depths against the PON PDS baseline,
    then applies Bayesian change-point detection to the ordered PDS ratio.
    """

    def __init__(self, detector: BayesianChangePointDetector | None = None) -> None:
        self.detector = detector or BayesianChangePointDetector()

    def deconvolve(
        self, sample_pds: pd.DataFrame, pon_pds: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[ChangePoint]]:
        """Normalize PDS depth signal and call structural shifts.

        Parameters
        ----------
        sample_pds
            Sample PDS signal with ``chrom``, ``pos0``, and ``depth``.
        pon_pds
            PON PDS signal with ``chrom``, ``pos0``, ``mean_depth``, and
            ``sd_depth``.

        Returns
        -------
        tuple
            Normalized signal dataframe and change-point calls.
        """

        if sample_pds.empty or pon_pds.empty:
            return pd.DataFrame(), []
        merged = sample_pds.merge(pon_pds, on=["chrom", "pos0"], how="inner")
        if merged.empty:
            return pd.DataFrame(), []
        eps = np.finfo(float).eps
        merged["pds_ratio"] = merged["depth"] / np.maximum(merged["mean_depth"], eps)
        merged["z_score"] = (merged["depth"] - merged["mean_depth"]) / np.maximum(
            merged["sd_depth"], eps
        )
        merged = merged.sort_values(["chrom", "pos0"]).reset_index(drop=True)
        breakpoints = self.detector.detect(
            merged["pos0"].to_numpy(dtype=np.int64),
            merged["pds_ratio"].to_numpy(dtype=float),
        )
        return merged, breakpoints
