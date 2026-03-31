"""
Tests for calibration metrics (NSE, KGE, IoU, F1, composite_score).
"""

import numpy as np
import pytest
from calibration.metrics import CalibrationMetrics


class TestNashSutcliffe:
    """Tests for Nash-Sutcliffe Efficiency (NSE)"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_perfect_match(self):
        obs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        assert self.metrics.nash_sutcliffe(obs, obs) == pytest.approx(1.0)

    def test_mean_prediction(self):
        obs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        sim = np.full(5, np.mean(obs))
        assert self.metrics.nash_sutcliffe(obs, sim) == pytest.approx(0.0)

    def test_poor_prediction(self):
        obs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        sim = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        assert self.metrics.nash_sutcliffe(obs, sim) < 0.0

    def test_nan_handling(self):
        obs = np.array([1.0, np.nan, 3.0, 4.0])
        sim = np.array([1.0, 2.0, np.nan, 4.0])
        # Only indices 0 and 3 are valid in both
        nse = self.metrics.nash_sutcliffe(obs, sim)
        assert not np.isnan(nse)

    def test_all_nan(self):
        obs = np.array([np.nan, np.nan])
        sim = np.array([np.nan, np.nan])
        assert self.metrics.nash_sutcliffe(obs, sim) == -np.inf

    def test_constant_observation(self):
        obs = np.array([5.0, 5.0, 5.0])
        sim = np.array([5.0, 5.0, 5.0])
        # denominator = 0, numerator = 0 → perfect match
        assert self.metrics.nash_sutcliffe(obs, sim) == 1.0

    def test_constant_observation_bad_sim(self):
        obs = np.array([5.0, 5.0, 5.0])
        sim = np.array([1.0, 2.0, 3.0])
        # denominator = 0, numerator > 0 → -inf
        assert self.metrics.nash_sutcliffe(obs, sim) == -np.inf


class TestKlingGupta:
    """Tests for Kling-Gupta Efficiency (KGE)"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_perfect_match(self):
        obs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        assert self.metrics.kling_gupta(obs, obs) == pytest.approx(1.0)

    def test_poor_prediction(self):
        obs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        sim = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        assert self.metrics.kling_gupta(obs, sim) < 1.0

    def test_nan_handling(self):
        obs = np.array([1.0, np.nan, 3.0, 4.0])
        sim = np.array([1.0, 2.0, np.nan, 4.0])
        kge = self.metrics.kling_gupta(obs, sim)
        assert not np.isnan(kge)

    def test_all_nan(self):
        obs = np.array([np.nan, np.nan])
        sim = np.array([np.nan, np.nan])
        assert self.metrics.kling_gupta(obs, sim) == -np.inf


class TestIoU:
    """Tests for Intersection over Union (IoU)"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_perfect_overlap(self):
        mask = np.ones((10, 10), dtype=np.uint8)
        assert self.metrics.iou(mask, mask) == pytest.approx(1.0)

    def test_no_overlap(self):
        obs = np.zeros((10, 10), dtype=np.uint8)
        sim = np.ones((10, 10), dtype=np.uint8)
        assert self.metrics.iou(obs, sim) == pytest.approx(0.0)

    def test_partial_overlap(self):
        obs = np.zeros((10, 10), dtype=np.uint8)
        sim = np.zeros((10, 10), dtype=np.uint8)
        obs[:5, :5] = 1
        sim[:5, :5] = 1
        sim[5:, 5:] = 1
        # intersection = 25, union = 25 + 25 = 50
        assert self.metrics.iou(obs, sim) == pytest.approx(0.5)

    def test_both_empty(self):
        mask = np.zeros((10, 10), dtype=np.uint8)
        assert self.metrics.iou(mask, mask) == 0.0

    def test_shape_mismatch(self):
        obs = np.ones((10, 10), dtype=np.uint8)
        sim = np.ones((8, 8), dtype=np.uint8)
        # Should crop to (8, 8) and compute IoU
        assert self.metrics.iou(obs, sim) == pytest.approx(1.0)

    def test_depth_grid_thresholding(self):
        obs = np.array([[0.5, 0.05], [0.5, 0.0]], dtype=np.float64)
        sim = np.array([[0.5, 0.05], [0.5, 0.0]], dtype=np.float64)
        # With threshold=0.1: obs = [[1,0],[1,0]], sim = [[1,0],[1,0]]
        assert self.metrics.iou(obs, sim, depth_threshold=0.1) == pytest.approx(1.0)


class TestF1Score:
    """Tests for F1 Score"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_perfect_match(self):
        mask = np.ones((10, 10), dtype=np.uint8)
        assert self.metrics.f1_score(mask, mask) == pytest.approx(1.0)

    def test_no_true_positives(self):
        obs = np.zeros((10, 10), dtype=np.uint8)
        sim = np.ones((10, 10), dtype=np.uint8)
        assert self.metrics.f1_score(obs, sim) == 0.0


class TestCompositeScore:
    """Tests for composite score"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_perfect_metrics(self):
        metrics = {"nse": 1.0, "kge": 1.0, "iou": 1.0}
        weights = {"nse": 0.4, "kge": 0.3, "iou": 0.3}
        assert self.metrics.composite_score(metrics, weights) == pytest.approx(1.0)

    def test_negative_metrics_clamped(self):
        metrics = {"nse": -5.0, "kge": -2.0, "iou": 0.5}
        weights = {"nse": 0.4, "kge": 0.3, "iou": 0.3}
        score = self.metrics.composite_score(metrics, weights)
        # NSE and KGE clamped to 0, only IoU contributes
        assert score == pytest.approx(0.3 * 0.5)

    def test_missing_weight(self):
        metrics = {"nse": 0.9, "iou": 0.8}
        weights = {"nse": 0.5}
        score = self.metrics.composite_score(metrics, weights)
        # iou not in weights, so only nse contributes
        assert score == pytest.approx(0.5 * 0.9)


class TestTimingError:
    """Tests for timing error (issue #7 fix)"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_same_peak_time(self):
        obs = np.array([1.0, 5.0, 3.0])
        sim = np.array([2.0, 6.0, 2.0])
        ts = np.array([0.0, 1.0, 2.0])
        assert self.metrics.timing_error(obs, sim, ts) == pytest.approx(0.0)

    def test_different_peak_time(self):
        obs = np.array([1.0, 5.0, 3.0])
        sim = np.array([3.0, 2.0, 6.0])
        ts = np.array([0.0, 1.0, 2.0])
        assert self.metrics.timing_error(obs, sim, ts) == pytest.approx(1.0)

    def test_with_nan(self):
        obs = np.array([1.0, np.nan, 5.0, 3.0])
        sim = np.array([2.0, 3.0, np.nan, 6.0])
        ts = np.array([0.0, 1.0, 2.0, 3.0])
        # After mask: obs=[1.0, 3.0], sim=[2.0, 6.0], ts=[0.0, 3.0]
        # obs peak at idx 1 (3.0), sim peak at idx 1 (6.0) → same time
        error = self.metrics.timing_error(obs, sim, ts)
        assert error == pytest.approx(0.0)

    def test_all_nan(self):
        obs = np.array([np.nan, np.nan])
        sim = np.array([np.nan, np.nan])
        ts = np.array([0.0, 1.0])
        assert self.metrics.timing_error(obs, sim, ts) == np.inf


class TestValidateMetrics:
    """Tests for metric validation"""

    def setup_method(self):
        self.metrics = CalibrationMetrics()

    def test_valid_metrics(self):
        metrics = {"nse": 0.85, "kge": 0.80, "iou": 0.70}
        valid, msg = self.metrics.validate_metrics(metrics)
        assert valid is True

    def test_nan_metric(self):
        metrics = {"nse": np.nan, "kge": 0.80, "iou": 0.70}
        valid, msg = self.metrics.validate_metrics(metrics)
        assert valid is False

    def test_nse_above_one(self):
        metrics = {"nse": 1.05, "kge": 0.80, "iou": 0.70}
        valid, msg = self.metrics.validate_metrics(metrics)
        assert valid is False

    def test_iou_out_of_range(self):
        metrics = {"nse": 0.85, "kge": 0.80, "iou": 1.5}
        valid, msg = self.metrics.validate_metrics(metrics)
        assert valid is False
