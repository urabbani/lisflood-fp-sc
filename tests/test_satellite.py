"""
Tests for satellite data preprocessing module.
"""

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from typing import Dict, Any, Tuple

import numpy as np
import rasterio

from satellite.preprocessor import SatelliteInundation


class TestMNDWI:
    """Tests for MNDWI water index calculation"""

    def setup_method(self):
        self.preprocessor = SatelliteInundation()

    def test_mndwi_basic(self):
        green = np.array([[0.5, 0.3], [0.6, 0.4]], dtype=np.float64)
        swir = np.array([[0.2, 0.1], [0.3, 0.2]], dtype=np.float64)
        result = self.preprocessor.mndwi_water_index(green, swir)
        assert result.shape == (2, 2)
        # MNDWI returns binary mask (values > threshold)
        assert result.dtype == np.uint8

    def test_mndwi_division_by_zero(self):
        """Both green and swir are zero → should return 0 not NaN"""
        green = np.array([[0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
        swir = np.array([[0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
        result = self.preprocessor.mndwi_water_index(green, swir)
        assert result.shape == (2, 2)
        assert not np.any(np.isnan(result))
        assert not np.any(np.isinf(result))

    def test_mndwi_partial_zero(self):
        """Some cells have zero denominator"""
        green = np.array([[0.5, 0.0], [0.3, 0.4]], dtype=np.float64)
        swir = np.array([[0.2, 0.0], [0.1, 0.1]], dtype=np.float64)
        result = self.preprocessor.mndwi_water_index(green, swir)
        assert not np.any(np.isnan(result))
        assert not np.any(np.isinf(result))

    def test_mndwi_threshold(self):
        """Default threshold -0.1 should classify correctly"""
        green = np.array([[0.3, 0.5], [0.1, 0.4]], dtype=np.float64)
        swir = np.array([[0.5, 0.2], [0.2, 0.3]], dtype=np.float64)
        result = self.preprocessor.mndwi_water_index(green, swir, threshold=-0.1)
        assert result.dtype == np.uint8
        assert set(np.unique(result)).issubset({0, 1})


class TestSARBackscatter:
    """Tests for SAR backscatter thresholding"""

    def setup_method(self):
        self.preprocessor = SatelliteInundation()

    def test_sar_basic(self):
        vv = np.array([[-20.0, -18.0], [-5.0, 0.0]], dtype=np.float64)
        result = self.preprocessor.sar_backscatter_threshold(vv, threshold=-15.0)
        assert result.shape == (2, 2)
        assert result.dtype == np.uint8
        # -20 and -18 are below -15 threshold → water (1)
        assert result[0, 0] == 1
        assert result[0, 1] == 1
        # -5 and 0 are above -15 threshold → not water (0)
        assert result[1, 0] == 0
        assert result[1, 1] == 0

    def test_sar_no_water(self):
        vv = np.array([[0.0, 5.0], [10.0, 15.0]], dtype=np.float64)
        result = self.preprocessor.sar_backscatter_threshold(vv, threshold=-15.0)
        assert np.sum(result) == 0


class TestThresholdDepth:
    """Tests for depth thresholding"""

    def setup_method(self):
        self.preprocessor = SatelliteInundation()

    def test_threshold_basic(self):
        depth = np.array([[0.5, 0.05], [0.3, 0.0]], dtype=np.float64)
        result = self.preprocessor.threshold_depth(depth, threshold=0.1)
        assert result[0, 0] == 1  # 0.5 > 0.1
        assert result[0, 1] == 0  # 0.05 <= 0.1
        assert result[1, 0] == 1  # 0.3 > 0.1
        assert result[1, 1] == 0  # 0.0 <= 0.1

    def test_threshold_custom(self):
        depth = np.array([[0.5, 0.3]], dtype=np.float64)
        result = self.preprocessor.threshold_depth(depth, threshold=0.4)
        assert result[0, 0] == 1
        assert result[0, 1] == 0


class TestSmoothMask:
    """Tests for mask smoothing"""

    def setup_method(self):
        self.preprocessor = SatelliteInundation()

    def test_smooth_removes_small_patches(self):
        # Create mask with a single-pixel patch
        mask = np.zeros((20, 20), dtype=np.uint8)
        mask[5, 5] = 1  # Single pixel
        result = self.preprocessor.smooth_mask(mask, min_area=10)
        assert result[5, 5] == 0  # Should be removed

    def test_smooth_keeps_large_patches(self):
        mask = np.zeros((20, 20), dtype=np.uint8)
        mask[5:15, 5:15] = 1  # 100 pixels
        result = self.preprocessor.smooth_mask(mask, min_area=10)
        assert np.sum(result) > 0  # Should keep most of it


class TestValidateMask:
    """Tests for mask validation"""

    def setup_method(self):
        self.preprocessor = SatelliteInundation()

    def test_validate_empty_mask(self):
        mask = np.zeros((10, 10), dtype=np.uint8)
        stats = self.preprocessor.validate_mask(mask)
        assert stats["flooded_pixels"] == 0
        assert stats["flood_percentage"] == 0.0

    def test_validate_full_mask(self):
        mask = np.ones((10, 10), dtype=np.uint8)
        stats = self.preprocessor.validate_mask(mask)
        assert stats["flooded_pixels"] == 100
        assert stats["flood_percentage"] == 100.0

    def test_validate_partial_mask(self):
        mask = np.zeros((10, 10), dtype=np.uint8)
        mask[:5, :5] = 1  # 25 pixels
        stats = self.preprocessor.validate_mask(mask)
        assert stats["flooded_pixels"] == 25
        assert stats["flood_percentage"] == 25.0
