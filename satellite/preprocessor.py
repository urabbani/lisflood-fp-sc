"""
Satellite data preprocessing and alignment.
Handles flood mask extraction and grid alignment.
"""

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from typing import Dict, Any, Tuple


class SatelliteInundation:
    """Process and align satellite-derived inundation data"""
    
    def __init__(self, model_crs: str = "EPSG:32642",
                 model_resolution: float = 30.0,
                 model_bounds: Tuple[float, float, float, float] = None):
        """
        Initialize satellite preprocessor.
        
        Args:
            model_crs: Coordinate reference system (default: UTM Zone 42N)
            model_resolution: Grid resolution in meters
            model_bounds: (min_x, max_x, min_y, max_y) in model CRS
        """
        self.model_crs = model_crs
        self.model_resolution = model_resolution
        self.model_bounds = model_bounds
    
    def load_flood_mask(self, file_path: str) -> np.ndarray:
        """
        Load satellite flood mask or depth grid.
        
        Args:
            file_path: Path to GeoTIFF file
        
        Returns:
            2D numpy array
        """
        with rasterio.open(file_path) as src:
            mask = src.read(1)
            self.src_crs = src.crs
            self.src_transform = src.transform
            self.src_bounds = src.bounds
            self.src_shape = mask.shape
        
        print(f"   Loaded flood mask: {mask.shape}")
        return mask
    
    def align_to_model_grid(self, satellite_mask: np.ndarray,
                          model_extent: Tuple[float, float, float, float] = None) -> np.ndarray:
        """
        Reproject satellite mask to model grid.
        
        Args:
            satellite_mask: Original satellite flood mask
            model_extent: (min_x, max_x, min_y, max_y) in model CRS
        
        Returns:
            Aligned mask
        """
        extent = model_extent or self.model_bounds
        
        if extent is None:
            raise ValueError("Model extent must be provided")
        
        # Calculate target grid dimensions
        x_range = extent[1] - extent[0]
        y_range = extent[3] - extent[2]
        n_cols = int(x_range / self.model_resolution)
        n_rows = int(y_range / self.model_resolution)
        
        # Create target transform
        dst_transform = from_origin(
            extent[1],  # top-right x
            extent[3],  # top-right y
            self.model_resolution,
            self.model_resolution
        )
        
        # Create destination array
        aligned_mask = np.zeros((n_rows, n_cols), dtype=satellite_mask.dtype)
        
        # Reproject
        reproject(
            source=satellite_mask,
            destination=aligned_mask,
            src_transform=self.src_transform,
            src_crs=self.src_crs,
            dst_transform=dst_transform,
            dst_crs=self.model_crs,
            resampling=Resampling.nearest
        )
        
        print(f"   Aligned to model grid: {aligned_mask.shape}")
        return aligned_mask
    
    def threshold_depth(self, depth_grid: np.ndarray,
                       threshold: float = 0.1) -> np.ndarray:
        """
        Convert depth grid to binary flood mask.
        
        Args:
            depth_grid: 2D array of water depths (meters)
            threshold: Minimum depth to consider flooded (meters)
        
        Returns:
            Binary mask (0=dry, 1=flooded)
        """
        mask = (depth_grid > threshold).astype(np.uint8)
        
        flooded_pixels = np.sum(mask)
        total_pixels = mask.size
        flood_percentage = (flooded_pixels / total_pixels) * 100
        
        print(f"   Thresholded at {threshold}m: {flood_percentage:.1f}% flooded")
        
        return mask
    
    def mndwi_water_index(self, green_band: np.ndarray,
                          swir_band: np.ndarray,
                          threshold: float = -0.1) -> np.ndarray:
        """
        Compute Modified Normalized Difference Water Index (MNDWI).
        
        MNDWI = (Green - SWIR) / (Green + SWIR)
        
        Args:
            green_band: Green reflectance band
            swir_band: SWIR reflectance band
            threshold: MNDWI threshold for water classification
        
        Returns:
            Binary water mask
        """
        denominator = green_band + swir_band
        with np.errstate(divide='ignore', invalid='ignore'):
            mndwi = np.where(denominator != 0, (green_band - swir_band) / denominator, 0.0)
        water_mask = (mndwi > threshold).astype(np.uint8)
        
        print(f"   MNDWI range: {mndwi.min():.3f} to {mndwi.max():.3f}")
        
        return water_mask
    
    def sar_backscatter_threshold(self, vv_band: np.ndarray,
                                threshold: float = -15.0) -> np.ndarray:
        """
        Extract water from Sentinel-1 SAR backscatter.
        
        Water typically has low backscatter (dark in SAR images).
        
        Args:
            vv_band: VV polarized backscatter (dB)
            threshold: Backscatter threshold (dB, lower = more water)
        
        Returns:
            Binary water mask
        """
        water_mask = (vv_band < threshold).astype(np.uint8)
        
        print(f"   SAR backscatter range: {vv_band.min():.1f} to {vv_band.max():.1f} dB")
        print(f"   Threshold: {threshold} dB")
        
        return water_mask
    
    def smooth_mask(self, mask: np.ndarray, 
                   iterations: int = 1,
                   min_area: int = 10) -> np.ndarray:
        """
        Smooth flood mask by removing small isolated patches.
        
        Args:
            mask: Binary flood mask
            iterations: Number of morphological operations
            min_area: Minimum area (in pixels) to keep
        
        Returns:
            Smoothed mask
        """
        from scipy import ndimage
        
        smoothed = mask.copy()
        
        # Morphological closing (fill small holes)
        for i in range(iterations):
            smoothed = ndimage.binary_closing(smoothed, iterations=1)
            smoothed = ndimage.binary_opening(smoothed, iterations=1)
        
        # Remove small patches
        labeled, num_features = ndimage.label(smoothed)
        
        for label in range(1, num_features + 1):
            area = np.sum(labeled == label)
            if area < min_area:
                smoothed[labeled == label] = 0
        
        return smoothed.astype(np.uint8)
    
    def validate_mask(self, mask: np.ndarray) -> Dict[str, Any]:
        """
        Validate flood mask quality.
        
        Args:
            mask: Binary flood mask
        
        Returns:
            Dict with validation metrics
        """
        stats = {
            "total_pixels": mask.size,
            "flooded_pixels": np.sum(mask),
            "flood_percentage": (np.sum(mask) / mask.size) * 100,
            "mean_patch_size": 0,
            "edge_pixels": 0
        }
        
        # Calculate mean patch size
        from scipy import ndimage
        labeled, num_features = ndimage.label(mask)
        
        if num_features > 0:
            patch_sizes = [np.sum(labeled == i) for i in range(1, num_features + 1)]
            stats["mean_patch_size"] = np.mean(patch_sizes)
            stats["num_patches"] = num_features
        
        # Calculate edge pixels (vertical + horizontal transitions separately)
        vertical_edges = np.sum(mask[1:, :] != mask[:-1, :])
        horizontal_edges = np.sum(mask[:, 1:] != mask[:, :-1])
        stats["edge_pixels"] = vertical_edges + horizontal_edges
        
        print(f"   Flood percentage: {stats['flood_percentage']:.1f}%")
        print(f"   Number of patches: {stats.get('num_patches', 0)}")
        print(f"   Mean patch size: {stats.get('mean_patch_size', 0):.1f} pixels")
        
        return stats
