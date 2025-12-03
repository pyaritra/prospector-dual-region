"""
Prospector Dual-Region SED Fitting

A custom extension for Prospector that enables simultaneous SED fitting
of two spatially resolved regions with shared (blended) photometry.
"""

__version__ = "0.1.0"
__author__ = "Your Name"

from .dual_region_fitting import (
    DualRegionLikelihood,
    build_observations_dual_region,
    build_dual_models,
    run_dual_region_fit,
)

__all__ = [
    "DualRegionLikelihood",
    "build_observations_dual_region",
    "build_dual_models",
    "run_dual_region_fit",
]
