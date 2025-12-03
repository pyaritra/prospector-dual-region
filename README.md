# Prospector Dual-Region SED Fitting
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.17805754)

A custom extension for [Prospector](https://github.com/bd-j/prospector) that enables simultaneous SED fitting of two spatially resolved regions with shared (blended) photometry in some bands.

**Use case**: When you have high-resolution photometry (e.g., HST) that resolves two regions separately, but low-resolution photometry (e.g., WISE) where both regions are blended together.


- Simultaneous fitting of two independent stellar population models
- Adds constraint that blended bands equal sum of model predictions

## Installation

### Prerequisites

First, install Prospector and its dependencies:

```bash
# Create conda environment
conda create -n prospector python=3.9
conda activate prospector

# Install dependencies
pip install astro-prospector
pip install python-fsps
pip install dynesty sedpy h5py astropy

# Install FSPS stellar population models
# Follow instructions at: https://dfm.io/python-fsps/current/installation/
```

### Install This Package

```bash
git clone https://github.com/YOUR_USERNAME/prospector-dual-region.git
cd prospector-dual-region
pip install -e .
```

Or just copy the files directly into your working directory.

## Quick Start

```python
from prospector_dual_region import run_dual_region_fit
import numpy as np

# Your HST photometry for Region A
hst_filters_a = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w']
hst_maggies_a = np.array([1.5e-6, 2.3e-6])  # flux in maggies (Jy/3631)
hst_unc_a = np.array([1.5e-7, 2.3e-7])

# Your HST photometry for Region B
hst_filters_b = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w']
hst_maggies_b = np.array([0.8e-6, 1.1e-6])
hst_unc_b = np.array([1.0e-7, 1.5e-7])

# Your WISE photometry (blended - sum of both regions)
wise_filters = ['wise_w1', 'wise_w2']
wise_maggies = np.array([5.0e-6, 4.5e-6])  # Total observed flux
wise_unc = np.array([5.0e-7, 5.0e-7])

# Run the fit
output = run_dual_region_fit(
    hst_filters_a, hst_maggies_a, hst_unc_a,
    hst_filters_b, hst_maggies_b, hst_unc_b,
    wise_filters, wise_maggies, wise_unc,
    redshift=2.817,
    outfile='my_dual_fit.h5',
    nlive_init=500
)

# Extract results
print("Region A best-fit parameters:", output['bestfit_params_a'])
print("Region B best-fit parameters:", output['bestfit_params_b'])
```

See `examples/exampl_usage.py`for additional info

## How It Works

### The Problem

Standard SED fitting tools fit one model to one set of observations. But what if you have:

- **Resolved bands**: Two regions with separate photometry (e.g., HST)
- **Blended bands**: Both regions combined in one measurement (e.g., WISE)

You need to fit both regions simultaneously with the constraint:

```
flux_model_A(WISE) + flux_model_B(WISE) = flux_observed(WISE)
```

### The Solution

This package implements a custom likelihood function:

```python
ln L_total = ln L_A(resolved) + ln L_B(resolved) + ln L(blended)
```

where:
- `ln L_A(resolved)`: Chi-square fit of Region A model to Region A resolved photometry
- `ln L_B(resolved)`: Chi-square fit of Region B model to Region B resolved photometry  
- `ln L(blended)`: Chi-square comparing observed blended flux to sum of both model predictions

The parameter space is `θ = [θ_A, θ_B]` where each region has independent stellar population parameters.

## Citation

**This work:**
```bibtex
@software{your_name_2024_prospector_dual,
  author       = {Aritra Aich},
  title        = {Prospector Dual-Region SED Fitting},
  month        = dec,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.17805754},
  url          = {https://doi.org/10.5281/zenodo.17805754}
}
```

**Prospector:**
```bibtex
@ARTICLE{2021ApJS..254...22J,
   author = {{Johnson}, Benjamin D. and {Leja}, Joel and 
             {Conroy}, Charlie and {Speagle}, Joshua S.},
   title = "{Stellar Population Inference with Prospector}",
   journal = {ApJS},
   year = 2021,
   volume = {254},
   pages = {22},
   doi = {10.3847/1538-4365/abef67}
}
```

## It currently works for


- **Lensed galaxies**: Multiple lensed images resolved in HST but blended in WISE/Spitzer
- **Spatially resolved galaxies**: Different regions (e.g., bulge + disk) with partial resolution

## Examples

Coming soon: Example notebooks for:
- [ ] Lensed galaxy with multiple images
- [ ] Star-forming clumps in a high-z galaxy
- [ ] Post-merger system with two nuclei

## Contributing

Contributions welcome! Some ideas:

- [ ] Extend to >2 regions
- [ ] Support for Prospector v2.X Observation classes

## Contact

- **Issues**: Please open an issue on GitHub

## Acknowledgments

This pkg builds on Prospector(https://github.com/bd-j/prospector)

