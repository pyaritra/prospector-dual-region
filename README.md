# Prospector Dual-Region SED Fitting

A custom extension for [Prospector](https://github.com/bd-j/prospector) that enables simultaneous SED fitting of two spatially resolved regions with shared (blended) photometry in some bands.

**Use case**: When you have high-resolution photometry (e.g., HST) that resolves two regions separately, but low-resolution photometry (e.g., WISE) where both regions are blended together.

## Features

- ✅ Simultaneous fitting of two independent stellar population models
- ✅ Enforces constraint that blended bands equal sum of model predictions
- ✅ Full Bayesian inference with nested sampling (dynesty)
- ✅ Compatible with all standard Prospector features (flexible SFHs, nebular emission, dust, etc.)
- ✅ Returns posteriors for both regions' physical properties

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

See `examples/example_usage.py` for a more detailed example.

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

## Documentation

- **[docs/README_DUAL_REGION.md](docs/README_DUAL_REGION.md)**: Detailed documentation covering:
  - Data preparation and format requirements
  - Customization options (priors, SFH types, nebular emission)
  - Performance optimization
  - Troubleshooting guide

- **[examples/example_usage.py](examples/example_usage.py)**: Complete example showing:
  - How to load and format your data
  - How to run the fit
  - How to extract physical properties
  - How to make diagnostic plots

## Citation

If you use this code in your research, please cite both this repository and the original Prospector paper:

**This work:**
```
[Your citation format - could add a Zenodo DOI]
```

**Prospector:**
```bibtex
@ARTICLE{2021ApJS..254...22J,
   author = {{Johnson}, Benjamin D. and {Leja}, Joel and {Conroy}, Charlie and {Speagle}, Joshua S.},
   title = "{Stellar Population Inference with Prospector}",
   journal = {ApJS},
   year = 2021,
   volume = {254},
   pages = {22},
   doi = {10.3847/1538-4365/abef67}
}
```

## Use Cases

This tool is designed for scenarios like:

- **Lensed galaxies**: Multiple lensed images resolved in HST but blended in WISE/Spitzer
- **Merging galaxies**: Two distinct components in optical but blended in IR
- **Spatially resolved galaxies**: Different regions (e.g., bulge + disk) with partial resolution
- **Binary AGN/QSO**: Resolved optically but blended in IR

## Examples

Coming soon: Example notebooks for:
- [ ] Lensed galaxy with multiple images
- [ ] Star-forming clumps in a high-z galaxy
- [ ] Post-merger system with two nuclei

## Contributing

Contributions welcome! Some ideas:

- [ ] Extend to >2 regions
- [ ] Add support for spectroscopy
- [ ] Implement joint priors between regions
- [ ] Add convenience plotting functions
- [ ] Support for Prospector v2.X Observation classes

## Performance

**Typical runtimes** (on modern CPU):
- Test fit (`nlive_init=300`): ~1-2 hours
- Production fit (`nlive_init=500`): ~4-8 hours
- High-accuracy fit (`nlive_init=1000`): ~12-24 hours

Runtime scales approximately as:
- Linear with `nlive_init`
- ~3× slower than single-region fits (due to doubled SPS calls)

## License

MIT License (same as Prospector)

## Contact

- **Issues**: Please open an issue on GitHub
- **Questions**: [Your contact or discussion forum]

## Acknowledgments

This package builds on [Prospector](https://github.com/bd-j/prospector) by Ben Johnson, Joel Leja, Charlie Conroy, and Joshua Speagle.

Developed for spatially resolved SED fitting of gravitationally lensed galaxies and star-forming regions.
