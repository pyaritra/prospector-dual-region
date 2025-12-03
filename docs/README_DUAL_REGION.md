# Dual-Region SED Fitting with Prospector

Custom implementation for fitting two spatially resolved regions with shared (blended) mid-IR photometry.

## Problem Statement

You have:
- **HST photometry**: Resolved measurements for Region A and Region B separately
- **WISE photometry**: Blended measurements where both regions contribute
- **Goal**: Fit stellar population models to both regions simultaneously with the constraint that their combined WISE flux matches the observed blended flux

## Solution

This implementation creates a custom likelihood function that:

1. Fits Region A model to Region A HST data
2. Fits Region B model to Region B HST data  
3. Enforces that `Model_A(WISE) + Model_B(WISE) = Observed(WISE)`

## Files

- **dual_region_prospector.py**: Main implementation with custom likelihood class
- **example_usage.py**: Template showing how to use it with your data

## Key Components

### DualRegionLikelihood Class

Custom likelihood function that handles:
- Two independent stellar population models (one per region)
- Two independent SPS (stellar population synthesis) objects
- Separate observation sets for resolved and blended data
- Combined likelihood calculation with blending constraint

### Installation Requirements

```bash
pip install astro-prospector
pip install fsps python-fsps
pip install sedpy
pip install dynesty
```

You'll also need FSPS stellar population models installed. See:
https://prospect.readthedocs.io/en/latest/installation.html

## Usage

### 1. Prepare Your Data

Convert your photometry to "maggies" (flux units of Jy/3631):

```python
maggies = 10**(-0.4 * AB_magnitudes)
```

Organize into the format expected by the functions:

```python
# HST data for Region A
hst_filters_a = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w', ...]
hst_maggies_a = np.array([...])
hst_unc_a = np.array([...])

# HST data for Region B
hst_filters_b = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w', ...]
hst_maggies_b = np.array([...])
hst_unc_b = np.array([...])

# WISE data (blended - total from both regions)
wise_filters = ['wise_w1', 'wise_w2', 'wise_w3', 'wise_w4']
wise_maggies = np.array([...])  # Total observed flux
wise_unc = np.array([...])
```

### 2. Run the Fit

```python
from dual_region_prospector import run_dual_region_fit

output = run_dual_region_fit(
    hst_filters_a, hst_maggies_a, hst_unc_a,
    hst_filters_b, hst_maggies_b, hst_unc_b,
    wise_filters, wise_maggies, wise_unc,
    redshift=2.817,  # Your galaxy redshift
    outfile='my_fit_results.h5',
    nlive_init=500,
    dlogz_init=0.01
)
```

### 3. Extract Results

```python
# Best-fit parameters for each region
params_a = output['bestfit_params_a']
params_b = output['bestfit_params_b']

# Access the sampling results
results = output['sampling']
weights = np.exp(results['logwt'] - results['logz'][-1])
samples = results['samples']

# Get physical properties
model_a = output['model_a']
model_a.set_parameters(params_a)
# Access model properties like mass, SFR, etc.
```

## How It Works

### Parameter Space

The combined parameter vector is:
```
theta = [params_region_A, params_region_B]
```

For a standard delayed-tau SFH, each region has ~5-7 free parameters:
- log(M*): stellar mass
- log(τ): SFH timescale
- log(age): stellar age  
- dust2: dust attenuation
- log(Z): metallicity
- ...

Total: ~10-14 parameters

### Likelihood Calculation

```
ln P(data | theta) = ln P(HST_A | theta_A) + ln P(HST_B | theta_B) 
                     + ln P(WISE | theta_A + theta_B)
```

Where:
- `ln P(HST_A | theta_A)`: Normal chi-square for Region A HST bands
- `ln P(HST_B | theta_B)`: Normal chi-square for Region B HST bands
- `ln P(WISE | theta_A + theta_B)`: Chi-square comparing observed WISE to sum of model predictions

### Priors

Each region has independent priors on its parameters. You can customize these in the `build_dual_models()` function:

```python
# Example: Different mass ranges for each region
from prospect.models import priors

model_params_a["logmass"]["prior"] = priors.TopHat(mini=9.0, maxi=11.0)
model_params_b["logmass"]["prior"] = priors.TopHat(mini=8.0, maxi=10.0)
```

## Customization Options

### 1. Change Star Formation History

Modify `build_dual_models()` to use different SFH templates:

```python
from prospect.models.templates import TemplateLibrary

# Options: "parametric_sfh", "continuity_sfh", "nonparametric_sfh"
model_params_a = TemplateLibrary["continuity_sfh"]
```

### 2. Add Nebular Emission

```python
model_params_a["add_neb"] = {"N": 1, "isfree": False, "init": True}
```

### 3. Add AGN Component

If one region has AGN, you can use `AGNSpecModel` instead of `SpecModel`.

### 4. Adjust Sampling Parameters

```python
output = run_dual_region_fit(
    ...,
    nlive_init=1000,      # More live points = better sampling (slower)
    dlogz_init=0.005,     # Smaller = more accurate stopping
    sample='rwalk',       # Sampling algorithm
    bound='multi',        # Bounding method
)
```

## Performance Notes

- **Typical runtime**: 2-12 hours depending on parameter space and nlive_init
- Each likelihood call requires:
  - 2 SPS predictions (one per region)
  - 2 additional WISE predictions
  - ~3x slower than single-region fits

- **Memory**: ~1-2 GB for typical parameter space

## Limitations

1. **Assumes independent stellar populations**: The two regions are treated as completely separate stellar systems. If they have correlated star formation histories, this might not be appropriate.

2. **No spatial information**: The model doesn't use any spatial information beyond the resolved/unresolved designation of bands.

3. **Linear blending**: Assumes simple flux addition for WISE bands. Doesn't account for complex PSF effects or scattered light.

4. **Computational cost**: ~3x slower than single-region fits due to doubled model evaluations.

## Troubleshooting

### Issue: Fit is very slow
- Reduce `nlive_init` (try 300-400 for testing)
- Increase `dlogz_init` (try 0.05 for testing)
- Check that priors aren't too broad

### Issue: WISE constraint not satisfied
- Check that WISE fluxes are in maggies, not AB mags
- Verify WISE flux is actually the sum of both regions
- Check filter names match sedpy conventions

### Issue: Poor convergence
- Increase `nlive_init` (try 800-1000)
- Tighten priors based on preliminary fits
- Check for bad data points (set mask to False)

### Issue: Import errors
- Make sure FSPS is properly installed: `python -c "import fsps"`
- Check Prospector version: `python -c "import prospect; print(prospect.__version__)"`

## References

- Johnson et al. 2021, ApJS, 254, 22 (Prospector paper)
- Leja et al. 2017, ApJ, 837, 170 (SED fitting methodology)
- Conroy et al. 2009, ApJ, 699, 486 (FSPS)

## Contact

For issues specific to this dual-region implementation, check the comments in the code.
For general Prospector questions, see: https://github.com/bd-j/prospector/issues
