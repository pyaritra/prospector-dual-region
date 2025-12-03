"""
Example usage of dual-region Prospector fitting with your actual data.

This shows how to adapt the dual_region_prospector.py script to your
specific HST + WISE dataset.
"""

import numpy as np
from prospector_dual_region import run_dual_region_fit
import matplotlib.pyplot as plt


def load_your_data():
    """
    Replace this with your actual data loading code.
    
    Returns
    -------
    Dictionary with all photometric data for both regions
    """
    
    # YOUR CODE HERE: Load your HST and WISE photometry
    # This is just a template showing the expected format
    
    data = {
        # Region A
        'filters_a': ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 
                      'wfc3_uvis_f606w', 'wfc3_uvis_f814w'],
        'maggies_a': np.array([...]),  # Your HST photometry in maggies
        'maggies_unc_a': np.array([...]),  # Uncertainties
        
        # Region B
        'filters_b': ['wfc3_uvis_f275w', 'wfc3_uvis_f336w',
                      'wfc3_uvis_f606w', 'wfc3_uvis_f814w'],
        'maggies_b': np.array([...]),
        'maggies_unc_b': np.array([...]),
        
        # WISE (blended)
        'filters_wise': ['wise_w1', 'wise_w2', 'wise_w3', 'wise_w4'],
        'maggies_wise': np.array([...]),  # Total WISE flux
        'maggies_unc_wise': np.array([...]),
        
        # Galaxy info
        'redshift': 2.817  # Your galaxy redshift
    }
    
    return data


def magnitudes_to_maggies(mags, mag_unc):
    """
    Convert AB magnitudes to maggies (Jy/3631).
    
    Parameters
    ----------
    mags : array
        AB magnitudes
    mag_unc : array
        Magnitude uncertainties
        
    Returns
    -------
    maggies : array
        Fluxes in maggies
    maggies_unc : array
        Uncertainties in maggies
    """
    maggies = 10**(-0.4 * mags)
    # Uncertainty propagation
    maggies_unc = maggies * mag_unc * np.log(10) / 2.5
    return maggies, maggies_unc


def run_your_fit():
    """
    Main function to run the dual-region fit on your data.
    """
    
    # Load data
    data = load_your_data()
    
    # Run the fit
    print("Running dual-region SED fit...")
    output = run_dual_region_fit(
        # Region A data
        hst_filters_a=data['filters_a'],
        hst_phot_a=data['maggies_a'],
        hst_unc_a=data['maggies_unc_a'],
        
        # Region B data
        hst_filters_b=data['filters_b'],
        hst_phot_b=data['maggies_b'],
        hst_unc_b=data['maggies_unc_b'],
        
        # WISE data (blended)
        wise_filters=data['filters_wise'],
        wise_phot=data['maggies_wise'],
        wise_unc=data['maggies_unc_wise'],
        
        # Galaxy properties
        redshift=data['redshift'],
        
        # Output
        outfile='my_galaxy_dual_fit.h5',
        
        # Sampling parameters (adjust as needed)
        nlive_init=500,  # More live points = better sampling but slower
        dlogz_init=0.01,  # Stopping criterion (smaller = more accurate)
    )
    
    return output


def plot_results(output):
    """
    Make diagnostic plots of the fit results.
    """
    from prospect.plotting import pretty
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Get samples
    results = output['sampling']
    weights = np.exp(results['logwt'] - results['logz'][-1])
    samples = results['samples']
    
    nfree_a = output['nfree_a']
    
    # Plot corner plot for Region A parameters (first few)
    # You would typically use corner.py for this
    ax = axes[0, 0]
    ax.set_title('Region A: Mass vs Age')
    # Extract mass and age parameters for region A
    # (adjust indices based on your model)
    # ax.scatter(samples[:, 0], samples[:, 1], c=weights, alpha=0.3)
    
    # Similar for Region B
    ax = axes[0, 1]
    ax.set_title('Region B: Mass vs Age')
    
    # Plot photometric fits
    ax = axes[1, 0]
    ax.set_title('Region A Photometry')
    # Plot observed vs model
    
    ax = axes[1, 1]
    ax.set_title('WISE Blended Photometry')
    # Show that sum of models matches observed WISE
    
    plt.tight_layout()
    plt.savefig('dual_region_fit_results.png', dpi=150)
    print("Saved results plot to dual_region_fit_results.png")


def extract_physical_properties(output):
    """
    Extract physical properties (masses, SFRs, etc.) for both regions.
    """
    # Get best-fit parameters
    theta_a = output['bestfit_params_a']
    theta_b = output['bestfit_params_b']
    
    model_a = output['model_a']
    model_b = output['model_b']
    
    # Set parameters
    model_a.set_parameters(theta_a)
    model_b.set_parameters(theta_b)
    
    # Extract physical properties
    # (These depend on your model parameterization)
    props_a = {
        'log_mass': theta_a[model_a.free_params.index('logmass')],
        'dust2': theta_a[model_a.free_params.index('dust2')],
        # Add other parameters
    }
    
    props_b = {
        'log_mass': theta_b[model_b.free_params.index('logmass')],
        'dust2': theta_b[model_b.free_params.index('dust2')],
        # Add other parameters
    }
    
    print("\n=== Region A Properties ===")
    for key, val in props_a.items():
        print(f"{key}: {val:.3f}")
    
    print("\n=== Region B Properties ===")
    for key, val in props_b.items():
        print(f"{key}: {val:.3f}")
    
    return props_a, props_b


if __name__ == "__main__":
    
    # Run the fit
    output = run_your_fit()
    
    # Extract physical properties
    props_a, props_b = extract_physical_properties(output)
    
    # Make plots
    plot_results(output)
    
    print("\nFit complete! Check output file and plots.")


# TIPS FOR CUSTOMIZATION:
# 
# 1. FILTER NAMES: Use sedpy filter names. Common ones:
#    HST: 'wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w'
#    WISE: 'wise_w1', 'wise_w2', 'wise_w3', 'wise_w4'
#    Full list: from sedpy import observate; observate.list_filters()
#
# 2. PRIORS: Customize priors in build_dual_models() function in dual_region_prospector.py
#    Example: Set different mass ranges for each region based on your expectations
#
# 3. SFH: The default is delayed-tau. You can change to other SFHs by modifying
#    the TemplateLibrary call in build_dual_models()
#
# 4. NEBULAR EMISSION: Add nebular emission by including 'add_neb': True in model_params
#
# 5. SAMPLING: Increase nlive_init for better sampling (but slower).
#    Decrease dlogz_init for more accurate evidence estimates.
