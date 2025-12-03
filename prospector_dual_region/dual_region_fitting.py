"""
Custom Prospector implementation for fitting two spatially resolved regions
with shared (blended) WISE photometry constraint.

This implements simultaneous SED fitting where:
- HST bands: resolved photometry for each region separately
- WISE bands: blended photometry with constraint that model_A + model_B = observed
"""

import numpy as np
from prospect.models import SpecModel
from prospect.sources import CSPSpecBasis
from prospect import prospect_args
from prospect.fitting import fit_model
from prospect.io import write_results as writer
import sedpy


class DualRegionLikelihood:
    """
    Custom likelihood function for two regions with partial blending.
    
    Parameters
    ----------
    obs_region_a : list of Observation objects
        Observations for region A (resolved HST only)
    obs_region_b : list of Observation objects  
        Observations for region B (resolved HST only)
    obs_wise : list of Observation objects
        WISE observations (blended, single Observation)
    model_a : SpecModel
        Model for region A
    model_b : SpecModel
        Model for region B
    sps_a : CSPSpecBasis
        SPS object for region A
    sps_b : CSPSpecBasis
        SPS object for region B
    """
    
    def __init__(self, obs_region_a, obs_region_b, obs_wise,
                 model_a, model_b, sps_a, sps_b):
        self.obs_region_a = obs_region_a
        self.obs_region_b = obs_region_b
        self.obs_wise = obs_wise
        self.model_a = model_a
        self.model_b = model_b
        self.sps_a = sps_a
        self.sps_b = sps_b
        
        # Total number of parameters
        self.ndim = model_a.ndim + model_b.ndim
        
    def __call__(self, theta, nested=False, verbose=False):
        """
        Compute log posterior probability for combined system.
        
        Parameters
        ----------
        theta : array-like
            Combined parameter vector [params_A, params_B]
        nested : bool
            If True, return only log-likelihood (for nested sampling)
        verbose : bool
            Print debug info
            
        Returns
        -------
        lnprob : float
            Log posterior probability (or log likelihood if nested=True)
        """
        # Split parameters
        theta_a = theta[:self.model_a.ndim]
        theta_b = theta[self.model_a.ndim:]
        
        # Update models with new parameters
        self.model_a.set_parameters(theta_a)
        self.model_b.set_parameters(theta_b)
        
        # Check priors
        lnp_prior_a = self.model_a.prior_product(theta_a, nested=nested)
        lnp_prior_b = self.model_b.prior_product(theta_b, nested=nested)
        
        if (not np.isfinite(lnp_prior_a)) or (not np.isfinite(lnp_prior_b)):
            return -np.inf
        
        # Get model predictions for both regions
        spec_a, phot_a, extras_a = self.model_a.predict(theta_a, obs=self.obs_region_a, sps=self.sps_a)
        spec_b, phot_b, extras_b = self.model_b.predict(theta_b, obs=self.obs_region_b, sps=self.sps_b)
        
        # Also predict WISE photometry for both regions
        _, phot_wise_a, _ = self.model_a.predict(theta_a, obs=self.obs_wise, sps=self.sps_a)
        _, phot_wise_b, _ = self.model_b.predict(theta_b, obs=self.obs_wise, sps=self.sps_b)
        
        # Compute likelihoods
        # 1. Region A HST photometry
        lnlike_a = self._compute_phot_likelihood(
            self.obs_region_a[0], phot_a
        )
        
        # 2. Region B HST photometry  
        lnlike_b = self._compute_phot_likelihood(
            self.obs_region_b[0], phot_b
        )
        
        # 3. WISE blended photometry (sum of models)
        phot_wise_total = phot_wise_a + phot_wise_b
        lnlike_wise = self._compute_phot_likelihood(
            self.obs_wise[0], phot_wise_total
        )
        
        # Total likelihood
        lnlike = lnlike_a + lnlike_b + lnlike_wise
        
        if verbose:
            print(f"lnlike_A: {lnlike_a:.2f}, lnlike_B: {lnlike_b:.2f}, "
                  f"lnlike_WISE: {lnlike_wise:.2f}")
        
        # Return likelihood or posterior
        if nested:
            return lnlike
        else:
            return lnlike + lnp_prior_a + lnp_prior_b
    
    def _compute_phot_likelihood(self, obs, model_phot):
        """
        Compute Gaussian log-likelihood for photometry.
        
        Parameters
        ----------
        obs : Observation object
            Must have 'maggies', 'maggies_unc', and 'phot_mask' attributes
        model_phot : array
            Model photometry predictions
            
        Returns
        -------
        lnlike : float
            Log likelihood
        """
        # Apply mask (True = use this data point)
        mask = obs['phot_mask']
        
        # Get data and uncertainties
        data = obs['maggies'][mask]
        unc = obs['maggies_unc'][mask]
        model = model_phot[mask]
        
        # Chi-square
        chi2 = np.sum(((data - model) / unc)**2)
        
        # Log likelihood (Gaussian)
        lnlike = -0.5 * (chi2 + np.sum(np.log(2 * np.pi * unc**2)))
        
        return lnlike


def build_observations_dual_region(hst_filters_a, hst_phot_a, hst_unc_a,
                                    hst_filters_b, hst_phot_b, hst_unc_b,
                                    wise_filters, wise_phot, wise_unc):
    """
    Build Observation objects for dual region fitting.
    
    Parameters
    ----------
    hst_filters_a : list of str
        Filter names for region A (e.g., ['f606w', 'f814w'])
    hst_phot_a : array
        HST photometry for region A in maggies
    hst_unc_a : array
        HST uncertainties for region A in maggies
    (similar for region B)
    wise_filters : list of str
        WISE filter names (e.g., ['wise_w1', 'wise_w2'])
    wise_phot : array
        WISE photometry (blended, total) in maggies
    wise_unc : array
        WISE uncertainties in maggies
        
    Returns
    -------
    obs_a : list
        Observation list for region A
    obs_b : list  
        Observation list for region B
    obs_wise : list
        Observation list for WISE
    """
    from prospect.observation import Photometry
    
    # Region A HST
    obs_a_dict = {
        'filters': sedpy.observate.load_filters(hst_filters_a),
        'maggies': hst_phot_a,
        'maggies_unc': hst_unc_a,
        'phot_mask': np.ones(len(hst_phot_a), dtype=bool)
    }
    obs_a = [Photometry(**obs_a_dict)]
    
    # Region B HST
    obs_b_dict = {
        'filters': sedpy.observate.load_filters(hst_filters_b),
        'maggies': hst_phot_b,
        'maggies_unc': hst_unc_b,
        'phot_mask': np.ones(len(hst_phot_b), dtype=bool)
    }
    obs_b = [Photometry(**obs_b_dict)]
    
    # WISE (blended)
    obs_wise_dict = {
        'filters': sedpy.observate.load_filters(wise_filters),
        'maggies': wise_phot,
        'maggies_unc': wise_unc,
        'phot_mask': np.ones(len(wise_phot), dtype=bool)
    }
    obs_wise = [Photometry(**obs_wise_dict)]
    
    return obs_a, obs_b, obs_wise


def build_dual_models(redshift=None, **extras):
    """
    Build two independent stellar population models.
    
    Parameters
    ----------
    redshift : float
        Galaxy redshift (fixed)
    **extras : dict
        Additional configuration parameters
        
    Returns
    -------
    model_a : SpecModel
        Model for region A
    model_b : SpecModel
        Model for region B
    """
    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors
    
    # Get default parameter set (delayed-tau SFH)
    model_params_a = TemplateLibrary["parametric_sfh"]
    model_params_b = TemplateLibrary["parametric_sfh"].copy()
    
    # Fix redshift if provided
    if redshift is not None:
        model_params_a["zred"]["init"] = redshift
        model_params_a["zred"]["isfree"] = False
        model_params_b["zred"]["init"] = redshift
        model_params_b["zred"]["isfree"] = False
    
    # You can customize priors here
    # For example, different mass ranges for each region:
    # model_params_a["logmass"]["prior"] = priors.TopHat(mini=8.0, maxi=11.0)
    # model_params_b["logmass"]["prior"] = priors.TopHat(mini=7.0, maxi=10.0)
    
    model_a = SpecModel(model_params_a)
    model_b = SpecModel(model_params_b)
    
    return model_a, model_b


def run_dual_region_fit(hst_filters_a, hst_phot_a, hst_unc_a,
                        hst_filters_b, hst_phot_b, hst_unc_b,
                        wise_filters, wise_phot, wise_unc,
                        redshift=None, outfile='dual_region_fit.h5',
                        **fitting_kwargs):
    """
    Run the dual-region SED fit with blended WISE constraint.
    
    Parameters
    ----------
    [photometry parameters as in build_observations_dual_region]
    redshift : float
        Galaxy redshift
    outfile : str
        Output HDF5 file name
    **fitting_kwargs : dict
        Additional arguments for dynesty sampler
        
    Returns
    -------
    output : dict
        Fitting results
    """
    # Build observations
    obs_a, obs_b, obs_wise = build_observations_dual_region(
        hst_filters_a, hst_phot_a, hst_unc_a,
        hst_filters_b, hst_phot_b, hst_unc_b,
        wise_filters, wise_phot, wise_unc
    )
    
    # Build models
    model_a, model_b = build_dual_models(redshift=redshift)
    
    # Build SPS objects
    sps_a = CSPSpecBasis(zcontinuous=1)
    sps_b = CSPSpecBasis(zcontinuous=1)
    
    # Create custom likelihood
    dual_lnprob = DualRegionLikelihood(
        obs_a, obs_b, obs_wise,
        model_a, model_b,
        sps_a, sps_b
    )
    
    # Set up initial parameters
    # Concatenate initial guesses from both models
    initial_theta = np.concatenate([model_a.theta, model_b.theta])
    
    # Run fitting with dynesty
    # Note: We need to wrap the likelihood for dynesty
    import dynesty
    from dynesty import utils as dyfunc
    
    # Define prior transform for dynesty
    def prior_transform(u):
        """Transform unit cube to parameter space"""
        theta = np.zeros(dual_lnprob.ndim)
        
        # Region A parameters
        for i in range(model_a.ndim):
            pname = model_a.free_params[i]
            prior = model_a.config_dict[pname]['prior']
            theta[i] = prior.unit_transform(u[i])
        
        # Region B parameters  
        for i in range(model_b.ndim):
            j = i + model_a.ndim
            pname = model_b.free_params[i]
            prior = model_b.config_dict[pname]['prior']
            theta[j] = prior.unit_transform(u[j])
            
        return theta
    
    # Define likelihood for dynesty (just log-likelihood)
    def loglike(theta):
        return dual_lnprob(theta, nested=True)
    
    # Set default dynesty parameters
    dyn_kwargs = {
        'nlive_init': 500,
        'sample': 'rwalk',
        'bound': 'multi',
        'dlogz_init': 0.01,
    }
    dyn_kwargs.update(fitting_kwargs)
    
    print("Starting dual-region fit with dynesty...")
    print(f"Total parameters: {dual_lnprob.ndim}")
    print(f"Region A free parameters: {model_a.ndim}")
    print(f"Region B free parameters: {model_b.ndim}")
    
    sampler = dynesty.NestedSampler(
        loglike, prior_transform, dual_lnprob.ndim,
        **{k: v for k, v in dyn_kwargs.items() if k != 'dlogz_init'}
    )
    
    sampler.run_nested(dlogz=dyn_kwargs['dlogz_init'])
    results = sampler.results
    
    # Process results
    output = {
        'sampling': results,
        'obs_a': obs_a,
        'obs_b': obs_b, 
        'obs_wise': obs_wise,
        'model_a': model_a,
        'model_b': model_b,
        'ndim': dual_lnprob.ndim,
        'nfree_a': model_a.ndim,
        'nfree_b': model_b.ndim
    }
    
    # Get best-fit parameters
    weights = np.exp(results['logwt'] - results['logz'][-1])
    best_idx = np.argmax(results['logwt'])
    output['bestfit_params'] = results['samples'][best_idx]
    output['bestfit_params_a'] = output['bestfit_params'][:model_a.ndim]
    output['bestfit_params_b'] = output['bestfit_params'][model_a.ndim:]
    
    print("Fit complete!")
    print(f"Log evidence: {results['logz'][-1]:.2f}")
    
    return output


# Example usage
if __name__ == "__main__":
    
    # Example: Mock data
    # HST photometry for region A (F606W, F814W)
    hst_filters_a = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w']
    hst_phot_a = np.array([1.5e-6, 2.3e-6])  # maggies
    hst_unc_a = np.array([1.5e-7, 2.3e-7])
    
    # HST photometry for region B  
    hst_filters_b = ['wfc3_uvis_f606w', 'wfc3_uvis_f814w']
    hst_phot_b = np.array([0.8e-6, 1.1e-6])
    hst_unc_b = np.array([1.0e-7, 1.5e-7])
    
    # WISE photometry (blended)
    wise_filters = ['wise_w1', 'wise_w2']
    wise_phot = np.array([5.0e-6, 4.5e-6])  # Should be ~ sum of A+B in mid-IR
    wise_unc = np.array([5.0e-7, 5.0e-7])
    
    # Run fit
    output = run_dual_region_fit(
        hst_filters_a, hst_phot_a, hst_unc_a,
        hst_filters_b, hst_phot_b, hst_unc_b,
        wise_filters, wise_phot, wise_unc,
        redshift=0.8,
        outfile='test_dual_fit.h5',
        nlive_init=400,  # Smaller for testing
        dlogz_init=0.05
    )
    
    print("\nBest-fit parameters for Region A:")
    print(output['bestfit_params_a'])
    print("\nBest-fit parameters for Region B:")
    print(output['bestfit_params_b'])
