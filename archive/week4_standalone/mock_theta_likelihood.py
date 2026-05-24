"""
Mock likelihood for TEP-Cobaya testing
Computes chi2 based on theta* (100 * theta_s) compared to Planck target
"""

import numpy as np
from cobaya.likelihood import Likelihood

class mock_theta_likelihood(Likelihood):
    """
    Mock likelihood that compares computed theta* to Planck 2018 value.
    
    Planck 2018 theta* = 1.04109 ± 0.00030 (for 100*theta_s)
    """
    
    # Parameters that the likelihood takes
    params = {"theta_s_100": None}
    
    def initialize(self):
        """Set up the likelihood with target value and uncertainty"""
        # Planck 2018 TT,TE,EE+lowE best fit
        self.target_theta = self.target_theta_s_100  # 100 * theta_s
        self.sigma_theta = self.sigma_theta_s_100
        
        self.log.info(f"Mock theta* likelihood initialized")
        self.log.info(f"  Target: {self.target_theta:.5f} ± {self.sigma_theta:.5f}")
    
    def get_requirements(self):
        """What we need from the theory code"""
        return {"theta_s_100": None}
    
    def logp(self, **params_values):
        """
        Compute log-likelihood given derived parameters from theory.
        
        This is called by Cobaya with the derived parameters.
        """
        # Get theta* from the theory
        theta_s_100 = self.provider.get_param("theta_s_100")
        
        # Compute chi2
        delta = theta_s_100 - self.target_theta
        chi2 = (delta / self.sigma_theta) ** 2
        
        # Return log-likelihood (neglecting constant)
        loglike = -0.5 * chi2
        
        return loglike
    
    def get_can_provide_params(self):
        """Parameters this likelihood can provide"""
        return []
