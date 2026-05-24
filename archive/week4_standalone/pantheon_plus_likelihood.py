"""
Pantheon+ SNe Ia likelihood for Cobaya
Implements TEP distance modulus calculation
"""

import numpy as np
import os
from cobaya.likelihood import Likelihood


class pantheon_plus(Likelihood):
    """
    Pantheon+ SNe Ia likelihood with TEP cosmology support.
    
    Uses Pantheon+ compilation of 1701 SNe Ia covering z = 0.001 to 2.26.
    Includes SH0ES calibration for absolute magnitude.
    
    Data files expected in:
    - data/raw/pantheon_plus_shoes.dat
    - data/raw/pantheon_plus_shoes.cov
    """
    
    # No required parameters - we compute distances from theory
    # Theory must provide angular diameter distance
    
    def initialize(self):
        """Load Pantheon+ data and covariance matrix."""
        
        # Find data files
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, "data", "raw")
        
        data_file = os.path.join(data_dir, "pantheon_plus_shoes.dat")
        cov_file = os.path.join(data_dir, "pantheon_plus_shoes.cov")
        
        if not os.path.exists(data_file):
            raise FileNotFoundError(f"Pantheon+ data file not found: {data_file}")
        if not os.path.exists(cov_file):
            raise FileNotFoundError(f"Pantheon+ covariance file not found: {cov_file}")
        
        self.log.info(f"Loading Pantheon+ data from {data_file}")
        
        # Load SNe data
        # Pantheon+ format (columns from header):
        # Col 3: zHD (redshift) 
        # Col 9: m_b_corr (corrected apparent magnitude)
        # Col 10: m_b_corr_err_DIAG (error)
        # Col 11: MU_SH0ES (distance modulus from SH0ES)
        # Col 12: MU_SH0ES_ERR_DIAG (distance modulus error)
        self.z_sn, self.mu_sn, self.dmu_sn = [], [], []
        
        with open(data_file, 'r') as f:
            lines = f.readlines()
            
        # Skip header (first line)
        for line in lines[1:]:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                self.z_sn.append(float(parts[2]))      # Col 3: zHD
                # Use distance modulus directly from SH0ES calibration
                self.mu_sn.append(float(parts[10]))    # Col 11: MU_SH0ES
                self.dmu_sn.append(float(parts[11]))   # Col 12: MU_SH0ES_ERR_DIAG
        
        self.z_sn = np.array(self.z_sn)
        self.mu_sn = np.array(self.mu_sn)
        self.dmu_sn = np.array(self.dmu_sn)
        self.n_sn = len(self.z_sn)
        
        self.log.info(f"Loaded {self.n_sn} SNe Ia")
        self.log.info(f"Redshift range: [{self.z_sn.min():.4f}, {self.z_sn.max():.4f}]")
        self.log.info(f"Distance modulus range: [{self.mu_sn.min():.2f}, {self.mu_sn.max():.2f}]")
        
        # Load covariance matrix
        # Pantheon+ format: first line is N, then N^2 values follow
        self.log.info(f"Loading covariance from {cov_file}")
        with open(cov_file, 'r') as f:
            lines = f.readlines()
        
        # First line is number of SNe (should match)
        n_cov = int(lines[0].strip())
        assert n_cov == self.n_sn, f"Covariance size {n_cov} doesn't match data size {self.n_sn}"
        
        # Read remaining lines as floats
        cov_elements = []
        for line in lines[1:]:
            if line.strip():
                cov_elements.append(float(line.strip()))
        
        self.cov = np.array(cov_elements).reshape(self.n_sn, self.n_sn)
        
        # Add diagonal statistical uncertainties
        self.cov_stat_sys = self.cov.copy()
        for i in range(self.n_sn):
            self.cov_stat_sys[i, i] += self.dmu_sn[i]**2
        
        # Compute inverse covariance
        self.inv_cov = np.linalg.inv(self.cov_stat_sys)
        
        self.log.info(f"Covariance matrix shape: {self.cov.shape}")
        
    def get_requirements(self):
        """What we need from the theory code."""
        # We need angular diameter distances at SNe redshifts
        return {
            "angular_diameter_distance": {"z": list(self.z_sn)}
        }
    
    def _compute_tep_distance_modulus(self, z, Da):
        """
        Compute distance modulus for TEP cosmology.
        
        For TEP: D_L = D_A * (1 + z_TEP)^2
        where z_TEP includes temporal transport contribution.
        
        In our CLASS implementation, the angular diameter distance
        from CLASS already includes TEP effects through H_TEP(z).
        
        For standard interpretation:
        mu = 5 * log10(D_L / 10pc) = 5 * log10(D_A * (1+z)^2) - 5
        
        With TEP, we use the distance from CLASS which has H_TEP built in.
        """
        # D_L in Mpc
        D_L = Da * (1 + z)**2
        
        # Distance modulus: mu = 5*log10(D_L) + 25 (when D_L in Mpc)
        mu = 5.0 * np.log10(D_L) + 25.0
        
        return mu
    
    def logp(self, **params_values):
        """
        Compute log-likelihood for Pantheon+ SNe.
        
        This is called by Cobaya with parameters from theory.
        """
        # Get angular diameter distances from theory
        Da_list = []
        for z in self.z_sn:
            Da = self.provider.get_angular_diameter_distance(z)
            # Ensure Da is a scalar (not an array)
            if hasattr(Da, '__len__'):
                Da = float(Da[0]) if len(Da) > 0 else float(Da)
            else:
                Da = float(Da)
            Da_list.append(Da)
        
        Da = np.array(Da_list)
        
        # Compute theoretical distance moduli
        mu_theory = self._compute_tep_distance_modulus(self.z_sn, Da)
        
        # Compute residuals: data - theory
        residuals = self.mu_sn - mu_theory
        
        # Chi2 calculation: chi2 = res^T C^{-1} res
        chi2 = np.dot(residuals, np.dot(self.inv_cov, residuals))
        
        # Log-likelihood (neglecting constant)
        loglike = -0.5 * chi2
        
        return loglike
    
    def get_can_provide_params(self):
        """Parameters this likelihood can provide."""
        return ["chi2_pantheon"]


# For testing without Cobaya
if __name__ == "__main__":
    print("Testing Pantheon+ likelihood...")
    
    # Create instance
    like = pantheon_plus({})
    like.initialize()
    
    print(f"\nTest complete!")
    print(f"Number of SNe: {like.n_sn}")
    print(f"Redshift range: [{like.z_sn.min():.4f}, {like.z_sn.max():.4f}]")
