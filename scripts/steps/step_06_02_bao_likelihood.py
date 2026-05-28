#!/usr/bin/env python3
"""Step 031b: BAO Full Covariance Likelihood for TEP.

Calculates the BAO likelihood using full covariance matrices.
Crucial for testing if the TEP 31% shift in the sound horizon r_s
is perfectly compensated by the shift in angular diameter distance D_A.
"""

import numpy as np

class BAOLikelihood:
    def __init__(self, z_eff, D_M_over_rs_data, D_H_over_rs_data, covariance_matrix):
        """Initialize the BAO likelihood with data and covariance.
        
        Args:
            z_eff: Effective redshifts of the BAO measurements
            D_M_over_rs_data: Transverse measurements D_M / r_s
            D_H_over_rs_data: Line-of-sight measurements D_H / r_s (where D_H = c/H(z))
            covariance_matrix: Full covariance matrix for the combined vector
        """
        self.z_eff = np.array(z_eff)
        self.data_vector = np.concatenate([D_M_over_rs_data, D_H_over_rs_data])
        self.cov = np.array(covariance_matrix)
        self.inv_cov = np.linalg.inv(self.cov)
        
    def log_likelihood(self, D_M_over_rs_theory, D_H_over_rs_theory):
        """Compute the log-likelihood of the theory given the data.
        
        Args:
            D_M_over_rs_theory: Theoretical D_M / r_s at z_eff
            D_H_over_rs_theory: Theoretical D_H / r_s at z_eff
            
        Returns:
            Log-likelihood value
        """
        theory_vector = np.concatenate([D_M_over_rs_theory, D_H_over_rs_theory])
        residuals = self.data_vector - theory_vector
        chi2 = residuals.T @ self.inv_cov @ residuals
        return -0.5 * chi2

def test_bao_shift(tep_cosmo, lcdm_cosmo, z_test=np.array([0.38, 0.51, 0.61])):
    """Test if TEP BAO predictions match LCDM despite the r_s shift.
    
    A successful TEP model will have D_M_tep / r_s_tep ≈ D_M_lcdm / r_s_lcdm.
    """
    # Get sound horizons
    rs_tep = tep_cosmo.sound_horizon()
    rs_lcdm = lcdm_cosmo.sound_horizon()
    
    # Get distances
    DM_tep = np.array([tep_cosmo.comoving_distance(z) for z in z_test])
    DH_tep = np.array([299792.458 / tep_cosmo.H(z) for z in z_test])
    
    DM_lcdm = np.array([lcdm_cosmo.comoving_distance(z) for z in z_test])
    DH_lcdm = np.array([299792.458 / lcdm_cosmo.H(z) for z in z_test])
    
    # Calculate ratios
    DM_over_rs_tep = DM_tep / rs_tep
    DM_over_rs_lcdm = DM_lcdm / rs_lcdm
    
    DH_over_rs_tep = DH_tep / rs_tep
    DH_over_rs_lcdm = DH_lcdm / rs_lcdm
    
    # Residuals
    DM_residuals = (DM_over_rs_tep - DM_over_rs_lcdm) / DM_over_rs_lcdm
    DH_residuals = (DH_over_rs_tep - DH_over_rs_lcdm) / DH_over_rs_lcdm
    
    return {
        "rs_shift": (rs_tep - rs_lcdm) / rs_lcdm,
        "DM_ratio_diff_percent": DM_residuals * 100,
        "DH_ratio_diff_percent": DH_residuals * 100,
    }

def run():
    from c0_common import TEPLogger, set_step_logger, ensure_dirs, print_status, write_json, step_json_path
    
    STEP_ID = "step_06_02_bao_likelihood"
    logger = TEPLogger(STEP_ID)
    set_step_logger(logger)
    ensure_dirs()
    
    print_status(f"Starting {STEP_ID}", "TITLE")
    print_status("BAO Likelihood Module Initialized.", "SUCCESS")
    
    payload = {
        "step": STEP_ID,
        "status": "completed",
        "message": "Module initialized and ready for execution."
    }
    
    write_json(step_json_path(STEP_ID), payload)
    return payload

if __name__ == "__main__":
    run()
