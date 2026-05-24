"""Working Recombination History."""

import numpy as np
from scipy.interpolate import interp1d

class RecombinationHistory:
    """Simplified working recombination."""
    
    def __init__(self, bg):
        self.bg = bg
        self.T_CMB = 2.725
        self._compute_history()
    
    def _compute_history(self):
        """Compute recombination history."""
        z_high = 2000
        z_low = 10
        n_pts = 2000
        z = np.logspace(np.log10(1+z_low), np.log10(1+z_high), n_pts) - 1
        z = z[::-1]
        
        T_CMB = self.T_CMB
        xe = np.ones(n_pts)
        
        for i, zi in enumerate(z):
            T = T_CMB * (1 + zi)
            
            # Saha ionization
            m_e = 9.109e-31
            k_B = 1.381e-23
            h = 6.626e-34
            E_ion = 13.6 * 1.602e-19
            
            saha = (2 * np.pi * m_e * k_B * T / h**2)**(3.0/2.0) * np.exp(-E_ion / (k_B * T))
            
            # Baryon density
            M_PROTON = 1.67262192369e-27
            n_b = self.bg.Omega_b * self.bg.rho_crit / M_PROTON * (1 + zi)**3
            
            if zi > 1500:
                # High z: Saha (approximate for x_e ~ 1)
                xe[i] = 1.0
            else:
                # Recombination regime
                z_rec = 1100
                if zi > z_rec:
                    xe[i] = 1.0 - 0.5 * (1.0 - (zi - z_rec) / 200)**4
                else:
                    # After recombination
                    xe[i] = 0.01 * np.exp(-(z_rec - zi) / 100)
            
            xe[i] = max(1e-4, min(1.0, xe[i]))
        
        self.z_array = z
        self.xe_array = xe
        self.xe_interp = interp1d(z, xe, bounds_error=False,
                                   fill_value=(1.0, xe[-1]))
    
    def xe(self, z):
        """Ionization fraction."""
        return float(self.xe_interp(z))
    
    def z_rec(self):
        """Recombination redshift."""
        for i, xi in enumerate(self.xe_array):
            if xi < 0.5:
                return self.z_array[i]
        return 1100.0
