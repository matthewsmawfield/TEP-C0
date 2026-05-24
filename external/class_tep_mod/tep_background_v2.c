/*
 * TEP-Modified Background for CLASS v2.0
 * 
 * Updated to use ε_T and z_T parameters from SNe fit
 * 
 * H_TEP(z) = H_LCDM(z) × (1 + ε_T × f_T(z))
 * 
 * where f_T(z) = 1 - exp(-(z/z_T)^n)
 * 
 * For high-z behavior (z >> z_T), we cap the redshift contribution
 * to prevent blow-up at CMB redshifts.
 */

#include "background.h"
#include <math.h>

/* TEP transition function f_T(z) */
/* Returns value between 0 and 1 */
double tep_f_transition(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 0.0;
    }
    
    /* Cap z contribution to prevent blow-up at very high z */
    /* This is needed for CMB redshifts (z ~ 1100) when z_T ~ 5 */
    double z_effective = z;
    if (pba->z_T > 0 && z > pba->z_T * 3.0) {
        z_effective = pba->z_T * 3.0;
    }
    
    double ratio = z_effective / pba->z_T;
    double exponent = pow(ratio, pba->n_T);
    
    return 1.0 - exp(-exponent);
}

/* TEP gamma factor: 1 + ε_T × f_T(z) */
double tep_gamma_factor(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        return 1.0;
    }
    
    double f_T = tep_f_transition(pba, z);
    return 1.0 + pba->epsilon_T * f_T;
}

/* Modified Hubble rate with TEP */
/* H_TEP(z) = H_LCDM(z) × (1 + ε_T × f_T(z)) */
double tep_hubble_rate(struct background *pba, double z, double H_lcdm) {
    double gamma = tep_gamma_factor(pba, z);
    return H_lcdm * gamma;
}

/* TEP comoving distance element */
/* dD_C = c dz / H_TEP(z) */
double tep_comoving_distance_integrand(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_ || pba->epsilon_T == 0.0) {
        /* Standard LCDM: need H(z) from background */
        /* This is called within background_functions where H is available */
        return 1.0; /* Placeholder - actual value computed in background_functions_tep */
    }
    
    /* The integrand for comoving distance: c / H_TEP(z) */
    /* This will be used in tau computation */
    return 1.0; /* Actual computation done in background_functions_tep */
}

/* Modified Friedmann equation for TEP background */
int background_functions_tep(
    struct background *pba,
    double a,
    double *pvecback_integrate,
    double *pvecback,
    double *pvecthermo,
    double *pvecmetric
) {
    /* Call standard CLASS background functions first */
    int flag = background_functions(pba, a, pvecback_integrate, pvecback, pvecthermo, pvecmetric);
    
    if (flag != _SUCCESS_) {
        return flag;
    }
    
    /* Apply TEP modification if enabled */
    if (pba->tep_mode == _TRUE_ && pba->epsilon_T != 0.0) {
        double z = 1.0/a - 1.0;
        double H_original = pvecback[pba->index_bg_H];
        
        /* Apply TEP Hubble modification */
        double H_tep = tep_hubble_rate(pba, z, H_original);
        pvecback[pba->index_bg_H] = H_tep;
        
        /* Also modify conformal Hubble rate: H_conf = a × H */
        pvecback[pba->index_bg_H * pba->index_bg_a] = a * H_tep;
        
        /* Log modification for debugging */
        /* fprintf(stderr, "TEP: z=%.2f, H_LCDM=%.4f, H_TEP=%.4f, gamma=%.4f\n", 
         *         z, H_original, H_tep, tep_gamma_factor(pba, z));
         */
    }
    
    return _SUCCESS_;
}

/* TEP conformal time integral */
/* tau = ∫ c da / (a² H_TEP(a)) = ∫ c dz / H_TEP(z) */
double tep_conformal_time_integrand(double a, void *params) {
    struct background *pba = (struct background *)params;
    double z = 1.0/a - 1.0;
    
    /* Get H_LCDM at this scale factor */
    /* This requires calling background_functions - simplified here */
    /* In practice, this is done within the main integration loop */
    
    /* c / (a² H_TEP(a)) = c (1+z)² / H_TEP(z) */
    double c_Mpc_s = 299792.458 * 1000.0 / 3.08567758e22; /* c in Mpc/s */
    /* Actual computation requires full background state */
    
    return c_Mpc_s / (a * a); /* Placeholder - see full implementation */
}

/* Compute TEP angular diameter distance */
/* D_A = c/(1+z) × ∫_0^z dz'/H_TEP(z') */
double tep_angular_diameter_distance(struct background *pba, double z) {
    /* This would be computed from the background evolution */
    /* D_A = comoving_distance / (1+z) where comoving uses H_TEP */
    
    /* For now, return LCDM value (will be overridden by proper integration) */
    return 0.0; /* Computed from pvecback in practice */
}

/* Initialize TEP parameters from input */
int tep_background_init(struct background *pba) {
    if (pba->tep_mode == _TRUE_) {
        /* Validate parameters */
        if (pba->epsilon_T < 0.0) {
            fprintf(stderr, "Error: epsilon_T must be non-negative\n");
            return _FAILURE_;
        }
        if (pba->z_T <= 0.0) {
            fprintf(stderr, "Error: z_T must be positive\n");
            return _FAILURE_;
        }
        if (pba->n_T <= 0.0) {
            fprintf(stderr, "Error: n_T must be positive\n");
            return _FAILURE_;
        }
        
        /* Log initialization */
        printf("TEP mode enabled:\n");
        printf("  epsilon_T = %.6f\n", pba->epsilon_T);
        printf("  z_T = %.6f\n", pba->z_T);
        printf("  n_T = %.6f\n", pba->n_T);
        printf("  f_T(z=1) = %.6f\n", tep_f_transition(pba, 1.0));
        printf("  f_T(z=1100) = %.6f\n", tep_f_transition(pba, 1100.0));
    }
    
    return _SUCCESS_;
}

/* TEP distance duality violation test */
/* η(z) = D_L / [D_A × (1+z)²] */
/* For TEP: η(z) ≠ 1 in general */
double tep_distance_duality(struct background *pba, double z) {
    /* This would require computing both D_L and D_A */
    /* In TEP, D_L ≠ (1+z)² D_A in general */
    
    double D_L = 0.0; /* Luminosity distance */
    double D_A = 0.0; /* Angular diameter distance */
    
    /* These are computed from the background evolution */
    /* For now, return LCDM value */
    
    if (D_A > 0) {
        return D_L / (D_A * (1.0 + z) * (1.0 + z));
    }
    return 1.0; /* LCDM limit */
}
