/*
 * TEP-Modified Background for CLASS
 * 
 * This file modifies the background Hubble rate to include TEP effects:
 * H_TEP(z) = H_LCDM(z) * Gamma_TEP(z)
 * 
 * where Gamma_TEP(z) = exp(Sigma_0 * c/H0 * ln(1+z))
 */

#include "background.h"
#include <math.h>

/* TEP parameters - added to background structure */
/* In background.h, add to struct background:
 *   short tep_mode;
 *   double Sigma_0;
 *   double A_env;
 */

/* TEP path enhancement factor */
double tep_gamma_factor(struct background *pba, double z) {
    if (pba->tep_mode == _FALSE_) {
        return 1.0;
    }
    
    double c_kms = 299792.458;  /* Speed of light in km/s */
    double ln_gamma = pba->Sigma_0 * c_kms / pba->H0 * log(1.0 + z);
    return exp(ln_gamma);
}

/* Modified Hubble rate with TEP */
double tep_hubble_rate(struct background *pba, double z, double H_lcdm) {
    double gamma = tep_gamma_factor(pba, z);
    return H_lcdm * gamma;
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
    /* Call standard CLASS background functions */
    int flag = background_functions(pba, a, pvecback_integrate, pvecback, pvecthermo, pvecmetric);
    
    if (flag != _SUCCESS_) {
        return flag;
    }
    
    /* Apply TEP modification if enabled */
    if (pba->tep_mode == _TRUE_) {
        double z = 1.0/a - 1.0;
        double H_original = pvecback[pba->index_bg_H];
        pvecback[pba->index_bg_H] = tep_hubble_rate(pba, z, H_original);
    }
    
    return _SUCCESS_;
}
