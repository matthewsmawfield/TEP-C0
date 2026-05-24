/*
 * TEP Background Header File
 * 
 * Add TEP parameters to CLASS background structure
 */

#ifndef __TEP_BACKGROUND__
#define __TEP_BACKGROUND__

/* TEP parameters to add to struct background in background.h:
 * 
 *   short tep_mode;       // _TRUE_ or _FALSE_
 *   double epsilon_T;     // TEP coupling amplitude (fit from SNe: ~0.174)
 *   double z_T;          // Characteristic redshift (fit from SNe: ~5.0)
 *   double n_T;          // Power-law index (default: 1.0)
 * 
 * These should be added to the background structure in include/background.h
 */

/* TEP function declarations */
double tep_f_transition(struct background *pba, double z);
double tep_gamma_factor(struct background *pba, double z);
double tep_hubble_rate(struct background *pba, double z, double H_lcdm);
int background_functions_tep(struct background *pba, double a, 
                           double *pvecback_integrate, double *pvecback,
                           double *pvecthermo, double *pvecmetric);
int tep_background_init(struct background *pba);

/* Default TEP parameter values */
#define TEP_EPSILON_T_DEFAULT 0.0
#define TEP_Z_T_DEFAULT 1.0
#define TEP_N_T_DEFAULT 1.0

/* Input string mappings for TEP parameters */
#define TEP_MODE_STRING "tep_mode"
#define TEP_EPSILON_T_STRING "tep_epsilon_T"
#define TEP_Z_T_STRING "tep_z_T"
#define TEP_N_T_STRING "tep_n_T"

#endif /* __TEP_BACKGROUND__ */
