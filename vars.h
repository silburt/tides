//
//  vars.h
//  
//
//  Created by Ari Silburt on 2015-04-13.
//
//

#ifndef _vars_h
#define _vars_h

// A.S. variables added
double  K;              /**<tau_a/tau_e>*/
int     _N;             /**<# of planets>*/
int     tides_on;       /**<Parameter to control if tides on/off>**/
int     tide_force;     /**<If == 1, implement tides as forces, not a' & e'>**/
double  tide_delay;     /**<Lag time (in years) to turn on tides after**/
int     mig_forces;     /**<Parameter to control if migration on/off>**/
double  afac;           /**<Factor to increase a by of all planets>**/
int     p_suppress;     /**<If = 1, then suppress all printing>**/
double* tau_a;          /**< Migration timescale in years for all particles */
double* tau_e;          /**< Eccentricity damping timescale in years for all particles */
double* lambda;         /**<Resonant angle>**/
double* omega;          /**<argument of periapsis>**/
double* t_mig;          /**<Migration timescale calc according to Goldreich & Schlichting (2014)>**/
double* t_damp;
double* mu_a;
double* en;             /**<mean motion array - for pendulum energy>**/
double* term1;          /**<gold & schlich**/
double* term2a;          /**<gold & schlich**/
double* coeff2;
char* c;
int* phi_i;
int tide_print;         /**<print message when tides are turned on>**/
char txt_file[80];

double* tidetau_e;       /**<tau_e/a if calculating tides as forces via Papaloizou & Larwood>**/
double* tidetau_a;

#endif