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
double* expmigfac;      /**<exponential migration factor>**/
int tautide_force_calc=0; /**<switch to calc tau_a and tau_e for tidal forces (post migration)>**/
char* Keplername;
int* phi_i;
int tide_print;         /**<print message when tides are turned on>**/
char txt_file[80];
int N_ini;              /**<initial number of planets - i.e. was there a collision?>**/
int collision_print = 0;/**<message to print if there's a collision>**/

double* tidetauinv_e;   /**<1/tau_e if calculating tides as forces via Papaloizou & Larwood>**/
double* tidetau_a;

#endif
