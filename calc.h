//
//  migration.h
//  
//
//  Created by Ari Silburt on 2015-04-14.
//
//

#ifndef ____migration__
#define ____migration__

#include <stdio.h>

void migration(char* sysname, double* tau_a, double* t_mig, double* t_damp, double *expmigfac, int* phi_i, double* max_t_mig, double* P, int i, double RT, double Ms, double mp, double iptmig_fac, double a, double afac, int p_suppress);

void assignk2Q(double* k2, double* Q, double k2fac, double rp);

void calc_tidetau(double* tau_a, double* tau_e, double K, double Qp, double mp, double rp, double Ms, double e_default, double a_default, char* sysname, int i, int p_suppress);

void k2fac_check(char* sysname, double* k2fac);

void calcsemi(double* a, double Ms, double P);

#endif /* defined(____migration__) */
