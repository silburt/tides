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

void migration(double* tau_a, double* t_mig, double* t_damp, double *expmigfac, int* phi_i, double* max_t_mig, double* P, int i, double RT, double Ms, double mp, double iptmig_fac, double a, double afac, int p_suppress);

void assignQp(double* Qp, double Qpfac, double rp);

void special_cases(char* sysname, int i, double* mig_fac);

void calc_tidetau(double* tau_a, double* tau_e, double Qp, double mp, double rp, double Ms, double e_default, double a_default, char* sysname, int i, int p_suppress);

void printwrite(int i, char* txt_file, double a,double P,double e,double mp,double rp,double Qp,double tau_a,double t_mig, double t_damp,double afac,int p_suppress);

#endif /* defined(____migration__) */
