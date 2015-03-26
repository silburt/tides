//
//  assignparams.h
//  
//
//  Created by Ari Silburt on 2015-01-07.
//
//

#ifndef _assignparams_h
#define _assignparams_h

void assignparams(double *Qp, double Qpfac, double mp, double rp, double* T,double* t_mig_var, double Ms,char *txt_file, double a, double a_f, double P, double migspeed_fac);

void special_cases(char* sysname, int i, double* mig_fac);

#endif
