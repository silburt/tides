//
//  readplanets.h
//  
//
//  Created by Ari Silburt on 2015-01-07.
//
//

#ifndef _readplanets_h
#define _readplanets_h

void readplanets(char *sysname, char *charac_txt, int *char_pos, int *_N, double *Ms, double *Rs, double *mp, double *rp, double* P, int p_suppress);

void extractplanets(int* char_pos, double *mp, double *rp, double* P, int p_suppress);

void naming(char* sysname, char* txt, double K, double iptmig_fac, double e_ini, double Qpfac, int tide_force);

void printwrite(int i, char* txt_file, double a,double P,double e,double mp,double rp,double Qp,double tau_a,double t_mig, double t_damp,double afac,int p_suppress);

#endif
