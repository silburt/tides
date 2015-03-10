//
//  readplanets.h
//  
//
//  Created by Ari Silburt on 2015-01-07.
//
//

#ifndef _readplanets_h
#define _readplanets_h

void readplanets(char *sysname, char *charac_txt, int *char_pos, int *_N, double *Ms, double *Rs, double *rho, double *inc, double *mp, double *rp, double* P, double *dt, double timefac, int p_suppress);

void extractplanets(int* char_pos, double *rho, double *inc, double *mp, double *rp, double* P, int p_suppress);

void calcsemi(double* a, double Ms, double P);

#endif
