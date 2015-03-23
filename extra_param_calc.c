//
//  extra_param_calc.c
//  
//
//  Created by Ari Silburt on 2015-03-23.
//
//

//This file just stores all the additional ways to calculate orbital parameters
//to clean up my problem.c without losing all the precious info.

#include <stdio.h>

double a = -mu/( v*v - 2.*mu/r );			// semi major axis
const double cosE = (a - r)/(a*e);

double cosf = (1 - e*e)/(e - e*e*cosE) - 1/e;
double cosf = (a*(1 - e*e) - r)/(r*e);
double sinf = sin(f);

//Eccentric Anomaly
double terme = sqrt((1+e)/(1-e));
double const E = 2*atan(tan(f/2)/terme);
const double cosf = (cos(E) - e)/(1 - e*cos(E));
double a = r/(1 - e*cos(E));
double n = sqrt(mu/(a*a*a));

double r_new = a*(1 - e*cosE);