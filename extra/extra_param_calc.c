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



//Tidal forces (Mignard), implemented in Rodriguez (2013)
double fx[N];
double fy[N];
double fz[N];
//radius of planet must be in AU for units to work out since G=1, [t]=yr/2pi, [m]=m_star
struct particle com = particles[0];
const double Rs5 = pow(com.r*0.00464913,5); //Rs from Solar Radii to AU
const double Qpstar = 0.028/1e6;        //stellar k/Q, Wu & Murray (2003)
double GM2 = G*com.m*com.m;
for(int i=1;i<N;i++){
    struct particle* p = &(particles[i]);
    const double m = p->m;
    const double mu = G*(com.m + m);
    const double rp = p->r*0.00464913;  //Rp from Solar Radii to AU
    const double Qp = p->Qp;
    
    const double dvx = p->vx-com.vx;
    const double dvy = p->vy-com.vy;
    const double dvz = p->vz-com.vz;
    const double dx = p->x-com.x;
    const double dy = p->y-com.y;
    const double dz = p->z-com.z;
    const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
    const double r = sqrt ( dx*dx + dy*dy + dz*dz );
    const double vr = (dx*dvx + dy*dvy + dz*dvz);
    const double a = -mu/( v*v - 2.*mu/r );
    
    //Mignard tidal forces (1979), implemented in Rodriguez (2013)
    double r2 = r*r;
    double r10 = pow(r2,5);         //planet-star distance, AU
    const double rp2 = pow(rp,2);
    double n = sqrt( mu/(a*a*a) );
    double Gmp2 = G*m*m;
    double kdt_p = Qp/n;                    //tidal lag time, planet
    double kdt_s = Qpstar/n;
    double coeffp = -3*kdt_p*GM2*rp2*rp2*rp/r10;   //planet
    double coeffs = 3*kdt_s*Gmp2*Rs5/r10;  //star
    fx[i] = (coeffp - coeffs)*(2*dx*vr + r2*dvx);  //assumes rot. speed of planet/star = 0
    fy[i] = (coeffp - coeffs)*(2*dy*vr + r2*dvy);
    fz[i] = (coeffp - coeffs)*(2*dz*vr + r2*dvz);
}

for(int i=1;i<N;i++){
    struct particle* p = &(particles[i]);
    const double m = p->m;
    const double mratio = (com.m + m)/(com.m * m);
    if(i==1 && tide_print == 0){
        printf("\n fx=%.16f, fy=%.16f,mratio=%f, p->ax=%.16f, p->ay=%.16f \n",fx[i]*mratio,fy[i]*mratio,mratio, p->ax, p->ay);
        tide_print = 1;
    }
    p->ax += mratio*fx[i];
    p->ay += mratio*fy[i];
    p->az += mratio*fz[i];
    for(int j=1;j<N;j++){
        if(j != i){
            p->ax += fx[j]/com.m;
            p->ay += fy[j]/com.m;
            p->az += fz[j]/com.m;
        }
    }
}
//print message
if(tide_print == 0 && p_suppress == 0){
    printf("\n ***Tides (forces) have just been turned on at t=%f years***\n",t);
    tide_print = 1;
}
