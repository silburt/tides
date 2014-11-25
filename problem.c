/**
 * @file 	problem.c
 * @brief 	Example problem: forced migration of GJ876.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Willy Kley <kley@uni-tuebingen.de>
 * @detail 	This example applies dissipative forces to two
 * bodies orbiting a central object. The forces are specified
 * in terms of damping timescales for the semi-major axis and
 * eccentricity. This mimics planetary micration in a protostellar disc. 
 * The example reproduces the study of Lee & Peale (2002) on the 
 * formation of the planetary system GJ876. For a comparison, 
 * see figure 4 in their paper. The IAS15 integrator is used 
 * because the forces are velocity dependent.
 * Special thanks goes to Willy Kley for helping me to implement
 * the damping terms as actual forces. 
 *
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "tools.h"
#include "problem.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
void problem_migration_forces();
double* Qp;     /**< Q'=k_2/Q used for tides. p for prime> */
double* radius;
double* m;
double* a;
double* e;
double* w;
double* f;

#ifdef OPENGL
extern int display_wire;
#endif 	// OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= 1e-2*2.*M_PI;		// in year/(2*pi)
	boxsize 	= 3;			// in AU
	//tmax		= 4.5e4*2.*M_PI;	// in year/(2*pi)
    tmax        = 50000;
    
    K           = 100;      //tau_a/tau_e ratio. I.e. Lee & Peale (2002)
    T           = 2.*M_PI*20000.0;  //tau_a;
    t_mig[0]    = 30000.; //migration times
    t_mig[1]    = 35000.; //migration damp out over 5000 years.
#ifdef OPENGL
	display_wire 	= 1;			
#endif 	// OPENGL
	init_box();

	// Initial conditions
	// Parameters are those of Lee & Peale 2002, Figure 4. 
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1.1;			// This is a sub-solar mass star
	particles_add(star);
    
    int _N = 1; //number of planets
    m      = calloc(sizeof(double),_N);
    a      = calloc(sizeof(double),_N);
    e      = calloc(sizeof(double),_N);
    w      = calloc(sizeof(double),_N);
    f      = calloc(sizeof(double),_N);
    //Fix later, planet 1 and 2
    m[0]=1e-3;
    a[0]=0.05;
    e[0]=0.15;
    w[0]=125.*(M_PI/180.);
    f[0]=M_PI/4.;
    
    //m[1]=1.89e-5;
    //a[1]=0.44;
    //e[1]=0.05;
    //w[1]=5.*M_PI/4.;
    //f[1]=0;
    
    for(int i=0;i<_N;i++){
        struct particle p = tools_init_orbit2d(star.m,m[i],a[i],e[i],w[i],f[i]);
        particles_add(p);
    }

/*
	struct particle p1;		// Planet 1
	p1.x 	= 0.25;	p1.y = 0;	p1.z = 0;
	p1.ax 	= 0;	p1.ay = 0; 	p1.az = 0;
	p1.m  	= 0.56e-3;
	p1.vz 	= 0;
	p1.vx 	= 0;	p1.vy = sqrt(G*(star.m+p1.m)/p1.x);;
	particles_add(p1); 
	
	struct particle p2;		// Planet 2
	p2.x 	= 0.52;	p2.y = 0; 	p2.z = 0;
	p2.ax 	= 0;	p2.ay = 0; 	p2.az = 0;
	p2.m  	= 1.89e-3;
	p2.vz 	= 0;
	p2.vx 	= 0;	p2.vy = sqrt(G*(star.m+p2.m)/p2.x);
	particles_add(p2); 
*/
 
    tau_a  = calloc(sizeof(double),N);  //migration of semi-major axis
	tau_e  = calloc(sizeof(double),N);  //migration (damp) of eccentricity
    Qp     = calloc(sizeof(double),N);  //Q' = k_2/Q
    radius = calloc(sizeof(double),N);

    //Need to find a better way to arrange this
    //tau_a[2] = T;	// Migration timescale of planet 2 is 20000 years.
	//tau_e[2] = tau_a[2]/K;      // Eccentricity damping timescale.
    Qp[1]    = 0.35/1.0e4;     // = k_2/Q, k_2 mercury~0.5 (Padovan et al.), Jupiter=Gavrilov et al.
    //Qp[2]    = 0.0127/.5e3;
    radius[0]=0.1;             //in Rs
    radius[1]=0.05;
    //radius[2]=0.04;
    
	problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.
#ifndef INTEGRATOR_WH			// The WH integrator assumes a heliocentric coordinate system.
	tools_move_to_center_of_momentum();  		
#endif // INTEGRATOR_WH

	system("rm -v orbits.txt"); // delete previous output file
}

void problem_migration_forces(){
    //ramp down the migration force
    if (t > t_mig[0] && t < t_mig[1]) {
        tau_a[2] = T + (t - t_mig[0])*(2.*M_PI*200000.0 - T)/(t_mig[1] - t_mig[0]);
        tau_e[2] = tau_a[2]/K;
    } else if(t > t_mig[1]){
        tau_a[2]=0.;
        tau_e[2]=0.;
    }
    
    struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
	for(int i=1;i<N;i++){
		if (tau_e[i]!=0||tau_a[i]!=0){
			struct particle* p = &(particles[i]);
			const double dvx = p->vx-com.vx;
			const double dvy = p->vy-com.vy;
			const double dvz = p->vz-com.vz;

			if (tau_a[i]!=0){ 	// Migration
				p->ax -=  dvx/(2.*tau_a[i]);
				p->ay -=  dvy/(2.*tau_a[i]);
				p->az -=  dvz/(2.*tau_a[i]);
			}
			if (tau_e[i]!=0){ 	// Eccentricity damping
				const double mu = G*(com.m + p->m);
				const double dx = p->x-com.x;
				const double dy = p->y-com.y;
				const double dz = p->z-com.z;

				const double hx = dy*dvz - dz*dvy; 
				const double hy = dz*dvx - dx*dvz;
				const double hz = dx*dvy - dy*dvx;
				const double h = sqrt ( hx*hx + hy*hy + hz*hz );
				const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
				const double r = sqrt ( dx*dx + dy*dy + dz*dz );
				const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
				const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
				const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
				const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
				const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
				const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
				p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
                //printf("a=%f,e=%f,",a,e);
			}
		}
		com = tools_get_center_of_mass(com,particles[i]);
	}
}

void problem_inloop(){
}

void problem_output(){
    struct particle com = particles[0];
    //Tides
    for(int i=1;i<N;i++){
        struct particle* par = &(particles[i]);
        const double m = par->m;
        const double mu = G*(com.m + m);
        //radius of planet must be in AU for units to work out since G=1, [t]=yr/2pi, [m]=m_star
        const double rp = radius[i]*0.00464913; //converts from Rs to AU units
        
        const double dvx = par->vx-com.vx;
        const double dvy = par->vy-com.vy;
        const double dvz = par->vz-com.vz;
        const double dx = par->x-com.x;
        const double dy = par->y-com.y;
        const double dz = par->z-com.z;
        
        const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
        const double r = sqrt ( dx*dx + dy*dy + dz*dz );
        const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
        const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
        const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
        const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
        
        double e = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
        double a = -mu/( v*v - 2.*mu/r );			// semi major axis
        double n = sqrt(mu/(a*a*a));
        double w = atan2(ey,ex);              // seems right upon check -  Fundamentals of Astrodynamics and
                                              //Applications, by Vallado, 2007
        if(ey < 0.) w = 2*M_PI + w;
        const double rdote = dx*ex + dy*ey + dz*ez;
        const double cosf = rdote/(e*r);
        double f = acos(cosf);                  // seems right upon check
        if(vr < 0.) f = 2*M_PI - f;
        double const sinf = sin(f);
        double const sinwf = sin(w+f);
        double const coswf = cos(w+f);
        
        //Eccentric Anomaly
        //double terme = sqrt((1+e)/(1-e));
        //double const E = 2*atan(tan(f/2)/terme);
        
        //Tides
        /*
        const double R5 = rp*rp*rp*rp*rp;
        const double factor = sqrt(mu*(com.m+m)*(com.m+m))*R5*Qp[i]/m;
        const double factor2 = pow(a,11.0/2.0);
        const double da = -21.0*dt*factor*e*e/factor2;              //Tidal change for a
        const double de = -(21.0/2.0)*dt*factor*e/(factor2*a);      //Tidal change for e
        if(i==1 && t < 5.)printf("INI tides: da=%.15f,de=%.15f,a=%f,e=%f \n",da,de,a,e);
        if(i==1 && t > 49990.)printf("FINI tides: da=%.15f,de=%.15f,a=%f,e=%f \n",da,de,a,e);
        a += da;
        e += de;
        if(e < 0.) e=1e-5;
        */
        
        
        //Re-update coords.
        const double r_new = a*(1 - e*e)/(1 + e*cosf);
        n = sqrt(mu/(a*a*a));
        par->x = r_new*coswf + com.x;
        par->y = r_new*sinwf + com.y;
        const double term = n*a/sqrt(1-e*e);
        const double rdot = term*e*sinf;
        const double rfdot = term*(1 + e*cosf);
        par->vx = rdot*coswf - rfdot*sinwf + com.vx;
        par->vy = rdot*sinwf + rfdot*coswf + com.vy;
        
        com = tools_get_center_of_mass(com,particles[i]);  //Does this need to happen before updating par?
        
        //if(i==1 && t<5.)printf("DELTA INI: x,y,vx,vy=%.15f,%.15f,%.15f,%.15f \n ", r_new*coswf - dx, r_new*sinwf -dy, rdot*coswf - rfdot*sinwf -dvx, rdot*sinwf + rfdot*coswf -dvy);
        //if(i==1 && t > 9990.)printf("DELTA FINI: x,y,vx,vy=%.15f,%.15f,%.15f,%.15f \n ", r_new*coswf - dx, r_new*sinwf -dy, rdot*coswf - rfdot*sinwf -dvx, rdot*sinwf + rfdot*coswf -dvy);
    }

    
	if(output_check(10000.*dt)){
		output_timing();
	}
	if(output_check(40.)){
        //A.S. - append orbits in orbits.txt. Ordering of outputs goes:
        //time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
		output_append_orbits("orbits.txt");
#ifndef INTEGRATOR_WH
		tools_move_to_center_of_momentum();  			// The WH integrator assumes a heliocentric coordinate system.
#endif // INTEGRATOR_WH
	}
}

void problem_finish(){
}
