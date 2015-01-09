/**
Code developed by Ari Silburt to evolve Kepler planets under the influence of tides, with an 
initial migration to put planets into resonance.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "main.h"
#include "tools.h"
#include "problem.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"
#include "../examples/tides/readplanets.h"
#include "../examples/tides/assignparams.h"

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double *lambda; /**<Resonant angle>**/
double *omega;  /**<argument of periapsis>**/
int tide_print; /**<print message when tides are turned on>**/
void problem_migration_forces();

#ifdef OPENGL
extern int display_wire;
#endif 	// OPENGL

void problem_init(int argc, char* argv[]){
	/* Setup constants */
    //dt = (dt is calc in readplanets.c), unit is yr/2PI
	boxsize 	= 3;                // in AU
	//tmax		= 1e7*2.*M_PI;      // in year/(2*pi)
    tmax        = 100000.;
    
    K           = 100;              //tau_a/tau_e ratio. I.e. Lee & Peale (2002)
    T           = 2.*M_PI*300000.;  //tau_a, typical timescale=20,000 years;
    t_mig       = 20000.;           //Begin damping migration at t_mig.
    t_damp      = 30000.;           //length of migration damping. Afterwards, no migration.
    tide_forces = 1;                //If ==0, then no tidal forces on planets.
    tide_delay  = 0.;               //Lag time after which tidal forces are turned on. Requires tide_forces=1!!
    mig_forces  = 1;                //If ==0, no migration.
    afac        = 1.03;              //Factor to increase 'a' of OUTER planets by.
    char c[20]  = "TESTP10J";           //System being investigated
    txt_file    = "runs/orbits_temp2.txt";           //Where to store orbit outputs
    
#ifdef OPENGL
	display_wire 	= 1;			
#endif 	// OPENGL
	init_box();
    
    //Delete previous file if it exists.
    char sys_arg[50] = "rm -v ";
    strcat(sys_arg,txt_file);
    system(sys_arg); // delete previous output file
    tide_print = 0;
    
    // Initial conditions
    printf("You have chosen: %s \n",c);
    double Ms,Rs,a,rho,inc,mp,rp,tau_atemp,Qp_temp;
    int char_val, _N;
    
    //**Initial eccentricity**
    const double f=0., w=M_PI/2., e=0.1;
    readplanets(c,txt_file,&char_val,&_N,&Ms,&Rs,&a,&rho,&inc,&mp,&rp,&dt);
    struct particle star; //Star MUST be the first particle added.
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = Ms;
	particles_add(star);
    
    //Extra slot for star
    tau_a  = calloc(sizeof(double),_N+1);  //migration of semi-major axis
	tau_e  = calloc(sizeof(double),_N+1);  //migration (damp) of eccentricity
    lambda = calloc(sizeof(double),_N+1);  //resonant angle for each planet
    omega = calloc(sizeof(double),_N+1);  //argument of periapsis for each planet
    
    struct particle p = tools_init_orbit2d(Ms, mp, a, e, w, f);
    p.r = rp;
    assignparams(&tau_atemp,&Qp_temp,mp,rp,T,txt_file);
    p.Qp=Qp_temp;
    particles_add(p);
    printf("System Properties: # planets=%d, Rs=%f, Ms=%f \n",_N, Rs, Ms);
    printf("Planet 1: a=%f,mp=%f,rp=%f,Qp=%f \n",a,mp,rp,Qp_temp);
    
    for(int i=1;i<_N;i++){
        extractplanets(&char_val,&a,&rho,&inc,&mp,&rp);
        a *= afac;      //Increase 'a' of outer planets by afac
        struct particle p = tools_init_orbit2d(Ms, mp, a, e, w, f);
        p.r = rp;
        assignparams(&tau_atemp,&Qp_temp,mp,rp,T,txt_file);
        p.Qp=Qp_temp;
        tau_a[i+1]=tau_atemp;
        tau_e[i+1]=tau_atemp/K;
        
        particles_add(p);
        printf("Planet %i: a=%f,mp=%f,rp=%f,Qp=%f,tau_a=%f,afac=%f, \n",i+1,a,mp,rp,Qp_temp,tau_a[i+1],afac);
    }
    
	problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.
#ifndef INTEGRATOR_WH			// The WH integrator assumes a heliocentric coordinate system.
	tools_move_to_center_of_momentum();  		
#endif // INTEGRATOR_WH

}

void problem_migration_forces(){
    if(mig_forces==1){
        //ramp down the migration force (by increasing the migration timescale)
        double t_mig2 = t_mig + t_damp; //when migration ends
        if (t > t_mig && t < t_mig2) {
            tau_a[2] = T + (t - t_mig)*(1000000.0 - T)/(t_mig2 - t_mig);
            tau_e[2] = tau_a[2]/K;
        } else if(t > t_mig2){
            tau_a[2]=0.;
            tau_e[2]=0.;
        }
        struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
        for(int i=1;i<N;i++){ // N = _N + 1 = total number of planets + star
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
                }
            }
            com = tools_get_center_of_mass(com,particles[i]);
        }
    }
    
}

void problem_inloop(){
}

void problem_output(){
    //Calculate Orbital Elements
    struct particle com = particles[0];
    for(int i=1;i<N;i++){
        struct particle* par = &(particles[i]);
        const double m = par->m;
        const double mu = G*(com.m + m);
        //radius of planet must be in AU for units to work out since G=1, [t]=yr/2pi, [m]=m_star
        const double rp = par->r*0.00464913;       //Rp from Solar Radii to AU
        const double Qp = par->Qp;

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
        //double a = -mu/( v*v - 2.*mu/r );			// semi major axis
        //const double cosE = (a - r)/(a*e);
        //double cosf = (1 - e*e)/(e - e*e*cosE) - 1/e;
            
        // true anomaly + periapse (wiki, Fund. of Astrodyn. and App., by Vallado, 2007)
        const double rdote = dx*ex + dy*ey + dz*ez;
        double cosf = rdote/(e*r);
        //double cosf = (a*(1 - e*e) - r)/(r*e);
        if(cosf >= 1.) cosf = 1.;
        if(cosf <= -1.) cosf = -1.;
        //double const sinf = sin(f);
        double sinf = sqrt(1. - cosf*cosf);
        if(vr < 0.) sinf *= -1.;
        double const sinwf = dy/r;
        double const coswf = dx/r;
        double a = r*(1. + e*cosf)/(1. - e*e);
        double n;
            
        //Eccentric Anomaly
        //double terme = sqrt((1+e)/(1-e));
        //double const E = 2*atan(tan(f/2)/terme);
        //const double cosf = (cos(E) - e)/(1 - e*cos(E));
        //double a = r/(1 - e*cos(E));
        //double n = sqrt(mu/(a*a*a));
        
        //Tides
        if(tide_forces==1 && t > tide_delay){
            const double a2 = a*a;
            const double rp2 = rp*rp;
            const double R5a5 = rp2*rp2*rp/(a2*a2*a);
            const double GM3a3 = sqrt(G*com.m*com.m*com.m/(a2*a));
            const double de = -dt*(9.*M_PI*0.5)*Qp*GM3a3*R5a5*e/m;   //Tidal change for e
            const double da = 2.*a*e*de;                             //Tidal change for a
        
            a += da;
            e += de;
        
            //Re-update coords.
            const double r_new = a*(1. - e*e)/(1. + e*cosf);
            //const double r_new = a*(1 - e*cosE);
            
            const double x_new = r_new*coswf + com.x;
            const double y_new = r_new*sinwf + com.y;
            n = sqrt(mu/(a*a*a));
            
            const double term = n*a/sqrt(1.- e*e);
            const double rdot = term*e*sinf;
            const double rfdot = term*(1. + e*cosf);
            const double vx_new = rdot*coswf - rfdot*sinwf + com.vx;
            const double vy_new = rdot*sinwf + rfdot*coswf + com.vy;
            
            par->x = x_new;
            par->y = y_new;
            par->vx = vx_new;
            par->vy = vy_new;
            com = tools_get_center_of_mass(com,particles[i]);
            
            //Stop program if nan values being produced.
            if(x_new!=x_new || y_new!=y_new || vx_new!=vx_new ||vy_new!=vy_new){
                printf("\n cartesian before: dx=%f,dy=%f,dz=%f,ex=%f,ey=%f,ex=%f,r=%f,vx=%f,vy=%f,com.vx=%f,com.vy=%f,v=%f \n",dx,dy,dz,ex,ey,ez,r,par->vx,par->vy,com.vx,com.vy,v);
                printf("Orbital elements: mu=%f,e=%f,a=%f,cosf=%.22f,dt=%f,de=%f,da=%f,GM3a3=%f,R5a5=%f \n",mu,e,a,cosf,dt,de,da,GM3a3,R5a5);
                printf("\n cartesian after: x_new=%f,y_new=%f,vx_new=%f,vy_new=%f,term=%f,rdot=%f,rfdot=%f \n",x_new,y_new,vx_new,vy_new,term,rdot,rfdot);
                exit(0);
            }
            
            //print message
            if(tide_print == 0){
                printf("\n ***Tides have just been turned on at t=%f years***\n",t);
                tide_print = 1;
            }
            
        } else { n = sqrt(mu/(a*a*a));}     //Still need to calc this for period. 
        
        if(output_check(tmax/2500.)){
            omega[i] = atan2(ey,ex);
            //if(ey < 0.) omega[i] += 2*M_PI;
            double cosE = (a - r)/(a*e);
            double E;
            if(cosf > 1. || cosf < -1.){
                E = M_PI - M_PI*cosE;
            } else {
                E = acos(cosE);
            }
            if(vr < 0.) E = 2.*M_PI - E;
            double MA = E - e*sin(E);
            lambda[i] = MA + omega[i];
            double phi = 0.;     //resonant angle
            if(i>=2) phi = 2.*lambda[i] - lambda[i-1] - omega[i-1]; //2:1 resonance
            while(phi >= 2*M_PI) phi -= 2*M_PI;
            while(phi < 0.) phi += 2*M_PI;
                
            //output orbits in txt_file.
            FILE *append;
            append=fopen(txt_file, "a");
            //output order = time(yr/2pi),a(AU),e,P(days),arg. of peri., mean anomaly,
            //               eccentric anomaly, mean longitude, resonant angle
            fprintf(append,"%e\t%.10e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,a,e,365./n,omega[i],MA,E,lambda[i],phi);
            fclose(append);
            
            #ifndef INTEGRATOR_WH
            tools_move_to_center_of_momentum();  			// The WH integrator assumes a heliocentric coordinate system.
            #endif // INTEGRATOR_WH
        }
    }

	if(output_check(10000.*dt)){
		output_timing();
	}
}

void problem_finish(){
}
