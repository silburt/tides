/**
Code developed by Ari Silburt to evolve Kepler planets under the influence of tides, with an 
initial migration to put planets into resonance
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "input.h"
#include "main.h"
#include "tools.h"
#include "problem.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"
#include "../examples/tides/readplanets.h"
#include "../examples/tides/assignparams.h"

// A.S. variables added
double   K;           /**<tau_a/tau_e>*/
int      _N;          /**<# of planets>*/
int      tide_forces; /**<Parameter to control if tides on/off>**/
double   tide_delay;  /**<Lag time (in years) to turn on tides after**/
int      mig_forces;  /**<Parameter to control if migration on/off>**/
double   afac;        /**<Factor to increase a by of all planets>**/
int      p_suppress;  /**<If = 1, then suppress all printing>**/
double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double* lambda; /**<Resonant angle>**/
double* omega;  /**<argument of periapsis>**/
double* t_mig;  /**<Migration timescale calc according to Goldreich & Schlichting (2014)>**/
double* t_damp;
double* mu_a;
char* c;
int* phi_i;
int tide_print; /**<print message when tides are turned on>**/
char txt_file[80];
//

void problem_migration_forces();

#ifdef OPENGL
extern int display_wire;
#endif 	// OPENGL

void problem_init(int argc, char* argv[]){
    /* Setup constants */
    //dt = (dt is calc in readplanets.c), unit is yr/2PI
	boxsize 	= 3;                // in AU
    tmax        = input_get_double(argc,argv,"tmax",10000000.);  // in year/(2*pi)
    K           = 100;              //tau_a/tau_e ratio. I.e. Lee & Peale (2002)
    tide_forces = 1;                //If ==0, then no tidal forces on planets.
    mig_forces  = 1;                //If ==0, no migration.
    afac        = 1.06;             //Factor to increase 'a' of OUTER planets by.
    c           = argv[1];          //System being investigated, Must be first string after ./nbody!
    p_suppress  = 0;                //If = 1, suppress all print statements
    double RT   = 0.05;             //Resonance Threshold - if abs(P2/2*P1 - 1) < RT, then close enough to resonance
    double res  = 2.0;              //Resonance of interest: e.g. 2.0 = 2:1, 1.5 = 3:2, etc.
    
    //double efac = atof(argv[2]);          //added eccentricity
    double efac = 0.0;
    //double timefac = atof(argv[2]);          //timestep factor -> readplanets.c
    double timefac = 15.0;
    
#ifdef OPENGL
	display_wire 	= 1;			
#endif 	// OPENGL
	init_box();
    
    //Orbit outputs in txt_file
    char* dir = "runs/orbits_";
    char* ext = ".txt";
    strcat(txt_file, dir);
    strcat(txt_file, c);
    //char* underscore = "_";
    //strcat(txt_file, underscore);
    //char* c2    = argv[2];
    //strcat(txt_file, c2);
    strcat(txt_file, ext);
    
    //Delete previous file if it exists.
    char sys_arg[50] = "rm -v ";
    strcat(sys_arg,txt_file);
    system(sys_arg); // delete previous output file
    tide_print = 0;
    
    // Initial conditions
    if(p_suppress == 0) printf("You have chosen: %s \n",c);
    double Ms,Rs,a,rho,inc,mp,rp,P,Qp_temp;
    int char_val, _N;
    
    //Star
    readplanets(c,txt_file,&char_val,&_N,&Ms,&Rs,&a,&rho,&inc,&mp,&rp,&P,&dt,timefac,p_suppress);
    if(mig_forces == 0 && p_suppress == 0) printf("--> Migration is *off* \n");
    if(tide_forces == 0 && p_suppress == 0) printf("--> Tides are *off* \n");
    struct particle star; //Star MUST be the first particle added.
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = Ms;
	particles_add(star);
    
    //Write tidal info to file
    double tide_delay_output = 0;
    if(tide_forces == 1) tide_delay_output = tide_delay;
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%f \n",tide_delay_output);
    fclose(write);
    
    //Arrays, Extra slot for star, calloc sets values to 0 already.
    tau_a  = calloc(sizeof(double),_N+1);   //migration of semi-major axis
	tau_e  = calloc(sizeof(double),_N+1);   //migration (damp) of eccentricity
    lambda = calloc(sizeof(double),_N+1);   //resonant angle for each planet
    omega = calloc(sizeof(double),_N+1);    //argument of periapsis for each planet
    t_mig = calloc(sizeof(double),_N+1);
    t_damp = calloc(sizeof(double),_N+1);
    phi_i = calloc(sizeof(int),_N+1);       //phi index (for outputting resonance angles)
    mu_a = calloc(sizeof(double),_N+1);
    
    //Resonance vars
    double Period[_N],a_f; //a_f = final desired position of planet after mig
    Period[0] = 2.*M_PI*P/365.; //in yr/2pi (i.e. need to multiply by 2pi!!)
    double mig_fac,max_t_mig;   //automating length of tidal delay/migration
    //if N>3, planets repel each other and thus need to migrate longer.
    //if(_N > 2) mig_fac = 1.25 + 0.25*(_N - 2); else mig_fac=1.25;
    mig_fac=1.0;
     
    //**Initial eccentricity**
    double e=pow(mp/Ms, 0.3333333333) + efac;  //Goldreich & Schlichting (2014)
    double f=0., w=M_PI/2.;
    struct particle p = tools_init_orbit2d(Ms, mp, a, e, w, f);
    p.r = rp;
    double T=0.,t_mig_var=0.;
    a_f = a;
    assignparams(&Qp_temp,mp,rp,&T,&t_mig_var,Ms,txt_file,a,a_f,P);
    p.Qp=Qp_temp;
    particles_add(p);
    if(p_suppress == 0){
        printf("System Properties: # planets=%d, Rs=%f, Ms=%f \n",_N, Rs, Ms);
        printf("Planet 1: a=%f,P=%f,e=%f,mp=%f,rp=%f,Qp=%f,a'/a=%f,t_mig=%f \n",a,Period[0]*365./2./M_PI,e,mp,rp,Qp_temp,T,t_mig_var);
    }
    
    //deal with N>1 planets
    for(int i=1;i<_N;i++){
        extractplanets(&char_val,&a,&rho,&inc,&mp,&rp,&P,Ms);
        Period[i] = 2.*M_PI*P/365.; //in yr/2pi
        for(int k=0;k<i;k++){
            double delta = Period[i]/Period[k] - res; //calc if any "res" resonances (see var list at top)
            if(delta < RT && delta > 0.){ //current planet in res with inner one?
                double P_res = res*Period[k];
                double val = G*Ms*P_res*P_res/(4*M_PI*M_PI);
                a_f = pow(val,1./3.); //a_final of outer planet in order to be in resonance with inner
                phi_i[i] = i-k;   //how far off the resonance is (e.g. two planets away?)
                if(phi_i[i] > 1) mig_fac = 2.0; else mig_fac = 1.25; //if planet in between resonance, migrate a bit longer.
                if(p_suppress == 0)printf("%.0f:%.0f resonance for planets %i and %i, delta = %f \n",res,res-1,k+1,i+1,delta);
                break;      //can only be in a "res" resonance with one inner planet
            } else a_f = a; //if no resonance, migrate back to starting position
        }
        a *= afac;      //Increase 'a' of outer planets by afac
        e = pow(mp/Ms, 0.3333333333) + efac;
        struct particle p = tools_init_orbit2d(Ms, mp, a, e, w, f);
        p.r = rp;
        assignparams(&Qp_temp,mp,rp,&T,&t_mig_var,Ms,txt_file,a,a_f,P);
        p.Qp=Qp_temp;
        tau_a[i+1]=T;                           //migration rate
        tau_e[i+1]=T/K;                         //e_damping rate
        t_mig[i+1]=mig_fac*t_mig_var;           //length of time migrating for
        t_damp[i+1]=t_mig_var/3.;               //length of time damping migration out for
        particles_add(p);
        if(t_mig_var > t_mig[i]) max_t_mig = t_mig_var; //find max t_mig_var for tidal_delay
        mig_fac = 1.0;                          //reset
        if(p_suppress == 0) printf("Planet %i: a=%f,P=%f,e=%f,mp=%f,rp=%f,Qp=%f,a'/a=%f,t_mig=%f,t_damp=%f,afac=%f, \n",i+1,a/afac,Period[i]*365./2./M_PI,e,mp,rp,Qp_temp,tau_a[i+1],t_mig[i+1],t_damp[i+1],afac);
    }
    
    //tidal delay
    tide_delay = (mig_fac+1.0)*max_t_mig; //starts just after migration finishes
    if(tide_delay < 60000) tide_delay = 60000.;
    
	problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.
#ifndef INTEGRATOR_WH			// The WH integrator assumes a heliocentric coordinate system.
	tools_move_to_center_of_momentum();  		
#endif // INTEGRATOR_WH

}

void problem_migration_forces(){
    if(mig_forces==1){
        struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
        for(int i=1;i<N;i++){ // N = _N + 1 = total number of planets + star
                //ramp down the migration force (by increasing the migration timescale)
                double t_mig2 = t_mig[i] + t_damp[i]; //when migration ends
                if (t > t_mig[i] && t < t_mig2) {
                    tau_a[i] = tau_a[i] + (t - t_mig[i])*(3000000.0 - tau_a[i])/(t_mig2 - t_mig[i]);
                    tau_e[i] = tau_a[i]/K;
                } else if(t > t_mig2){
                    tau_a[i]=0.;
                    tau_e[i]=0.;
                }
            
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
        mu_a[i] = mu;
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

            //Stop program if nan values being produced.
            if(x_new!=x_new || y_new!=y_new || vx_new!=vx_new ||vy_new!=vy_new){
                printf("\n !!Failed run for: %s \n",c);
                printf("cartesian before: dx=%f,dy=%f,dz=%f,ex=%f,ey=%f,ez=%f,r=%f,vx=%f,vy=%f,com.vx=%f,com.vy=%f,v=%f \n",dx,dy,dz,ex,ey,ez,r,par->vx,par->vy,com.vx,com.vy,v);
                printf("Orbital elements: mu=%f,e=%f,a=%f,cosf=%.15f,sinf=%.15f,dt=%f,de=%f,da=%f,GM3a3=%f,R5a5=%f \n",mu,e,a,cosf,sinf,dt,de,da,GM3a3,R5a5);
                printf("\n cartesian after: x_new=%f,y_new=%f,vx_new=%f,vy_new=%f,term=%f,rdot=%f,rfdot=%f \n",x_new,y_new,vx_new,vy_new,term,rdot,rfdot);
                exit(0);
            }
            
            par->x = x_new;
            par->y = y_new;
            par->vx = vx_new;
            par->vy = vy_new;
            com = tools_get_center_of_mass(com,particles[i]);
            
            //print message
            if(tide_print == 0 && p_suppress == 0){
                printf("\n ***Tides have just been turned on at t=%f years***\n",t);
                tide_print = 1;
            }
            
        } else {
            n = sqrt(mu/(a*a*a));
        }     //Still need to calc this for period.
        
        int output_var=0;
        if(output_check(tmax/100000.)) output_var = 1;
        else if(t < 80000. && output_check(100.)) output_var = 1;
        if(output_var == 1){
            omega[i] = atan2(ey,ex);
            if(ey < 0.) omega[i] += 2*M_PI;
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
            double phi = 0., phi2 = 0., phi3 = 0.;     //resonant angles
            if(i>=2){//tailored for 2:1 resonance, between inner/outer planet
                phi = 2.*lambda[i] - lambda[i-phi_i[i-1]] - omega[i-phi_i[i-1]];
                phi2 = 2.*lambda[i] - lambda[i-phi_i[i-1]] - omega[i];
                phi3 = omega[i-phi_i[i-1]] - omega[i];
            }
            while(phi >= 2*M_PI) phi -= 2*M_PI;
            while(phi < 0.) phi += 2*M_PI;
            while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
            while(phi2 < 0.) phi2 += 2*M_PI;
            while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
            while(phi3 < 0.) phi3 += 2*M_PI;
            
            //output orbits in txt_file.
            FILE *append;
            append=fopen(txt_file, "a");
            //output order = time(yr/2pi),a(AU),e,P(days),arg. of peri., mean anomaly,
            //               eccentric anomaly, mean longitude, resonant angle, de/dt, 1.875/(n*mu^4/3*e)
            fprintf(append,"%e\t%.10e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,a,e,365./n,omega[i],MA,E,lambda[i],phi,phi2,phi3);
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
    free(tau_a);
    free(tau_e);
    free(lambda);
    free(omega);
    free(t_mig);
    free(t_damp);
    free(phi_i);
    free(mu_a);
}
