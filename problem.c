/**
Code developed by Ari Silburt to evolve Kepler planets under the influence of tides, with an 
initial migration to put planets into resonance. Updated code to work with newest version of rebound.

Some additional notes:
 - integrator_whfast_corrector = 0;  //symplectic correctors to make integration much more accurate. Slows integration considerably. Keep = 0
 - integrator_whfast_synchronize_manually = 0;	// kick/drift/kick (x and vx out of sync). Setting to 0 makes sure x, vx are synchronized at end of timestep. Setting to 1 means manual, i.e. you have to synchronize when using tides, outputting values, i.e. anything where x and vx have to be synchronized! You synchronize them by calling integrator_synchronize();
 
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
#include "integrator.h"
#include "integrator_whfast.h"
#include "../examples/tides/readplanets.h"
#include "../examples/tides/calc.h"
#include "../examples/tides/vars.h"

void problem_migration_forces();

#ifdef OPENGL
extern int display_wire;
#endif 	// OPENGL

void problem_init(int argc, char* argv[]){
    /* Setup constants */
	boxsize 	= 3;                // in AU
    integrator	= WHFAST;
    integrator_whfast_corrector = 0;
    integrator_whfast_synchronize_manually = 0;
    
    tmax        = 50000000.;  // in year/(2*pi)
    Keplername  = argv[1];          //Kepler system being investigated, Must be first string after ./nbody!
    p_suppress  = 0;                //If = 1, suppress all print statements
    double RT   = 0.06;             //Resonance Threshold - if abs(P2/2*P1 - 1) < RT, then close enough to resonance
    double timefac = 20.0;          //Number of kicks per orbital period (of closest planet)
    
    /* Migration constants */
    mig_forces  = 1;                //If ==0, no migration.
    K           = 100;              //tau_a/tau_e ratio. I.e. Lee & Peale (2002)
    e_ini       = 0.001;             //atof(argv[3]);    //initial eccentricity of the planets
    afac        = 1.01;             //Factor to increase 'a' of OUTER planets by.
    double iptmig_fac  = 1;         //reduction factor of inner planet's t_mig (lower value = more eccentricity)
    
    /* Tide constants */
    tides_on = 1;                   //If ==0, then no tidal torques on planets.
    tide_force = 0;                 //if ==1, implement tides as *forces*, not as e' and a'.
    double k2fac = 300;   //multiply k2 by this factor
    k2fac_check(Keplername,&k2fac); //For special systems, make sure that if k2fac is set too high, it's reduced.
    
#ifdef OPENGL
	display_wire 	= 1;			
#endif 	// OPENGL
	init_box();
    
    printf("%f k2fac \n \n \n",k2fac);
    
    //Naming
    naming(Keplername, txt_file, K, iptmig_fac, e_ini, k2fac, tide_force);
    
    // Initial vars
    if(p_suppress == 0) printf("You have chosen: %s \n",Keplername);
    double Ms,Rs,a,mp,rp,k2,Q,max_t_mig=0, P_temp;
    int char_val;
    
    //Star & Planet 1
    readplanets(Keplername,txt_file,&char_val,&_N,&Ms,&Rs,&mp,&rp,&P_temp,p_suppress);
    if(mig_forces == 0 && p_suppress == 0) printf("--> Migration is *off* \n");
    if(tides_on == 0 && p_suppress == 0) printf("--> Tides are *off* \n");
    struct particle star; //Star MUST be the first particle added.
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = Ms;
    star.r = Rs;
	particles_add(star);
    
    //timestep - units of 2pi*yr (required)
    dt = 2.*M_PI*P_temp/(365.*timefac);
    if(p_suppress == 0) printf("The timestep used for this simulation is (2pi*years): %f \n",dt);
    
    //Arrays, Extra slot for star, calloc sets values to 0 already.
    tau_a  = calloc(sizeof(double),_N+1);   //migration speed of semi-major axis
	tau_e  = calloc(sizeof(double),_N+1);   //migration (damp) speed of eccentricity
    lambda = calloc(sizeof(double),_N+1);   //resonant angle for each planet
    omega = calloc(sizeof(double),_N+1);    //argument of periapsis for each planet
    expmigfac = calloc(sizeof(double),_N+1);
    t_mig = calloc(sizeof(double),_N+1);
    t_damp = calloc(sizeof(double),_N+1);
    phi_i = calloc(sizeof(int),_N+1);       //phi index (for outputting resonance angles)
    if(tide_force == 1){
        //tidetau_a = calloc(sizeof(double),_N+1);
        tidetauinv_e = calloc(sizeof(double),_N+1);
    }
    double P[_N+1];       //array of period values, only needed locally
    P[1] = P_temp;
    
    //planet 1
    calcsemi(&a,Ms,P[1]);      //I don't trust archive values. Make it all consistent
    migration(Keplername,tau_a, t_mig, t_damp, &expmigfac[1], 0, &max_t_mig, P, 1, RT, Ms, mp, iptmig_fac, a, afac, p_suppress);
    struct particle p = tools_init_orbit2d(Ms, mp, a*afac, e_ini, 0, 0.);
    tau_e[1] = tau_a[1]/K;
    assignk2Q(&k2, &Q, k2fac, rp);
    p.Q=Q;
    p.k2=k2;
    p.r = rp;
    particles_add(p);
    
    //print/writing stuff
    printf("System Properties: # planets=%d, Rs=%f, Ms=%f \n",_N, Rs, Ms);
    printwrite(1,txt_file,a,P[1],e_ini,mp,rp,k2/Q,tau_a[1],t_mig[1],t_damp[1],afac,p_suppress);
    
    //outer planets (i=0 is star)
    for(int i=2;i<_N+1;i++){
        extractplanets(&char_val,&mp,&rp,&P[i],p_suppress);
        calcsemi(&a,Ms,P[i]);
        migration(Keplername,tau_a, t_mig, t_damp, &expmigfac[i], &phi_i[i], &max_t_mig, P, i, RT, Ms, mp, iptmig_fac, a, afac, p_suppress);
        struct particle p = tools_init_orbit2d(Ms, mp, a*afac, e_ini, 0, i*M_PI/4.);
        tau_e[i] = tau_a[i]/K;
        assignk2Q(&k2, &Q, k2fac, rp);
        p.Q = Q;
        p.k2 = k2;
        p.r = rp;
        particles_add(p);
        printwrite(i,txt_file,a,P[i],e_ini,mp,rp,k2/Q,tau_a[i],t_mig[i],t_damp[i],afac,p_suppress);
    }
    
    //tidal delay
    if(max_t_mig < 20000)tide_delay = 20000.; else tide_delay = max_t_mig + 20000.;    //Have at least 30,000 years grace before turning on tides.
    double tide_delay_output = 0;
    if(tides_on == 1) tide_delay_output = tide_delay;
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%f \n",tide_delay_output);
    fclose(write);
    tide_print = 0;
    
    problem_additional_forces = problem_migration_forces; 	//Set function pointer to add dissipative forces.
    integrator_force_is_velocitydependent = 1;
    if (integrator != WH){	// The WH integrator assumes a heliocentric coordinate system.
        tools_move_to_center_of_momentum();
    }

}

void problem_migration_forces(){
    if(mig_forces==1){
        int print = 1;
        struct particle com = particles[0]; // calculate migration forces with respect to center of mass;
        for(int i=1;i<N;i++){ // N = _N + 1 = total number of planets + star
                double t_migend = t_mig[i] + t_damp[i]; //when migration ends
                if (t > t_mig[i] && t < t_migend) { //ramp down the migration force (by increasing the migration timescale)
                    tau_a[i] *= expf(dt/expmigfac[i]);
                    tau_e[i] = tau_a[i]/K;
                } else if(t > t_migend){
                    tau_a[i]=0.;
                    tau_e[i]=0.;
                    double sum = 0;
                    for(int j=0;j<N;j++) sum += tau_a[j];
                    if(sum < 0.1){
                        mig_forces = 0; //turn migration loop off altogether, save calcs
                        if(p_suppress == 0 && print == 1) printf("\n\n **migration loop off at t=%f** \n\n",t);
                        print = 0;
                    }
                }
            
            if (tau_e[i]!=0||tau_a[i]!=0){
                struct particle* p = &(particles[i]);
                const double dvx = p->vx-com.vx;
                const double dvy = p->vy-com.vy;
                const double dvz = p->vz-com.vz;
                
                if (tau_a[i]!=0){ 	// Migration
                    double const term = 1/(2.*tau_a[i]);
                    p->ax -=  dvx*term;
                    p->ay -=  dvy*term;
                    p->az -=  dvz*term;
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
                    const double a = -mu/( v*v - 2.*mu/r );			// semi major axis, AU
                    
                    //TESTP5m*********************** to get them all starting off in line together
                    /**
                    if(a < 0.091432 && t > 500 && i==2){
                        mig_forces = 0;
                        if(p_suppress == 0) printf("\n\n **migration loop off (abrupt) at t=%f** \n\n",t);
                    }
                    **/
                    
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
    
    //This calculates the tau_e for tides as forces, now that we have the resonance eccentricities
    if(tautide_force_calc == 0 && mig_forces == 0. && tide_force == 1 && tides_on == 1){
        tautide_force_calc = 1; //only do this process once
        struct particle com = particles[0];
        for(int i=1;i<N;i++){
            struct particle* par = &(particles[i]);
            const double m = par->m;
            const double mu = G*(com.m + m);
            const double rp = par->r*0.00464913;       //Rp from Solar Radii to AU
            const double Qp = par->k2/(par->Q);
            
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
            const double e = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
            const double a = -mu/( v*v - 2.*mu/r );			// semi major axis, AU
            double a5r5 = pow(a/rp, 5);
            
            tidetauinv_e[i] = 1.0/(2./(9*M_PI)*(1./Qp)*sqrt(a*a*a/com.m/com.m/com.m)*a5r5*m);
            //tidetau_a[i] = tidetau_e[i]*K; //Dan uses a K factor instead.
            
            if(p_suppress == 0) printf("planet %i: tau_e,=%.1f Myr, a=%f,e=%f \n",i,1.0/(tidetauinv_e[i]*1e6),a,e);
        }
    }
    
    if(tides_on == 1 && tide_force == 1 && t > tide_delay){
        struct particle com = particles[0]; // calculate add. forces w.r.t. center of mass
        for(int i=1;i<N;i++){
            struct particle* p = &(particles[i]);
            const double dvx = p->vx - com.vx;
            const double dvy = p->vy - com.vy;
            const double dvz = p->vz - com.vz;
            
          //if(i==1){
                /*Papaloizou & Larwood (2000) - Unlike the orbital elements version of tides,
                 tidetau_e already induces an a' (Gold&Schlich2015), and so the tidetau_a is
                 an additional migration on top... don't think I need it. */
                /*
                if (tidetau_a[i] != 0.){
                    p->ax +=  -dvx/(2.*tidetau_a[i]);
                    p->ay +=  -dvy/(2.*tidetau_a[i]);
                    p->az +=  -dvz/(2.*tidetau_a[i]);
                }
                 */
                 
                //Papaloizou & Larwood (2000)
                if (tidetauinv_e[i] != 0. ){ 	// need h and e vectors for both types
                    const double dx = p->x-com.x;
                    const double dy = p->y-com.y;
                    const double dz = p->z-com.z;
                
                    const double rinv = 1/sqrt ( dx*dx + dy*dy + dz*dz );
                    const double vr = (dx*dvx + dy*dvy + dz*dvz)*rinv;
                    const double term = -2.*tidetauinv_e[i]*vr*rinv;
                    
                    p->ax += term*dx;
                    p->ay += term*dy;
                    p->az += term*dz;

                }
                com = tools_get_center_of_mass(com,particles[i]);
            //}
        }
        
        //print message
        if(tide_print == 0 && p_suppress == 0){
            printf("\n\n ***Tides (forces!) have just been turned on at t=%f years***\n\n",t);
            tide_print = 1;
        }
    }
}

void problem_inloop(){
}

void problem_output(){
    //conditions for entering loops
    int output_var=0;
    if(output_check(tmax/100000.)) output_var = 1; //Used to be 100,000
    else if(t < tide_delay && output_check(100.)) output_var = 1; //used to be 100
    int tide_go = 0;
    if(tides_on == 1 && tide_force == 0 && t > tide_delay) tide_go = 1;
    
    //**Main Loop**
    if(output_var == 1 || tide_go == 1){
        //Calculate Orbital Elements
        struct particle com = particles[0];
        for(int i=1;i<N;i++){
            struct particle* par = &(particles[i]);
            const double m = par->m;
            const double mu = G*(com.m + m);
            const double rp = par->r*0.00464913;       //Rp from Solar Radii to AU, G=1, [t]=yr/2pi, [m]=m_star
            const double Qp = par->k2/(par->Q);
            
            const double dvx = par->vx-com.vx;
            const double dvy = par->vy-com.vy;
            const double dvz = par->vz-com.vz;
            const double dx = par->x-com.x;
            const double dy = par->y-com.y;
            const double dz = par->z-com.z;
            
            const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
            const double r = sqrt ( dx*dx + dy*dy + dz*dz );
            const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
            const double rinv = 1/r;            //some extra terms to speed up code
            const double muinv = 1./mu;
            const double term1 = v*v-mu*rinv;
            const double term2 = r*vr;
            const double ex = muinv*( term1*dx - term2*dvx );
            const double ey = muinv*( term1*dy - term2*dvy );
            const double ez = muinv*( term1*dz - term2*dvz );
            double e = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
            
            // true anomaly + periapse (wiki, Fund. of Astrodyn. and App., by Vallado, 2007)
            const double rdote = dx*ex + dy*ey + dz*ez;
            double cosf = rdote/(e*r);
            if(cosf >= 1.) cosf = 1.;
            if(cosf <= -1.) cosf = -1.;
            double sinf = sqrt(1. - cosf*cosf);
            if(vr < 0.) sinf *= -1.;
            double sinwf = dy*rinv;
            double coswf = dx*rinv;
            double a = r*(1. + e*cosf)/(1. - e*e);
            double n;
            
            //Test for collision
            if(N < _N+1 && collision_print == 0){
                printf("\n\n system %s with e_ini=%f,e_now=%f had a collision!! \n\n",Keplername,e_ini,e);
                FILE *append;
                append=fopen("python_scripts/Kepler_e_coll.txt", "a");
                fprintf(append,"%s,%e,%e\n",Keplername,e_ini,t);
                fclose(append);
                collision_print = 1;
            }
            
            //Tides
            if(tide_go == 1){//For TESTP5m need && i==1
                const double a2 = a*a;
                const double rp2 = rp*rp;
                const double R5a5 = rp2*rp2*rp/(a2*a2*a);
                const double GM3a3 = sqrt(G*com.m*com.m*com.m/(a2*a));
                const double de = -dt*(9.*M_PI*0.5)*Qp*GM3a3*R5a5*e/m;   //Tidal change for e
                const double da = 2.*a*e*de;                             //Tidal change for a
                
                a += da;
                e += de;
                
                integrator_whfast_particles_modified = 1;
                
                /*
                //Test repulsion vs. tugging (change w to put planets out of res)
                if(repuls_v_tugg_on == 0 && i==1 && p_suppress == 0){
                    repuls_v_tugg_on = 1;
                    double angle_change = 0.34*(M_PI/180.);
                    double costemp = coswf;
                    double sintemp = sinwf;
                    coswf = costemp*cos(angle_change) - sintemp*sin(angle_change);
                    sinwf = sintemp*cos(angle_change) - costemp*sin(angle_change);
                    printf("\n  ***Testing Repulsion vs. tugging. Planets put out of resonance.*** \n");
                }*/
                
                //Re-update coords.
                const double r_new = a*(1. - e*e)/(1. + e*cosf);
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
                    printf("\n !!Failed run for: %s \n",Keplername);
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
                    printf("\n\n ***Tides (a', e') have just been turned on at t=%f years***\n\n",t);
                    tide_print = 1;
                }
                
            } else {
                n = sqrt(mu/(a*a*a)); //Still need to calc this for period.
            }
            
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
                if(i>1){//tailored for 2:1 resonance, between inner/outer planet
                    phi = 2.*lambda[i] - lambda[i-phi_i[i]] - omega[i-phi_i[i]];
                    phi2 = 2.*lambda[i] - lambda[i-phi_i[i]] - omega[i];
                    phi3 = omega[i-phi_i[i]] - omega[i];
                }
                while(phi >= 2*M_PI) phi -= 2*M_PI;
                while(phi < 0.) phi += 2*M_PI;
                while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
                while(phi2 < 0.) phi2 += 2*M_PI;
                while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
                while(phi3 < 0.) phi3 += 2*M_PI;
                
                //integrator_synchronize(); //if synchronize_manual = 1, then call this before each output to synchronize x and vx.
                
                //output orbits in txt_file.
                FILE *append;
                append=fopen(txt_file, "a");
                //output order = time(yr/2pi),a(AU),e,P(days),arg. of peri., mean anomaly,
                //               eccentric anomaly, mean longitude, resonant angle, de/dt, phi1     phi2     phi3
                fprintf(append,"%e\t%.10e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,a,e,365./n,omega[i],MA,E,lambda[i],phi,phi2,phi3);
                fclose(append);
                
                if (integrator != WH){	// The WH integrator assumes a heliocentric coordinate system.
                    tools_move_to_center_of_momentum();
                }
            }
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
    if(tide_force == 1){
        //free(tidetau_a);
        free(tidetauinv_e);
    }
}
