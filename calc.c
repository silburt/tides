//
//  migration.c
//  
//
//  Created by Ari Silburt on 2015-04-14.
//
//
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<math.h>
#include<stddef.h>
#include<time.h>
#include "calc.h"
#include "readplanets.h"
#include "../../src/main.h"

void migration(char* sysname, double* tau_a, double* t_mig, double* t_damp, double *expmigfac, int* phi_i, double* max_t_mig, double* P, int i, double RT, double Ms, double mp, double iptmig_fac, double a, double afac, int p_suppress){
    //Resonance vars
    double mig_fac = 1.0;   //automating length of tidal delay/migration
    double Pfac = 2.*M_PI/365.; //converts period to yr/2pi
    double mintau_a = 5000.; //minimum a'/a
    double mintau_fac = 6.0; //absolute min is 3.75 (Gold&Schlich), but use 6 to be safe
    double min_rel_speed = 0.75;    //was 0.6
    double a_f = a; //default - if no resonance, migrate back to starting position
    int special_flag = 0; //for some special cases, do not reduce t_mig of inner planet
    
    //Goldreich & Schlichting (2014), mig rate for 2:1 resonance, units = yr/2pi.
    double n = 365.*2*M_PI/P[i];  //units = 2Pi/yr
    double mu43 = pow(mp/Ms,4./3.);
    tau_a[i] = mintau_fac/(n*mu43);
    if(tau_a[i] < mintau_a) tau_a[i] = mintau_a;
    
    for(int k=1;k<i;k++){
        double delta = P[i]/P[k] - 2.0; //calc if any 2:1 resonances
        if(delta < RT && delta > 0.){ //check for res with inner plannet
            double P_res = 2.0*Pfac*P[k];    //2*period of inner planet (units of yr/2pi)
            double val = G*Ms*P_res*P_res/(4*M_PI*M_PI);
            a_f = pow(val,1./3.); //a_final of outer planet in order to be in resonance with inner
            *phi_i = i-k;   //how far off the resonance is (e.g. two planets away?)
            if(*phi_i > 1) mig_fac = 1.75; else mig_fac = 1.2; //if planet in between resonance, migrate a bit longer.
            
            //double rel_speed = 0.75;    //*relative* migration velocity (key is relative)
            //if(rel_speed*tau_a[k]/(1. - rel_speed) > 5.0/(n*mu43)){ //condition for certain capture
            double rel_speed = 1/((n*mu43*tau_a[k])/mintau_fac + 1.); //fastest mig speed with guaranteed capture
            if(rel_speed < min_rel_speed) rel_speed = min_rel_speed;    //don't want it to be too fast.
            tau_a[i] = rel_speed*(tau_a[k]);    //set outer migration rate to rel_speed*tau_a[k]
            special_cases(sysname,i,k,&special_flag);
            if(special_flag == 0){
                double red_fac = 0.65;
                t_mig[k] *= red_fac*iptmig_fac;  //inner planet migrates for much less time
                t_damp[k] *= red_fac*iptmig_fac;
            }
            printf("** a/a' (outer) = %f a/a' (inner) ** (guarantees migration whilst in resonance) \n",rel_speed);
            if(p_suppress == 0)printf("2:1 resonance for planets %i and %i, delta = %f \n",k,i,delta);
            break;      //can only be in a "res" resonance with one inner planet
        }
    }
    
    //migration timescale
    t_mig[i] = tau_a[i]*mig_fac*(a*afac - a_f)/a_f;  //length of time migrate for, units = yr/2pi
    if(t_mig[i] + t_damp[i] > *max_t_mig) *max_t_mig = t_mig[i] + t_damp[i]; //find max t_mig_var for tidal_delay
    
    //migration damping timescale - need min damp time or weird eccentricity effects ensue
    double damp_fac = 3.0;
    t_damp[i] = t_mig[i]/damp_fac; //Need to damp minimum over a libration timescale.
    *expmigfac = t_damp[i]/log(2000000./tau_a[i]);
    
    //The amount of distance covered from the exp damp decay is equivalent to t_equiv travelling at tau_a[i].
    //Since we want inner planet to end up at its initial position, need to subtract this from mig_fac
    /*
    if(i > 1){
        int kk2 = i - *phi_i;
        double t_equiv = *expmigfac*(1 - exp(-t_damp[kk2]/ *expmigfac));
        t_mig[kk2] -= t_equiv;
        if(t_mig[kk2] < 0) t_mig[kk2] = 0;
    }
    */
}

void assignQp(double* Qp, double Qpfac, double rp){
    /*
     Q (Goldreich & Soter, 1966):
     Mercury <  190
     Venus <    17
     Earth =    13 (MacDonald, 1964)
     Mars <     26
     Jupiter ~  (1 to 2)e5
     Saturn ~   (6 to 7)e4
     Uranus >   7.2e4
     Neptune >  7.2e4
     
     Q & k_2 (Yoder, 1995 --> But actually from Solar System Dynamics)
     Values in brackets have been estimated (see pg. 166, Ch. 4.10)
     Mercury=   (100),  (0.1), apparently 0.5 (Padovan et al.)
     Venus  =   (100),  0.25
     Earth  =   12,     0.299
     Mars   =   86      0.14
     
     k_2 of Giant planets (Gavrilov & Zharkov, 1977)
     Jupiter=   0.379
     Saturn =   0.341
     Uranus =   0.104
     Neptune=   0.127
     */
    
    //Assign Qp = k2/Q
    if(rp > 2*0.009156 && rp < 0.1){//Lee/Fabrycky/Lin distringuish Rp > 2R_E as mini-Neptune vs. Super-Earth
        *Qp = Qpfac*1./(2.2e4); //Lowest Neptune value (Lee/Fabrycky/Lin 2013)
    } else if(rp >= 0.1){
        *Qp = Qpfac*1./(5.4e4); //Lowest Saturn value (Lee/Fabrycky/Lin 2013)
    } else if(rp <0.005){
        *Qp = Qpfac*10000.;
    } else{
        *Qp = Qpfac*1./40.;   //lowest Earth value (Lee/Fabrycky/Lin 2013)
    }
}

//Calculate tidal tau_a=a/a', tau_e=e/e' for Paploizou & Larwood (2000) version of tides.
void calc_tidetau(double* tau_a, double* tau_e, double K, double Qp, double mp, double rp, double Ms, double e_default, double a_default, char* sysname, int i, int p_suppress){
    
    FILE *f = fopen("reso/Kepler_ei.txt", "r");
    char temp[512];
    int line_num = 0, found_result=0, exit=0, char_pos=0;
    
    while(exit != 1){
        line_num += 1;
        fgets(temp, 512, f);      //get row of data from planets.txt
        if((strstr(temp, sysname)) != NULL){ //see if matches Kepler system name.
            char_pos=ftell(f);
            if(p_suppress == 0 && i==0) printf("Calc tidetau_e & tidetau_a. \n");
            found_result++;
            exit = 1;
        }
        if(line_num > 100) exit = 1;
    }
    if(f) fclose(f);
    
    
    double e,a;      //eccentricity once in resonance
    if(found_result == 0){
        printf("calc_tidetau: Cannot find %s! guestimate e_in=e_default, e_o=0.0925 \n",sysname);
        if(i==0){e = e_default;} else{ e=0.0925;}
        a = a_default;
    } else {
        const int numfields=7;
        int j=0;
        char *string = temp;
        double array[numfields];
        while (j<numfields){
            char *q = strsep(&string,",");
            array[j] = (atof(q));
            j++;
        }
        
        if(i == array[5]){ e = array[1]; a = array[2]; }    //inner res planet
        else if(i == array[6]){ e = array[3]; a = array[4]; }   //outer res planet
        else { e = e_default; a = a_default; }  //not a res planet
    }
    double a5r5 = pow(a/(rp*0.00464913), 5);
    *tau_e = 2./(9*M_PI)*(1./Qp)*sqrt(a*a*a/Ms/Ms/Ms)*a5r5*mp;
    //*tau_a = *tau_e/(2*e*e);
    *tau_a = *tau_e*K;      //Dan uses a K factor instead.
    printf("tau_e,tau_a,e = %f,%f,%f \n",*tau_e,*tau_a,e);
    
}

void printwrite(int i, char* txt_file, double a,double P,double e,double mp,double rp,double Qp,double tau_a,double t_mig,double t_damp,double afac,int p_suppress){
    
    if(p_suppress == 0) printf("Planet %i: a=%f,P=%f,e=%f,mp=%f,rp=%f,Qp=%f,a'/a=%f,t_mig=%f,t_damp=%f,afac=%f, \n",i,a,P,e,mp,rp,Qp,tau_a,t_mig,t_damp,afac);
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%.10f,%f,%f,%f,%f,%f\n", mp,rp,P,Qp,tau_a,t_mig+t_damp);
    fclose(write);
}

void special_cases(char* sysname, int i, int k, int* special_flag){
    if(strcmp(sysname, "Kepler-31") == 0 && i == 3 && k == 2){
        *special_flag = 1;
        printf("Special_flag=1 \n");
    }
}