#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<math.h>
#include<stddef.h>
#include<time.h>
#include "assignparams.h"
#include "../../src/main.h"

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

void assignparams(double* Qp, double Qpfac, double mp, double rp, double* T, double* t_mig_var, double Ms, char* txt_file, double a, double a_f, double P, double migspeed_fac){
    
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
    
    //Assign tau_migration values
    double n = 365.*2*M_PI/P;  //units = 2Pi/yr
    double mu43 = pow(mp/Ms,4./3.);
    /*Goldreich & Schlichting (2014), tau_mig rate for 2:1 resonance, units = yr/2pi. Min is 3.75, 
     but use 5.0 to be safe.*/
    *T = 5.00*migspeed_fac/(n*mu43);
    //*T = 10.00/(n*mu43);
    *t_mig_var = *T*(a - a_f)/a_f;  //length of time migrate for, units = yr/2pi
    if(a == a_f){//first planet
        *T = 0.;
        *t_mig_var = 0.;
    }
    
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%.10f,%f,%f,%f,%f,%f\n", mp,rp,P,*Qp,*T,*t_mig_var);
    fclose(write);
}

void special_cases(char* sysname, int i, double* mig_fac){
    if(strcmp(sysname, "Kepler-32") == 0 && i == 2) *mig_fac = 1.40;
    if(strcmp(sysname, "Kepler-11") == 0 && i == 4){ *mig_fac = 3.0; printf("mig_fac=%f \n",*mig_fac);}
}

//Calculate tau_a, tau_e for Paploizou & Larwood (2000) version of tides.
void calc_tidetau(double* tau_a, double* tau_e, double Qp, double mp, double rp, double Ms, double e_default, double a_default, char* sysname, int i, int p_suppress){
    
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
    *tau_a = *tau_e/(2*e*e);
    printf("tau_e,tau_a,e = %f,%f,%f \n",*tau_e,*tau_a,e);
    
}