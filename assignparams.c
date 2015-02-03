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

void assignparams(double* Qp, double mp, double rp, double* T, double* t_mig_var, double Ms, char* txt_file, double a, double a_f, double P){
    //Assign k2/Q
    double k2, Q;
    k2 = 0.1;
    //k2 = 1.;
    
    if(mp > 1e-4 && mp < 1e-3){//Uranus/Neptune Q
        //Q = (rand() %25 + 60)*1e3;
        Q = 7.2e4;
    } else if(mp >= 1e-3){//Saturn/Jupiter Q
        //Q = (rand() %30 + 5)*1e4;
        Q = 1e5;
    } else{//Earth or smaller Q
        //Q = rand() %80 + 10;
        Q = 10.;
    }
    *Qp = k2/Q;
    
    //Assign tau_migration values
    double n = 365.*2*M_PI/P;  //units = 2Pi/yr
    double mu43 = pow(mp/Ms,4./3.);
    //Goldreich & Schlichting (2014), tau_mig rate for 2:1 resonance, units = yr/2pi. Min is 3.75, but use 4.0 to be safe
    *T = 4.00/(n*mu43);
    *t_mig_var = *T*(a - a_f)/a_f;  //length of time migrate for, units = yr/2pi
    if(a == a_f){//first planet
        *T = 0.;
        *t_mig_var = 0.;
    }
    
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%f,%f,%f,%f,%f,%f\n", mp,rp,P,*Qp,*T,*t_mig_var);
    fclose(write);
}
