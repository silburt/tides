#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stddef.h>
#include<unistd.h>
#include "readplanets.h"
#include "../../src/main.h"

void readplanets(char* sysname, char* txt_file, int* char_pos, int* _N, double* Ms, double* Rs, double* mp, double* rp, double* P, int p_suppress){
    FILE *f = fopen("planets.txt", "r");
    char temp[512];
    int line_num = 0, found_result=0, exit=0;
    
    while(exit != 1){
        line_num += 1;
        fgets(temp, 512, f);      //get row of data from planets.txt
        if((strstr(temp, sysname)) != NULL){ //see if matches Kepler system name.
            *char_pos=ftell(f);
            if(p_suppress == 0) printf("A match found on line: %d, proceed with sim. \n\n", line_num);
            //printf("char_pos=%i",*char_pos);
            //printf("\n%s\n", temp);
            found_result++;
            exit = 1;
        }
    }
    if(found_result == 0) printf("\n Sorry, couldn't find a match.\n");
    if(f) fclose(f);
    
    //split temp into tokens using ',' delimeter, store in arrays
    //pointed to by datpoint.
    const int numfields=25;
    int i=0;
    char *string = temp;
    double array[numfields-4];
    while (i<numfields){
        char *q = strsep(&string,",");
        if(i>=4){
            array[i-4] = (atof(q));
        }
        i++;
    }
    
    *_N = array[0];                     //number of planets in system
    //*a = array[4];                      //semi-major axis (AU)
    *Ms = array[11];                    //Stellar mass
    *Rs = array[13];                    //Stellar radius
    *mp = array[15]*3e-6;               //planet mass (SOLAR units)
    *rp = array[18];                    //planet radius (SOLAR units)
    *P = array[1];                      //Period (days)
    
    if(*Ms == 0.){
        *Ms = pow(*Rs,1.00);
        if(p_suppress == 0) printf("--> calculated stellar mass \n");
    }
    
    double solar2earthRp = 109.21;
    double earth2solarMp = 3e-6;
    if(*mp == 0. && *rp < 0.04){//Weiss & Marcy 2014
        *mp = 2.69*pow(*rp*solar2earthRp,0.93)*earth2solarMp; //Solar mass units
        if(p_suppress == 0) printf("--> Planet 1 - Calculated planet mass (Earth-realm) \n");
    } else if(*mp == 0. && *rp >=0.04){//Jupiter scaling relation
        double r3 = pow(*rp*695800000.,3); //in meters
        *mp = 1330*r3/2e30;     //in solar mass (density = 1330 kg/m^3 * 4/3*pi = 5.554)
        if(p_suppress == 0) printf("--> Planet 1 - Calculated planet mass (Jovian-realm) \n");
    }
    
    //delete previous output file
    char sys_arg[100] = "rm -v ";
    strcat(sys_arg,txt_file);
    system(sys_arg);
    //write star characteristics to file. Planet characteristics come in assignparams.c
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%s,%f,%f,%i \n",sysname,*Ms,*Rs,*_N);
    fclose(write);
}

//*******************************************************************************//

void extractplanets(int* char_pos, double* mp, double* rp, double* P, int p_suppress){
    FILE *f = fopen("planets.txt", "r");
    char temp[512];
    fseek(f, *char_pos, SEEK_SET);
    fgets(temp, 512, f);
    *char_pos = ftell(f);
    
    if(f) {
        fclose(f);
    }
    //split temp into tokens using ',' delimeter, store in arrays
    //pointed to by datpoint.
    const int numfields=25;
    int i=0;
    char *string = temp;
    double array[numfields-4];
    while (i<numfields){
        char *p = strsep(&string,",");
        if(i>=4){
            array[i-4] = (atof(p));
            //printf("string=%f,%i \n",array[i-4],i);
        }
        i++;
    }
    //*a = array[4];          //semi-major axis (AU)
    *mp = array[15]*3e-6;   //planet mass (SOLAR units)
    *rp = array[18];        //planet radius (SOLAR units)
    *P = array[1];          //Period (days)
    
    double solar2earthRp = 109.21;
    double earth2solarMp = 3e-6;
    if(*mp == 0. && *rp < 0.04){//Weiss & Marcy 2014, Neptune-sized
        *mp = 2.69*pow(*rp*solar2earthRp,0.93)*earth2solarMp; //Solar mass units
        if(p_suppress == 0) printf("--> Calculated planet mass (Earth-realm) \n");
    } else if(*mp == 0. && *rp >=0.04){//Jupiter density scaling relation
        double r3 = pow(*rp*695800000.,3); //in meters
        *mp = 1330*r3/2e30;     //in solar mass
        if(p_suppress == 0) printf("--> Calculated planet mass (Jovian-realm) \n");
    }
    
}

void naming(char* sysname, char* txt, double K, double iptmig_fac, double e_ini, double Qpfac, int tide_force){
    char* dir = "runs/orbits_";
    char* ext = ".txt";
    strcat(txt, dir);
    strcat(txt, sysname);
    char* str = "_Qpfac";
    strcat(txt, str);
    char strQpfac[15];
    int Qpfactor = (int) Qpfac;
    sprintf(strQpfac, "%d", Qpfactor);
    strcat(txt, strQpfac);
    if(tide_force == 1){
        char* forcestring = "_tideF";   //implementing tides as forces
        strcat(txt, forcestring);
    }
    if(K != 100){
        char strK[15];
        int Kint = (int) K;
        sprintf(strK, "%d", Kint);
        strcat(txt, "_K");
        strcat(txt, strK);
    }
    if(iptmig_fac != 1){
        char strmig[15];
        int migint = (int) round(10*iptmig_fac);
        sprintf(strmig, "%d", migint);
        strcat(txt, "_migfac0.");
        strcat(txt, strmig);
    }
    if(e_ini != 0.01){
        char strmig[15];
        int migint = (int) round(100*e_ini);
        printf("migint = %i \n",migint);
        sprintf(strmig, "%d", migint);
        strcat(txt, "_ei0.");
        strcat(txt, strmig);
    }
    strcat(txt, ext);
}

void printwrite(int i, char* txt_file, double a,double P,double e,double mp,double rp,double Qp,double tau_a,double t_mig,double t_damp,double afac,int p_suppress){
    
    if(p_suppress == 0) printf("Planet %i: a=%f,P=%f,e=%f,mp=%f,rp=%f,Qp=%f,a'/a=%f,t_mig=%f,t_damp=%f,afac=%f, \n",i,a,P,e,mp,rp,Qp,tau_a,t_mig,t_damp,afac);
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%.10f,%f,%f,%f,%f,%f\n", mp,rp,P,Qp,tau_a,t_mig+t_damp);
    fclose(write);
}
