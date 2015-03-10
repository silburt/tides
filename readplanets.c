#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stddef.h>
#include<unistd.h>
#include "readplanets.h"
#include "../../src/main.h"

void readplanets(char* sysname, char* txt_file, int* char_pos, int* _N, double* Ms, double* Rs, double* rho, double* inc, double* mp, double* rp, double* P, double* dt, double timefac, int p_suppress){
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
    *rho = array[7];                    //density (g/cm**3)
    *inc = array[8];                    //inclination
    *Ms = array[11];                    //Stellar mass
    *Rs = array[13];                    //Stellar radius
    *mp = array[15]*3e-6;               //planet mass (SOLAR units)
    *rp = array[18];                    //planet radius (SOLAR units)
    *P = array[1];                      //Period (days)

    //timestep (P/11.), in year/(2*pi).
    //From Viswanath, Divakar 2002-03, need min. of 6*dt per P
    //*dt = 2.*M_PI*array[1]/(365.*11.);
    *dt = 2.*M_PI*array[1]/(365.*timefac);
    if(p_suppress == 0) printf("The timestep used for this simulation is (years/2pi): %f \n",*dt);
    
    if(*Ms == 0.){
        *Ms = pow(*Rs,1.25);
        if(p_suppress == 0) printf("--> calculated stellar mass \n");
    }
    
    double solar2earthRp = 109.21;
    double earth2solarMp = 3e-6;
    if(*mp == 0. && *rp < 0.04){//Weiss & Marcy 2014
        *mp = 2.69*pow(*rp*solar2earthRp,0.93)*earth2solarMp; //Solar mass units
        if(p_suppress == 0) printf("--> Planet 1 - Calculated planet mass (Earth-realm) \n");
    } else if(*mp == 0. && *rp >=0.04){//Jupiter scaling relation
        *mp = 0.001*(*rp/0.1);
        if(p_suppress == 0) printf("--> Planet 1 - Calculated planet mass (Jovian-realm) \n");
    }
    
    //delete previous output file
    char sys_arg[50] = "rm -v ";
    strcat(sys_arg,txt_file);
    system(sys_arg);
    //write star characteristics to file. Planet characteristics come in assignparams.c
    FILE *write;
    write=fopen(txt_file, "a");
    fprintf(write, "%s,%f,%f,%i,",sysname,*Ms,*Rs,*_N);
    fclose(write);
}

//*******************************************************************************//

void extractplanets(int* char_pos, double* rho, double* inc, double* mp, double* rp, double* P, int p_suppress){
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
    *rho = array[7];        //density (g/cm**3)
    *inc = array[8];        //inclination
    *mp = array[15]*3e-6;   //planet mass (SOLAR units)
    *rp = array[18];        //planet radius (SOLAR units)
    *P = array[1];          //Period (days)
    
    double solar2earthRp = 109.21;
    double earth2solarMp = 3e-6;
    if(*mp == 0. && *rp < 0.04){//Weiss & Marcy 2014, Neptune-sized
        *mp = 2.69*pow(*rp*solar2earthRp,0.93)*earth2solarMp; //Solar mass units
        if(p_suppress == 0) printf("--> Calculated planet mass (Earth-realm) \n");
    } else if(*mp == 0. && *rp >=0.04){//Jupiter scaling relation
        *mp = 0.001*(*rp/0.1);
        if(p_suppress == 0) printf("--> Calculated planet mass (Jovian-realm) \n");
    }
    
}

void calcsemi(double* a, double Ms, double P){
    double P_SI = P*24.*60.*60.; //Period in seconds
    double mass = Ms*1.989e30;
    double G_SI = 6.67e-11;
    double a3 = P_SI*P_SI*G_SI*mass/(4.*M_PI*M_PI);
    *a = pow(a3,1./3.)/1.496e11;     //in AU
}
