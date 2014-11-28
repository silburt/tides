#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stddef.h>

void readplanets(char *sysname, int *char_pos, int *_N, double *Ms, double *Rs, double *a, double *rho, double *inc, double *mp, double *rp);

void extractplanets(int *char_pos, double *a, double *rho, double *inc, double *mp, double *rp);


void readplanets(char *sysname, int *char_pos, int *_N, double *Ms, double *Rs, double *a, double *rho, double *inc, double *mp, double *rp){
    FILE *f = fopen("planets.csv", "r");
    char temp[512];
    int line_num = 0, found_result=0, exit=0;
    
    while(exit != 1){
        line_num += 1;
        fgets(temp, 512, f);
        if((strstr(temp, sysname)) != NULL){
            *char_pos=ftell(f);
            //printf("A match found on line: %d\n", line_num);
            //printf("char_pos=%i",*char_pos);
            //printf("\n%s\n", temp);
            found_result++;
            exit = 1;
        }
        
    }
    
    if(found_result == 0) {
        printf("\n Sorry, couldn't find a match.\n");
    }
    if(f) {
        fclose(f);
    }
    
    //split temp into tokens using ',' delimeter, store in arrays
    //pointed to by datpoint.
    const int numfields=25;
    int i=0;
    double array[numfields-4];
    char *p;
    p=strtok(temp,",");
    while (i<numfields){
        if(i>=4){
            array[i-4] = (atof(p));
            //printf("string=%f,%i \n",array[i-4],i);
        }
        p = strtok (NULL, ",");
        i++;
    }
    
    *_N = array[0];         //number of planets in system
    *a = array[4];          //semi-major axis (AU)
    *rho = array[7];        //density (g/cm**3)
    *inc = array[8];        //inclination
    *Ms = array[11];        //Stellar mass
    *Rs = array[13];        //Stellar radius
    *mp = array[15]*1e-6;   //planet mass (SOLAR units)
    *rp = array[18];        //planet radius (SOLAR units)
    
    if(*a==0. && array[11] != 0.){//many semi-major axis fields are empty. Calc
        double P = array[1]*24.*60.*60.; //Period in seconds
        double mass = array[11]*1.989e30;
        double G_SI = 6.67e-11;
        double calca = P*P*G_SI*mass/(4*M_PI*M_PI);
        *a = pow(calca,1./3.)/1.496e11;     //in AU
    }
}


void extractplanets(int *char_pos, double *a, double *rho, double *inc, double *mp, double *rp){
    FILE *f = fopen("planets.csv", "r");
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
    //char *array[numfields]; //store each field in its own slot
    double array[numfields-4];
    while (i<numfields){
        char *p = strsep(&string,",");
        //printf("p=%s\n",p);
        if(i>=4){
            array[i-4] = (atof(p));
            //printf("string=%f,%i \n",array[i-4],i);
        }
        //p = strsep (NULL, ",");
        i++;
    }
    *a = array[4];          //semi-major axis (AU)
    *rho = array[7];        //density (g/cm**3)
    *inc = array[8];        //inclination
    *mp = array[15]*1e-6;   //planet mass (SOLAR units)
    *rp = array[18];        //planet radius (SOLAR units)
    
    if(*a==0. && array[11] != 0.){
        double P = array[1]*24.*60.*60.; //Period in seconds
        double mass = array[11]*1.989e30;
        double G_SI = 6.67e-11;
        double calca = P*P*G_SI*mass/(4*M_PI*M_PI);
        *a = pow(calca,1./3.)/1.496e11;     //in AU
    }
    
}
