//For (Gold&Schlich)
term1[i] = 2*e*a*de/dt - da/dt;
coeff2[i] = 1.26*a*n;
term2a[i] = 2.38*e;

//Calculating deficit in migration due to resonance interaction (Gold&Schlich)
int GScalc = 0;
double val = 0;
double term2 = 0;
if(GScalc==1 && i>1){
    term2a[i-1] *= sin(phi);
    double term2b = 0.428*e*sin(phi2);
    coeff2[i-1] *= m/com.m;
    term2 = coeff2[i-1]*(term2a[i-1] - term2b);
    val = term1[i-1] - term2;
}

//Calculating Energy in pendulum model, j1=2,j2=-1,j4=-1 (8.6 in S.S.D.)
int ENcalc = 0;
double EN = 0,w02 = 0;
if(i>=2 && ENcalc == 1){
    double afs1 = 0.244190; //Table 8.5 S.S.D.
    double afd = -0.749964;
    double Cr = (m/com.m)*en[i-phi_i[i-1]]*afd;             //Eq. 8.32
    double Cs = (m/com.m)*en[i-phi_i[i-1]]*afs1;
    double cosphi = cos(phi);
    double wdot = 2*Cs + Cr*cosphi/e;                       //Eq. 8.30
    double epsdot = Cs*e*e + 0.5*Cr*e*cosphi;               //Eq. 8.31
    double phidot = 2*en[i] - (en[i-phi_i[i-1]] + epsdot) - wdot;
    w02 = -3*Cr*en[i-phi_i[i-1]]*e;                  //Eq. 8.47
    double sinhphi = sin(0.5*phi);
    EN = 0.5*phidot*phidot + 2*w02*sinhphi*sinhphi;  //Eq.8.48
}
