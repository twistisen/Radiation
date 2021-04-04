#define PI 3.14159265358979323846264338327950288419716939937510582
#include <cmath>
#define THREADS 256
#define BLOCKS 15360
#define thetaxfine 128
#define thetayfine 128
#define nfreq 240
#define delta_n 5
#define memsize 40

//#define N 150
// constants

const double pi = PI;



//const double m=510.998928*1e3;
const double meter = 1.97326979e-7;
const double re = 2.8179403*1e-15 / meter;
const double e = 0.0854;
const double m = 511e3;
const double alfa = 1.0 / 137.03599976;

//Laser params
const double eta = 1.0;
const double w0 = 1.55;
const double Estr = m * w0 * eta / e;
const double sigma = 3.0*2.0*pi / w0;
const double q = -e;
const double N = 5.0;


const double phase = -pi / 2.0;
const double poldegree = 1.0;


const double E0 = 30.0e9;
const double gam0 = E0 / m;
const double v0 = 1.0 - 0.5 / pow(gam0, 2.0);
const double tprim = 0.61e-3 / meter;
const double zpulseinit = tprim;
const int cores = 1;

const double L0 = pi / w0;


//const double wmax = 0.99*E0;
const double wmin = 0.01*E0;
const double wmax = 0.95*E0; // should be 0.8
const int rrstat = 0;
const int radiationstat = 1;
const int recoil = 1;

//const double dw = (wmax - wmin) / ((double)wfine*(double)cores);


//const double thetamax = 1.5*(1.0*eta + 1.0) * m/E0;
//const double thetamin = 1e-2*thetamax;
//const double phimin = 0.0;
//const double phimax = 1.0*pi;

const int pointsprim = 4000;
const int pointssec = 4000;
const int Ntot = 200000;
const int sections = 2;
const int wfine = nfreq;
const int freqiters = nfreq;


const double a0 = 5.29177e-11 / meter;
// SI 

//const double Z = 14.0;
//const double aScreen = 0.8853*a0*pow(1.0 + pow(Z, 2.0 / 3.0), -1.0 / 2.0);
//const double u = 0.106e-10 / meter;
//const double X0 = 9.37e-2 / meter;
//const double numberdensity = 5.011e28*pow(meter, 3.0);
//const double a[4] = { 2.1293*1e-10 / meter,2.5333*1e-10 / meter,0.8349*1e-10 / meter,0.3216*1e-10 / meter };
//const double B[4] = { 57.7748*pow(1e-10 / meter / (2.0*pi),2.0),16.4756*pow(1e-10 / meter / (2.0*pi),2.0),2.8796*pow(1e-10 / meter / (2.0*pi),2.0),0.3860*pow(1e-10 / meter / (2.0*pi),2.0) };
//const float aL = 5.431e-10 / meter; // changed

//const double r0 = 1.08e-10 / meter;
//const double dp = aL / pow(2.0, 1.5);
const double ISi = 172.0;
//const double U0 = 22.7;

//Ge
const double Z=32.0;
const double aScreen=0.8853*a0*pow(1.0+pow(Z,2.0/3.0),-1.0/2.0);
const double u=0.083e-10/meter;
const double X0=2.3e-2/meter;
const double numberdensity=4.39e28*pow(meter,3.0);
const double a[4]={2.4467*1e-10/meter,2.7015*1e-10/meter,1.6157*1e-10/meter,0.6009*1e-10/meter};
const double B[4]={55.893*pow(1e-10/meter/(2.0*pi),2.0),14.393*pow(1e-10/meter/(2.0*pi),2.0),2.4461*pow(1e-10/meter/(2.0*pi),2.0),0.3415*pow(1e-10/meter/(2.0*pi),2.0)};
const float aL=5.658e-10/meter; // changed
const double U0 = 35.73;
const double dp = aL / pow(2.0, 1.5);
const double convfac = 1.0;

const double R = 1.0*E0*dp / (2.0*U0);
const double c = m * pow(pi / alfa, 0.5) / pow(X0, 0.5);
const double cscat = c;
const int spinresolve = 0;
const int writestat = 1;


// end of constants