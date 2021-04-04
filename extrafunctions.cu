#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <thread>
#include <iomanip>
#include "functions.h"
//#include <armadillo>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuComplex.h>
#include <thrust/complex.h>
#include "constants.h"

using namespace std;
//using namespace arma;

double envelope(double x, double y, double z, double t) {
	return exp(-pow((t + z - zpulseinit), 2.0) / (2.0*sigma*sigma))+1.0*exp(-pow((t + z - zpulseinit-2.0*sigma), 2.0) / (0.5*sigma*sigma));
	//return 1.0;
}

double Fx(double x, double y, double z, double t) {

	return Estr * pow(cos(w0*(t + z)),1.0);
	//return Estr * cos(w0*(t + z - zpulseinit));
}

double Fy(double x, double y, double z, double t) {

	return poldegree*Estr * pow(cos(w0*(t + z)+phase), 1.0);
	//return Estr * cos(w0*(t + z - zpulseinit));
}

double ExDTplane(double x, double y, double z, double* rpass) {
	if (fabs(x) < *rpass) *rpass = fabs(x);

	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 2.0*pow(pi, 0.5)*e*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 0.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)))*2.0*x / (B[i] + pow(u, 2.0));
	}
	return out;
};

double dExdxDTplane(double x, double y, double z) {

	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 4.0*pow(pi, 0.5)*e*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 1.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)))*(1.0 - 2.0*pow(x, 2.0) / (B[i] + pow(u, 2.0)));
	}
	return out;
};

double UDTplane(double x, double y, double z, double* rpass) {
	if (abs(x) < *rpass) *rpass = abs(x);
	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 2.0*pow(pi, 0.5)*pow(e, 2.0)*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 0.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)));
	}
	return out;
};

double UDTplanesum(double x, double y, double z, double* rpass) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += UDTplane(x - i * dp, 0.0, 0.0, rpass);
	}
	return out;
}

double UDT(double x) {
	double *rpass = new double;
	*rpass = 10.0;
	double out = UDTplanesum(x, 0.0, 0.0, rpass) - UDTplanesum(dp / 2.0, 0.0, 0.0, rpass);
	delete rpass;
	return out;
}

double ExDTplanesum(double x, double y, double z, double* rpass) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += ExDTplane(x - i * dp, 0.0, 0.0, rpass);
	}
	return out;
}

double dExdxDTplanesum(double x, double y, double z) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += dExdxDTplane(x - i * dp, 0.0, 0.0);
	}
	return out;
}


double Ex(double x, double y, double z, double t) {

	
	//return Estr * cos(w0*(t + z));
	
	//double phi = w0 * (t + z)-2.0*2.0*pi;
	

	double phi = w0 * (t + z);
	double f = pow(sin(phi / (2.0*N)), 4.0);
	double fp = 2.0*pow(sin(phi / (2.0*N)), 3.0)*cos(phi / (2.0*N)) / N;
	return Estr * (f*sin(phi)-fp*cos(phi));
	//return Estr * sin(phi);
}

double Ey(double x, double y, double z, double t) {
	
	//double phi = w0 * (t + z) - 2.0*2.0*pi+phase;
	//return poldegree*Estr * sinh(phi) / (cosh(phi)*cosh(phi));

	//return Estr * sin(w0*(t + z));

	double phi = w0 * (t + z);
	double f = pow(sin(phi / (2.0*N)), 4.0);
	double fp = 2.0*pow(sin(phi / (2.0*N)), 3.0)*cos(phi / (2.0*N)) / N;
	return poldegree* Estr * (f*sin(phi+phase) - fp * cos(phi + phase));
	//return Estr * sin(phi + phase);
}

double Ez(double x, double y, double z, double t) {
	//return -(phifinal(x,y,z+potfine/2.0)-phifinal(x,y,z-potfine/2.0))/potfine;
	return 0.0;
}

double dExdx(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dExdxDTstringsum(x,y,z);
	return dExdxDTplanesum(x - z * z / (2.0*R), y, z);;
}

double dExdy(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dExdyDTstringsum(x,y,z);
	return 0.0;
}

double dEydy(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dEydyDTstringsum(x,y,z);
	return 0.0;
}

double Bx(double x, double y, double z, double t) {
	return 1.0*Ey(x, y, z, t);
}

double By(double x, double y, double z, double t) {
	return -1.0*Ex(x, y, z, t);
}

double Bz(double x, double y, double z) {
	return 0.0;
}
double Ex2(double x, double y, double z, double t, double Eparam) {
	return Eparam *( cos(w0*(t + z)));
}

double Ey2(double x, double y, double z, double t, double Eparam) {
	//return -(phifinal(x,y+potfine/2.0,z)-phifinal(x,y-potfine/2.0,z))/potfine;
	//return EyDTstringsum(x,y,z);
	return Eparam * pow(cos(w0*(t + z)+phase), 1.0);
}

double Ez2(double x, double y, double z, double t, double Eparam) {
	//return -(phifinal(x,y,z+potfine/2.0)-phifinal(x,y,z-potfine/2.0))/potfine;
	return 0.0;
}


double Bx2(double x, double y, double z, double t, double Eparam) {
	return Ey2(x, y, z, t, Eparam);
}

double By2(double x, double y, double z, double t, double Eparam) {
	return -Ex2(x, y, z, t, Eparam);
}

double Bz2(double x, double y, double z, double t, double Eparam) {
	return 0.0;
}

int func(double t, const double y[], double f[],
	void *params)
{
	double ptemp;
	double gammatemp;
	double *par = (double *)params;
	double gam = par[0];
	double *rpass = par + 1;

	//ptemp=pow(pow(y[0],2)+pow(y[1],2)+pow(y[2]+E*v0,2),0.5);
	gammatemp = y[6] + gam;
	double v0 = 1 - 0.5*pow(gam, -2.0);

	double crystalstat = 1.0;

	//double Fx=q*(Ex(t,y[3],y[4],y[5]+v0*t,y[7])-(y[2]+v0)*By(t,y[3],y[4],y[5]+v0*t,y[7]));
	//double Fy=q*(Ey(t,y[3],y[4],y[5]+v0*t,y[7])+(y[2]+v0)*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
	//double Fz=q*(y[0]*By(t,y[3],y[4],y[5]+v0*t,y[7])-y[1]*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));

	double Fx = q * (Ex(y[3], y[4], y[5] + v0 * t, t) - (y[2] + v0)*By(y[3], y[4], y[5] + v0 * t,t));
	double Fy = q * (Ey(y[3], y[4], y[5] + v0 * t,t) + (y[2] + v0)*Bx(y[3], y[4], y[5] + v0 * t,t));
	double Fz = q * Ez(y[3], y[4], y[5] + v0 * t,t) + q * (y[0] * By(y[3], y[4], y[5] + v0 * t,t) - y[1] * Bx(y[3], y[4], y[5] + v0 * t,t));

	
	//ADDING RADIATION REACTION
	/*if (rrstat == 1) {

		double chicur = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5)*gammatemp*e / (m*m);
		double baierfac = 1.0 / pow(1.0 + 4.8*(1.0 + chicur)*log(1.0 + 1.7*chicur) + 2.44*pow(chicur, 2.0), 2.0 / 3.0);
		baierfac = 1.0;
		Fx += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[0] * dExdx(y[3], y[4], y[5] + v0 * t) + y[1] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ex(y[3], y[4], y[5] + v0 * t, rpass)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[0] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));

		Fy += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[1] * dEydy(y[3], y[4], y[5] + v0 * t) + y[0] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ey(y[3], y[4], y[5] + v0 * t)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[1] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));

		Fz += baierfac * (-2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*(y[2] + v0)*(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
	}*/
	// DONE RR

	double gamprime = (y[0] * Fx + y[1] * Fy + (y[2] + v0)*Fz) / m;
	//gamprime=0.0;

	f[0] = Fx / (gammatemp*m) - gamprime / gammatemp * y[0];
	f[1] = Fy / (gammatemp*m) - gamprime / gammatemp * y[1];	
	f[2] = 1.0 / (gammatemp*m)*(Fz*(y[0] * y[0] + y[1] * y[1] + 1.0 / (gammatemp*gammatemp)) - Fx * y[0] - Fy * y[1]);

	f[3] = y[0];
	f[4] = y[1];	
	f[5] = 0.5*(pow(gam, -2.0) - pow(gammatemp, -2.0) - pow(y[0], 2.0) - pow(y[1], 2.0));

	


	f[6] = gamprime;


	return GSL_SUCCESS;
}

//int func2(double t, const double y[], double f[],
//	void *params)
//{
//	double ptemp;
//	double gammatemp;
//	double *par = (double *)params;
//	double gam = par[0];
//	double *rpass = par + 1;
//	double Eparamx = par[1];
//	double Eparamy = par[2];
//
//	//ptemp=pow(pow(y[0],2)+pow(y[1],2)+pow(y[2]+E*v0,2),0.5);
//	gammatemp = y[6] + gam;
//	double v0 = 1 - 0.5*pow(gam, -2.0);
//
//
//
//	//double Fx=q*(Ex(t,y[3],y[4],y[5]+v0*t,y[7])-(y[2]+v0)*By(t,y[3],y[4],y[5]+v0*t,y[7]));
//	//double Fy=q*(Ey(t,y[3],y[4],y[5]+v0*t,y[7])+(y[2]+v0)*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
//	//double Fz=q*(y[0]*By(t,y[3],y[4],y[5]+v0*t,y[7])-y[1]*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
//
//	double Fx = q * (Ex2(y[3], y[4], y[5] + v0 * t, t, Eparamx) - (y[2] + v0)*By2(y[3], y[4], y[5] + v0 * t, t, Eparamx));
//	double Fy = q * (Ey2(y[3], y[4], y[5] + v0 * t, t, Eparamy) + (y[2] + v0)*Bx2(y[3], y[4], y[5] + v0 * t, t, Eparamy));
//	double Fz = q * Ez2(y[3], y[4], y[5] + v0 * t, t, Eparamx) + q * (y[0] * By2(y[3], y[4], y[5] + v0 * t, t, Eparamx) - y[1] * Bx2(y[3], y[4], y[5] + v0 * t, t, Eparamy));
//
//
//	//ADDING RADIATION REACTION
//	//if (rrstat == 1) {
//	//	double chicur = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5)*gammatemp*e / (m*m);
//	//	double baierfac = 1.0 / pow(1.0 + 4.8*(1.0 + chicur)*log(1.0 + 1.7*chicur) + 2.44*pow(chicur, 2.0), 2.0 / 3.0);
//	//	//baierfac=1.0;
//	//	Fx += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[0] * dExdx(y[3], y[4], y[5] + v0 * t) + y[1] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ex(y[3], y[4], y[5] + v0 * t, rpass)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[0] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//
//	//	Fy += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[1] * dEydy(y[3], y[4], y[5] + v0 * t) + y[0] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ey(y[3], y[4], y[5] + v0 * t)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[1] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//
//	//	Fz += baierfac * (-2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*(y[2] + v0)*(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//	//}
//	// DONE RR
//
//	double gamprime = (y[0] * Fx + y[1] * Fy + (y[2] + v0)*Fz) / m;
//	//gamprime=0.0;
//
//	f[0] = Fx / (gammatemp*m) - gamprime / gammatemp * y[0];
//	f[1] = Fy / (gammatemp*m) - gamprime / gammatemp * y[1];
//	f[2] = Fz / (gammatemp*m) - gamprime / gammatemp * (y[2] + v0);
//
//	f[3] = y[0];
//	f[4] = y[1];
//	f[5] = y[2];
//	f[6] = gamprime;
//
//
//	return GSL_SUCCESS;
//}

void W(double* y, double gam, double* omega, double* dEdtspinupup1, double* dEdtspindowndown1, double* dEdtspinupdown1, double* dEdtspindownup1, double* dEdtspinupup2, double* dEdtspindowndown2, double* dEdtspinupdown2, double* dEdtspindownup2, double* radresult, double* thetas, double l_f, int curN, int sectionN, double t0) {

	ofstream orbit;
	char a[150];
	char b[150];
	strcpy_s(a, "orbitsec");
	sprintf_s(b, "%i", curN);
	strcat_s(a, b);
	sprintf_s(b, "-%i.txt", sectionN);
	strcat_s(a, b);
	
	if (writestat == 1) orbit.open(a);
	if (writestat == 1) orbit << scientific << setprecision(16);

	const gsl_rng_type * Trand;
	gsl_rng * r;
	gsl_rng_env_setup();

	Trand = gsl_rng_default;
	r = gsl_rng_alloc(Trand);
	gsl_rng_set(r, curN);

	double *params = new double[2];
	params[0] = gam;
	double *rpass = params + 1;
	double t = t0;
	int points = pointssec;
	double *trajec = new double[points * 8];
	//double *out = new double[freqiters]();
	double *finalout = new double[freqiters]();

	double gamtemp = y[6] + gam;
	trajec[0 * points] = t;
	trajec[1 * points] = y[0];
	trajec[2 * points] = y[1];
	trajec[3 * points] = y[2];
	trajec[4 * points] = y[3];
	trajec[5 * points] = y[4];
	trajec[6 * points] = y[5];
	trajec[7 * points] = y[6];

	double thetax_max = y[0];
	double thetax_min = y[0];
	double thetay_max = y[1];
	double thetay_min = y[1];

	gsl_odeiv2_system sys = { func, NULL, 7, params };
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-8, 1e-8, 1e-8);

	int i, status;
	double ti;
	for (i = 1; i <= points - 1; i++)
	{
		*rpass = aL;
		ti = t0+(double)i * l_f / ((double)points - 1);
		while (t < ti) {
			//status = gsl_odeiv_evolve_apply (e, c, s,&sys,&t,ti,&h, y);
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);
			if (status != GSL_SUCCESS)
			{
				printf("error, return value=%d\n", status);
				break;
			}
		}
		//   printf ("%.5e %.5e %.5e\n", t, y[3], y[4], y[5]);


		double scatvar = fmod(abs(y[3]), dp);
		double Gauss = dp / (u*pow(2.0*pi, 0.5))*(exp(-pow(scatvar, 2.0) / (2.0*u*u)) + exp(-pow(dp - scatvar, 2.0) / (2.0*u*u)));
		double nel = -1.0*dExdx(scatvar, 0.0, 0.0) / (4.0*pi*pow(alfa, 0.5)) + Z * numberdensity*Gauss;
		double sigsqncl = Gauss * pow(cscat / (m*(gam + y[6])), 2.0)*l_f / (double)pointssec;
		double sigsqel = 1.0*pi*pow(alfa, 2.0) / pow(m*(gam + y[6]), 2.0)*(log(2.0*m*pow(gam + y[6], 2.0) / ISi) - 1.0) * nel *l_f / (double)pointssec;
		//sigsqncl = pow(cscat / (m*(gam + y[6])), 2.0)*tprim / (double)pointssec;
		//sigsqel = 1.0*pi*pow(alfa, 2.0) / pow(m*(gam + y[6]), 2.0)*(log(2.0*m*pow(gam + y[6], 2.0) / ISi) - 1.0) * Z*numberdensity*tprim / (double)pointssec;
		//cout << sigsqncl / sigsqel << endl;
		double scatx = pow(1.0*sigsqel + 1.0*sigsqncl, 0.5);
		double scaty = scatx;

		scatx = 0.0*scatx*gsl_ran_gaussian(r, 1.0);
		scaty = 0.0*scaty*gsl_ran_gaussian(r, 1.0);
		//y[2] = y[2] - 0.5*(pow(scatx, 2.0) + pow(scaty, 2.0)) - y[0] * scatx - y[1] * scaty;
		y[0] = y[0] + scatx;
		y[1] = y[1] + scaty;
		y[2] = 0.5*(pow(gam, -2.0) - pow(gamtemp, -2.0) - pow(y[0], 2.0) - pow(y[1], 2.0));

		//if (0) {

		//	double gamold = y[6] + gam;
		//	//cout << scatx * gamold << endl;
		//	double deltagam = -2.0*alfa / (3.0*pi)*pow(gamold, 3.0)*(pow(scatx, 2.0) + pow(scaty, 2.0))*3.0 / 4.0;
		//	double vrat = 1.0 + deltagam / pow(gamold, 3.0);
		//	y[0] = y[0] * vrat;
		//	y[1] = y[1] * vrat;
		//	y[2] += deltagam / pow(gamold, 3.0);
		//	y[6] += deltagam;
		//}

		//double field = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5);
		//double chi = gamtemp * field*e / pow(m, 2.0);

		double Emek = UDT(y[3]) + 0.5*gamtemp*m*y[0] * y[0];
		
			if (writestat == 1 && i % 1 == 0) orbit << y[3] << " " << y[4] << " " << y[5] << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[6] + gam << " " << t << " " << " " << Emek << endl;
		


		trajec[i + 0 * points] = t;
		trajec[i + 1 * points] = y[0];
		trajec[i + 2 * points] = y[1];
		trajec[i + 3 * points] = y[2];
		trajec[i + 4 * points] = y[3];
		trajec[i + 5 * points] = y[4];
		trajec[i + 6 * points] = y[5];
		trajec[i + 7 * points] = y[6];



		//chi=0.01;
		//gamtemp=gam0;



		if (thetax_max < y[0]) thetax_max = y[0];
		if (thetay_max < y[1]) thetay_max = y[1];

		if (thetax_min > y[0]) thetax_min = y[0];
		if (thetay_min > y[1]) thetay_min = y[1];

	}//end time step

	if (writestat == 1) orbit.close();

	thetax_max += 3.0 / gamtemp;
	thetax_min -= 3.0 / gamtemp;
	thetay_max += 3.0 / gamtemp;
	thetay_min -= 3.0 / gamtemp;

	/*thetax_max = (thetax_max + thetax_min) / 2.0 + 0.1/gamtemp;
	thetax_min = (thetax_max + thetax_min) / 2.0 - 0.1 / gamtemp;
	thetay_max = (thetay_max + thetay_min) / 2.0 + 0.1 / gamtemp;
	thetay_min = (thetay_max + thetay_min) / 2.0 - 0.1 / gamtemp;*/



	thetax_max = (1.0*eta + 3.0) / gamtemp; // 1.0*eta+3.0
	thetax_min = -(1.0*eta + 3.0) / gamtemp;
	thetay_max = thetax_max;
	thetay_min = thetax_min;

	//thetay_min = 0.0;

	//if (thetax_min < tprim / (2.0*R)) thetax_min = tprim / (2.0*R);
	//if (thetax_max < thetax_min) thetax_max = thetax_min + 2.0 / gamtemp;

	thetas[0] = thetax_min;
	thetas[1] = thetax_max;
	thetas[2] = thetay_min;
	thetas[3] = thetay_max;

	gsl_odeiv2_driver_free(d);



	
	double dw = (wmax - wmin) / ((double)wfine);
	
	double *consts = new double[9];
	consts[0] = m;
	consts[1] = wmin;
	consts[2] = wmax;
	consts[3] = thetax_min;
	consts[4] = thetax_max;
	consts[5] = thetay_min;
	consts[6] = thetay_max;
	consts[7] = gamtemp;
	consts[8] = (double)pointssec;

	//calculator2(consts, trajec, dEdt, thetaxmin, thetaxmax, thetaymin, thetaymax);
	calculator3(consts, trajec, dEdtspinupup1, dEdtspindowndown1, dEdtspinupdown1, dEdtspindownup1, dEdtspinupup2, dEdtspindowndown2, dEdtspinupdown2, dEdtspindownup2);

	//double dw = (wmax - wmin) / ((double)wfine);

	for (int l = 0; l < wfine; l++) {
		omega[l] = wmin + dw * (double)l; //Sets the current frequency to be calculated
	}

	gsl_rng_free(r);

	delete[] trajec;
	delete[] finalout;
	delete[] params;

}



int myintmod(int a, int b) {
	int x;
	if (a > 0) x = a % b;
	else if (a < 0) x = a % b + b;
	else x = 0;
	return x;
};



double lnfac(int n) {
double out;
double x=(double)n;

if (n!=0) {
out=x*log(x)-x+0.5*log(2.0*pi*x)+1.0/(12.0*x)-1.0/( 360.0*pow(x,3.0) );
}
else out=0;

return out;
};

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
};

double mylag(int n, int alfa, double x) {
double out=0.0;
for (int i=0; i<=n; i++) {
out+=pow(-1.0,i)*factorial(n+alfa)/factorial(alfa+i)/factorial(n-i)*pow(x,(double)i)/factorial(i);
};

return out;
}

double trapz (double *t, double *f, int points) {
double dummy=0;
for (int i=0; i<points-1; i++) {
if ( (t[i+1]-t[i])>0.0 ) dummy+=(t[i+1]-t[i])*(f[i+1]+f[i])/2.0;
}
return dummy;
};

void mycross (double *a, double *b, double *c) {
c[1]=a[2]*b[3]-a[3]*b[2];
c[2]=a[3]*b[1]-a[1]*b[3];
c[3]=a[1]*b[2]-a[2]*b[1];
};

void normalize (double *a) {
double x=pow(a[1]*a[1]+a[2]*a[2]+a[3]*a[3],0.5);
a[1]=a[1]/x;
a[2]=a[2]/x;
a[3]=a[3]/x;
};

double mymod(double a, double b) {
	double x;
	if (a >= 0) x = fmod(a, b);
	else x = fmod(a, b) + b;
	return x;
};



int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

double fdot(double *x1, double *x2) {
double out=x1[0]*x2[0];
for (int i=1; i<4; i++) {
out-=x1[i]*x2[i];
}
return out;
}






double myinterp(double *a, double *b, double x, int points) {


	bool success = 0;

	if (x > a[points - 1]) {
		return b[points - 1];
	}

	if (x < a[0]) {
		return b[0];
	}

	else {
		for (int i = 0; i < points - 1; i++) {
			if (x >= a[i] && x < a[i + 1]) {
				success = 1;
				return b[i] + (b[i + 1] - b[i]) / (a[i + 1] - a[i])*(x - a[i]);
			}
		}
	}

	if (success == 0) {
		return 0;
	}
}