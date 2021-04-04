#define PI 3.14159265358979323846264338327950288419716939937510582
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf.h>
//#include <cmath>
#include <iostream>
#include <fstream>
//#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <random>
#include "functions.h"

//#include <cuComplex.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "constants.h"
#include <cusolverDn.h>
#include <time.h>
#include "math.h"
#include <thrust/device_vector.h>



using namespace std;





__host__ __device__ double mysinc(double x) {
	double out = 0.0;
	if (abs(x)> 1e-3) {
		out = sin(x) / x;
	}

	else out = 1.0;
	return out;
};

__host__ __device__ double myP(double x) {
	double out = 0.0;
	if (abs(x) > 1e-3) {
		out = (1.0-cos(x)) / x;
	}

	else out = x/2.0;
	return out;
};

__device__ thrust::complex<double> Gfunc(double a, double b, double L) {
	thrust::complex<double> out(0.0,0.0);
	thrust::complex<double> iunit(0.0, 1.0);
	if (a == 0.0 && b!=0.0) {
		out = (exp(iunit*b*L)*(1.0 - iunit * b*L) - 1.0) / (b*b);
	}
	else if (b == 0.0 && a != 0.0) {
		out = 1.0 / (iunit*a)*( (exp(iunit*a*L)-1.0)/(iunit*a)-L);
	}
	else if (a == 0.0 && b == 0.0) {
		out = L * L / 2.0;
	}
	else {
		out = 1.0 / (iunit*a)*( (exp(iunit*(a+b)*L)-1.0)/(iunit*(a+b))-(exp(iunit*b*L)-1.0)/(iunit*b) );
	}
	return out;
};

//__device__ double mysinc_dev(double x) {
//	double out = 0.0;
//	if (x != 0.0) {
//		out = sin(x) / x;
//	}
//
//	else out = 1.0;
//	return out;
//};

double lnfac(double x) {
double out=0.0;
if (x!=0.0) {
out=x*log(x)-x+0.5*log(2.0*pi*x)+1.0/(12.0*x)-1.0/(360.0*pow(x,3.0));
}
return out;
};

__device__ void mycross_dev(double *a, double *b, double *c) {
	c[1] = a[2] * b[3] - a[3] * b[2];
	c[2] = a[3] * b[1] - a[1] * b[3];
	c[3] = a[1] * b[2] - a[2] * b[1];
};

__device__ thrust::complex<double> M_dev(double *const_dev ,double* pveci, double* kvec, int ni, int nf, double kbi, double kbf, int n, double* MATi, double* MATf, double* VECi, double* VECf, int si, int sf, int sfot) {

	//arma::mat eigveci(MATi, 2 * N + 1, 2 * N + 1, false, true);
	double m = const_dev[0];
	double k0 = const_dev[1];
	double wmin = const_dev[2];
	double wmax = const_dev[3];
	double phimin = const_dev[4];
	double phimax = const_dev[5];
	double pxi = const_dev[6];
	double kbmin = const_dev[7];
	double kbmax = const_dev[8];
	double thetamax = const_dev[9];
	double e = const_dev[10];
	double pi = const_dev[11];
	
	int N = lround(const_dev[14]);

	double pxf = pveci[1] - kvec[1];
	double pzf = pveci[3] - kvec[3];
	
	double Ei = pow(2.0*pveci[1] *VECi[ni] + pveci[1] * pveci[1] + pveci[3]* pveci[3] + m * m, 0.5); // there seems to be a mistake here? pxi should be pveci[1]
	double Ef = pow(2.0*pxf*VECf[nf] + pxf * pxf + pzf* pzf + m * m, 0.5);



	/*double *pol1 = new double[4]();
	double *pol2 = new double[4]();*/

	double pol1[4];
	double pol2[4];


	/*mycross(kvec, pveci, pol1);
	normalize(pol1);

	mycross(kvec, pol1, pol2);
	normalize(pol2);*/

	double thetasq = (kvec[2]* kvec[2] + kvec[3]* kvec[3]) / (kvec[0]* kvec[0]);
	double theta = sqrt(thetasq);
	double cphi = kvec[2] / (kvec[0] * theta);
	double sphi = kvec[3] / (kvec[0] * theta);

	pol1[2] = sphi;
	pol1[3] = -cphi;

	pol2[1] = -theta;
	pol2[2] = cphi;
	pol2[3] = sphi;

	double *pol;

	/*if (sfot == 1) {
		pol = pol1;
	}
	else {
		pol = pol2;
	}*/

	double *pols[2];
	pols[0] = pol1;
	pols[1] = pol2;
	pol = pols[sfot];

	/*double *Pi = new double[4]();
	double *Pf = new double[4]();*/

	double Pi[4];
	double Pf[4];

	
	Pi[0] = 0.0;
	Pi[1] = 0.0;
	Pi[2] = 0.0;
	Pi[3] = 0.0;
	
	Pf[0] = 0.0;
	Pf[1] = 0.0;
	Pf[2] = 0.0;
	Pf[3] = 0.0;

	double prodsum = 0.0;

	int mimin, mimax;
	if (n >= 0) {
		mimin = -N;
		mimax = N - n;
	}
	else {
		mimin = -N - n;
		mimax = N;
	}
	for (int mi = mimin; mi <= mimax; mi++) {

		int mf = mi + n;
		//double cmi = eigveci(mi + N, ni);
		//double cmf = eigvecf(mf + N, nf);
		double cmi = MATi[mi + N + (2 * N + 1)*ni];
		double cmf = MATf[mf + N + (2 * N + 1)*nf];

		double py_i = k0 * (double)mi + kbi;
		double py_f = k0 * (double)mf + kbf;

		double px_itemp = pveci[1] + VECi[ni] - py_i * py_i / (2.0*Ei);
		double px_ftemp = pxf + VECf[nf] - py_f * py_f / (2.0*Ef);

		Pi[1] += cmi * cmf*px_itemp;
		Pi[2] += cmi * cmf*py_i;
		prodsum += cmi * cmf;

		Pf[1] += cmi * cmf* px_ftemp;
		Pf[2] += cmi * cmf* py_f;


	} // end ns sum

	Pi[1] = Pi[1] / (Ei + m);
	Pi[2] = Pi[2] / (Ei + m);
	Pi[3] = prodsum * pveci[3] / (Ei + m);

	Pf[1] = Pf[1] / (Ef + m);
	Pf[2] = Pf[2] / (Ef + m);
	Pf[3] = prodsum * pzf / (Ef + m);

	/*double *A = new double[4]();
	double *B = new double[4]();
	double *B2 = new double[4]();*/

	double A[4];
	double B[4];
	double B2[4];

	A[1] = Pi[1] + Pf[1];
	A[2] = Pi[2] + Pf[2];
	A[3] = Pi[3] + Pf[3];

	B2[1] = Pi[1] - Pf[1];
	B2[2] = Pi[2] - Pf[2];
	B2[3] = Pi[3] - Pf[3];

	mycross_dev(pol, B2, B);

	double realcomp = pol[1] * A[1] + pol[2] * A[2] + pol[3] * A[3];

	//delete[] pol1;
	//delete[] pol2;

	thrust::complex<double> out(0.0,0.0);

	if (si == 1 && sf == 1) out = thrust::complex<double>(realcomp, B[3]); // up-up case
	else if (si == 1 && sf == 0) out = thrust::complex<double>(-B[2], B[1]); // up-down
	else if (si == 0 && sf == 0) out = thrust::complex<double>(realcomp, -B[3]); // down-down
	else if (si == 0 && sf == 1) out = thrust::complex<double>(B[2], B[1]); // down-up
	
	/*if (si == 1 && sf == 1) {
		result[0] = realcomp;
		result[1] = B[3];
	}
	else if (si == 1 && sf == 0) {
		result[0]= -B[2];
		result[1] = B[1];
	}
	else if (si == 0 && sf == 0) {
		result[0] = realcomp;
		result[1] = -B[3];
	}
	else if (si == 0 && sf == 1) {
		result[0] = B[2];
		result[1] = B[1];
	}
*/
																		

	/*delete[] B2;
	delete[] B;
	delete[] A;*/
	/*delete[] Pi;
	delete[] Pf;*/

	//out = thrust::complex<double>(realcomp, B[1] + B[2] + B[3]);

	return out;
}



__device__ double mymod_dev(double a, double b) {
	double x;
	if (a >= 0) x = fmod(a, b);
	else x = fmod(a, b) + b;
	return x;
};


__global__ void myradKernel(double *const_dev, double *trajec_dev, double *result_dev) {
	
	double m = const_dev[0];
	double wmin = const_dev[1];
	double wmax = const_dev[2];
	double thetaxmin = const_dev[3];
	double thetaxmax = const_dev[4];
	double thetaymin = const_dev[5];
	double thetaymax = const_dev[6];
	double gam0 = const_dev[7];
	int points = round(const_dev[8]);

	
	thrust::complex<double> resultx(0.0, 0.0);
	thrust::complex<double> resulty(0.0, 0.0);
	thrust::complex<double> resultJ(0.0, 0.0);
	

	int id = threadIdx.x + blockIdx.x*blockDim.x;

	int a, b, c;
	a = id % thetaxfine;
	b = (id % (thetaxfine*thetayfine) - a)/ thetaxfine;
	c = (id - a - b * thetaxfine) / (thetaxfine*thetayfine);

	double w = (double)c/(double)nfreq*(wmax - wmin) + wmin;
	double thetay = (double)b/(double)thetayfine*(thetaymax - thetaymin) + thetaymin;
	double thetax = (double)a/(double)thetaxfine*(thetaxmax - thetaxmin) + thetaxmin;

	double nx = thetax;
	double ny = thetay;

	for (int i = 0; i < points - 1; i++) {
		double vx = trajec_dev[1 * points + i];
		double vy = trajec_dev[2 * points + i];
		double vz = trajec_dev[3 * points + i];
		double dt = trajec_dev[i + 1] - trajec_dev[i];

		double betaxdot = (trajec_dev[1 * points + i + 1] - trajec_dev[1 * points + i]) / dt;
		double betaydot = (trajec_dev[2 * points + i + 1] - trajec_dev[2 * points + i]) / dt;
		double betazdot = (trajec_dev[3 * points + i + 1] - trajec_dev[3 * points + i]) / dt;



		double t = trajec_dev[i];

		double x = trajec_dev[4 * points + i];
		double y = trajec_dev[5 * points + i];
		double z = trajec_dev[6 * points + i];
		double deltagam = trajec_dev[7 * points + i];

		double betaortho = sqrt(vx*vx + vy * vy);

		// Classical
		//(*exponent) = { 0.0,*omega*(-*nx**x - *ny**y - *z + 0.5*(pow(*northo,2) + pow(gam,-2))**t) };


		double Ecur = (gam0 + deltagam)*m;
		double Ep = Ecur - w;
		// Strong field
		thrust::complex<double> exponent(0.0, w*Ecur / (Ecur - w)*(-nx * x - ny * y - z + 0.5*(nx*nx+ny*ny + 1.0 / (gam0*gam0))*t));
		thrust::complex<double> faktor(1.0 / ((0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)*(0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)), 0.0);

		if (w < Ecur) {


			thrust::complex<double> v1(ny*(betaydot*(nx - vx) - betaxdot * (ny - vy)) - betaxdot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) + (nx - vx)*betazdot, 0.0);
			thrust::complex<double> v2((ny - vy)*betazdot - betaydot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) - nx * ((nx - vx)*betaydot - (ny - vy)*betaxdot), 0.0);

			double scale = sqrt((Ecur*Ecur + Ep * Ep) / (2.0*Ecur*Ecur)); // Strong field
			//scale = { 1.0,0.0 }; // Classical

			v1 = v1 * faktor*exp(exponent)*scale;
			v2 = v2 * faktor*exp(exponent)*scale;

			thrust::complex<double> vJ(nx*(betaxdot)+(ny)*(betaydot)+betazdot, 0.0);




			double scaleJ = w * m / (Ecur*Ecur)*0.70710678118;

			vJ = vJ * faktor*exp(exponent)*scaleJ;

			resultx += v1 * dt;
			resulty += v2 * dt;
			resultJ += vJ * dt;

		}

	}
		result_dev[id] = (thrust::norm(resultx)+ thrust::norm(resulty)+ thrust::norm(resultJ));

	//result_dev[id] = thrust::norm(resultx);
}




__global__ void myradKernelspin(double *const_dev, double *trajec_dev, double *result_dev) {

	double m = const_dev[0];
	double wmin = const_dev[1];
	double wmax = const_dev[2];
	double thetaxmin = const_dev[3];
	double thetaxmax = const_dev[4];
	double thetaymin = const_dev[5];
	double thetaymax = const_dev[6];
	double gam0 = const_dev[7];
	int points = round(const_dev[8]);
	double vorthosqavg = const_dev[9];




	thrust::complex<double> resultIx(0.0, 0.0);
	thrust::complex<double> resultIy(0.0, 0.0);
	thrust::complex<double> resultJ(0.0, 0.0);


	int id = threadIdx.x + blockIdx.x*blockDim.x;

	int a, b, c;
	a = id % thetaxfine;
	b = (id % (thetaxfine*thetayfine) - a) / thetaxfine;
	c = (id - a - b * thetaxfine) / (thetayfine*thetaxfine);	

	double w = (double)c / (double)nfreq*(wmax - wmin) + wmin;
	double thetay = (double)b / ((double)thetayfine-1.0)*(thetaymax - thetaymin) + thetaymin;
	double thetax = (double)a / ((double)thetaxfine-1.0)*(thetaxmax - thetaxmin) + thetaxmin;

	double w0 = 2.0*PI / (trajec_dev[points - 1] - trajec_dev[0]);
	double wp = w * gam0*m / (gam0*m - w);

	double Ecur = gam0 * m;
	double Ep = Ecur - w;

	double thetasq = thetax*thetax+thetay*thetay;
	double theta = sqrt(thetasq);
	
		

		double nx = thetax;
		double ny = thetay;

		for (int i = 0; i < points - 1; i++) {
			double vx = trajec_dev[1 * points + i];
			double vy = trajec_dev[2 * points + i];
			double vz = trajec_dev[3 * points + i];
			double dt = trajec_dev[i + 1] - trajec_dev[i];

			double betaxdot = (trajec_dev[1 * points + i + 1] - trajec_dev[1 * points + i]) / dt;
			double betaydot = (trajec_dev[2 * points + i + 1] - trajec_dev[2 * points + i]) / dt;
			double betazdot = (trajec_dev[3 * points + i + 1] - trajec_dev[3 * points + i]) / dt;



			double t = trajec_dev[i];

			double x = trajec_dev[4 * points + i];
			double y = trajec_dev[5 * points + i];
			double z = trajec_dev[6 * points + i];
			double deltagam = trajec_dev[7 * points + i];

			double betaortho = sqrt(vx*vx + vy * vy);

			// Classical
			//(*exponent) = { 0.0,*omega*(-*nx**x - *ny**y - *z + 0.5*(pow(*northo,2) + pow(gam,-2))**t) };


			double Ecur = (gam0 + deltagam)*m;
			double Ep = Ecur - w;
			// Strong field
			thrust::complex<double> exponent(0.0, w*Ecur / (Ecur - w)*(-nx * x - ny * y - z + 0.5*(nx*nx + ny * ny + 1.0 / (gam0*gam0))*t));
			thrust::complex<double> faktor(1.0 / ((0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)*(0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)), 0.0);

			if (w < Ecur) {


				thrust::complex<double> v1(ny*(betaydot*(nx - vx) - betaxdot * (ny - vy)) - betaxdot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) + (nx - vx)*betazdot, 0.0);
				thrust::complex<double> v2((ny - vy)*betazdot - betaydot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) - nx * ((nx - vx)*betaydot - (ny - vy)*betaxdot), 0.0);

				//double scale = sqrt((Ecur*Ecur + Ep * Ep) / (2.0*Ecur*Ecur)); // Strong field
				//scale = { 1.0,0.0 }; // Classical

				v1 = v1 * faktor*exp(exponent);
				v2 = v2 * faktor*exp(exponent);

				thrust::complex<double> vJ(nx*(betaxdot)+(ny)*(betaydot)+betazdot, 0.0);




				//double scaleJ = w * m / (Ecur*Ecur)*0.70710678118;

				vJ = vJ * faktor*exp(exponent);

				resultIx += v1 * dt;
				resultIy += v2 * dt;
				resultJ += vJ * dt;


			}
		}

		thrust::complex<double> si_u;
		thrust::complex<double> si_d;
		thrust::complex<double> sf_u;
		thrust::complex<double> sf_d;
		thrust::complex<double> polupx;
		thrust::complex<double> polupy;
		thrust::complex<double> polupz;
		thrust::complex<double> poldownx;
		thrust::complex<double> poldowny;
		thrust::complex<double> poldownz;
		thrust::complex<double> polx;
		thrust::complex<double> poly;
		thrust::complex<double> polz;

		
		//circular
		polupx = pow(2.0,-0.5)*thrust::complex<double>(1.0,0.0); // Left handed i.e. positive helicity
		polupy = pow(2.0, -0.5)*thrust::complex<double>(0.0, 1.0);
		polupz = pow(2.0, -0.5)*thrust::complex<double>(-nx, -ny);

		poldownx = pow(2.0, -0.5)*thrust::complex<double>(1.0, 0.0);
		poldowny = pow(2.0, -0.5)*thrust::complex<double>(0.0, -1.0);
		poldownz = pow(2.0, -0.5)*thrust::complex<double>(-nx, ny);


		//linear
		/*polupx = thrust::complex<double>(1.0, 0.0);
		polupy = thrust::complex<double>(0.0, 0.0);
		polupz = thrust::complex<double>(-nx, 0.0);

		poldownx = thrust::complex<double>(0.0, 0.0);
		poldowny = thrust::complex<double>(1.0, 0.0);
		poldownz = thrust::complex<double>(-ny, 0.0);*/
					
		double C1 = Ecur / (2.0*pow(Ecur*Ep, 0.5))*(pow((Ep+m)/(Ecur+m),0.5) + pow((Ecur + m) / (Ep + m), 0.5));
		double D1 = Ecur / (2.0*pow(Ecur*Ep, 0.5))*(pow((Ep + m) / (Ecur + m), 0.5) - pow((Ecur + m) / (Ep + m), 0.5));
		double D2 = w / (2.0*pow(Ecur*Ep, 0.5))*pow((Ecur + m) / (Ep + m), 0.5);
		result_dev[id] = 0.0;
		
		
		for (int d = 0; d < 8; d++) {
			if (d == 0) { //up-up
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(1.0 , 0.0);
				si_d = thrust::complex<double>(0.0, 0.0);
				sf_u = thrust::complex<double>(1.0, 0.0);
				sf_d = thrust::complex<double>(0.0, 0.0);

				polx = polupx;
				poly = polupy;
				polz = polupz;
				
			}

			if (d == 1) { //down-down
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(0.0, 0.0);
				si_d = thrust::complex<double>(1.0, 0.0);
				sf_u = thrust::complex<double>(0.0, 0.0);
				sf_d = thrust::complex<double>(1.0, 0.0);
				
				polx = polupx;
				poly = polupy;
				polz = polupz;
				
			}

			if (d == 2) { //up-down
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(1.0, 0.0);
				si_d = thrust::complex<double>(0.0, 0.0);
				sf_u = thrust::complex<double>(0.0, 0.0);
				sf_d = thrust::complex<double>(1.0, 0.0);
				polx = polupx;
				poly = polupy;
				polz = polupz;
				
			}

			if (d == 3) { //down-up
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(0.0, 0.0);
				si_d = thrust::complex<double>(1.0, 0.0);
				sf_u = thrust::complex<double>(1.0, 0.0);
				sf_d = thrust::complex<double>(0.0, 0.0);

				polx = polupx;
				poly = polupy;
				polz = polupz;
				
			}

			if (d == 4) { //up-up
			//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(1.0 , 0.0);
				si_d = thrust::complex<double>(0.0, 0.0);
				sf_u = thrust::complex<double>(1.0, 0.0);
				sf_d = thrust::complex<double>(0.0, 0.0);

				polx = poldownx;
				poly = poldowny;
				polz = poldownz;

			}

			if (d == 5) { //down-down
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(0.0, 0.0);
				si_d = thrust::complex<double>(1.0, 0.0);
				sf_u = thrust::complex<double>(0.0, 0.0);
				sf_d = thrust::complex<double>(1.0, 0.0);

				polx = poldownx;
				poly = poldowny;
				polz = poldownz;

			}

			if (d == 6) { //up-down
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(1.0, 0.0);
				si_d = thrust::complex<double>(0.0, 0.0);
				sf_u = thrust::complex<double>(0.0, 0.0);
				sf_d = thrust::complex<double>(1.0, 0.0);
				polx = poldownx;
				poly = poldowny;
				polz = poldownz;

			}

			if (d == 7) { //down-up
				//Spin along y
				/*si_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				si_d = thrust::complex<double>(0.0, -1.0 / sqrt(2.0));
				sf_u = thrust::complex<double>(1.0 / sqrt(2.0), 0.0);
				sf_d = thrust::complex<double>(0.0, 1.0 / sqrt(2.0));*/

				//Spin along z
				si_u = thrust::complex<double>(0.0, 0.0);
				si_d = thrust::complex<double>(1.0, 0.0);
				sf_u = thrust::complex<double>(1.0, 0.0);
				sf_d = thrust::complex<double>(0.0, 0.0);

				polx = poldownx;
				poly = poldowny;
				polz = poldownz;

			}

			thrust::complex<double> iunit(0.0, 1.0);
			thrust::complex<double> Bx = (D1*resultIx-(D1+D2)*nx*resultJ);
			thrust::complex<double> By = (D1*resultIy - (D1 + D2)*ny*resultJ);
			thrust::complex<double> Bz = 1.0*(0.0 - (D1 + D2)*resultJ);
			

				thrust::complex<double> Dx = conj(poly) * Bz - conj(polz) * By;
				thrust::complex<double> Dy = conj(polz) * Bx - conj(polx) * Bz;
				thrust::complex<double> Dz = 1.0*(conj(polx) * By - conj(poly) * Bx);

				thrust::complex<double> Ccomp = C1 * (conj(polx)*resultIx + conj(poly) * resultIy);
				thrust::complex<double> imagu(0.0, 1.0);

				result_dev[id+d* thetayfine*thetaxfine*nfreq] += thrust::norm((conj(sf_u)*si_u + conj(sf_d)*si_d)*Ccomp + conj(sf_u)*(imagu*Dz*si_u + (imagu*Dx + Dy)*si_d) + conj(sf_d)*((imagu*Dx - Dy)*si_u - imagu * Dz*si_d));

			
			result_dev[id + d * thetayfine*thetaxfine*nfreq] = result_dev[id+d* thetayfine*thetaxfine*nfreq] * w*w/(wp*wp);

		}
		

}


cudaError_t calculator2(double *consts, double *trajec, double *out, double thetaxmin, double thetaxmax, double thetaymin, double thetaymax) {
	cudaError_t cudaStatus;
	//double *consts = new double[9];

			int nDevices;
			cudaGetDeviceCount(&nDevices);
			for (int i = 0; i < nDevices; i++) {
				cudaDeviceProp prop;
				cudaGetDeviceProperties(&prop, i);
				printf("Device Number: %d\n", i);
				printf("  Device name: %s\n", prop.name);
				printf("  Memory Clock Rate (KHz): %d\n",
					prop.memoryClockRate);
				printf("  Memory Bus Width (bits): %d\n",
					prop.memoryBusWidth);
				printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
					2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);
			}
		

		cudaStatus = cudaSetDevice(0);
		//cudaDeviceReset();

	//int devID = 0;
	//int argc;
	//char argv[1000];
	//devID =findCudaDevice(argc, (const char **)argv);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}



	double *const_dev;
	double *result_dev;
	double *trajec_dev;

	cudaStatus = cudaMalloc((void**)&const_dev, 9 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&result_dev, thetaxfine*thetayfine*nfreq * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&trajec_dev, 8*pointsprim * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(const_dev, consts, 9 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 1 failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(trajec_dev, trajec, 8*pointsprim * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 2 failed!");
		goto Error;
	}

	

	myradKernel << <BLOCKS, THREADS >> > (const_dev, trajec_dev, result_dev);

	cudaStatus = cudaPeekAtLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "myKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}


	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching myKernel!\n", cudaStatus);
		goto Error;
	}


	double *result = new double[thetaxfine*thetayfine*nfreq]();

	cudaStatus = cudaMemcpy(result, result_dev, thetaxfine*thetayfine*nfreq * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 1 failed!");
		goto Error;
	}

	for (int i = 0; i < nfreq; i++) {
		for (int j = 0; j < thetayfine; j++) {
			for (int l = 0; l < thetaxfine; l++) {
				int id = i * thetaxfine*thetayfine + j * thetaxfine + l;
				out[i] += result[id]*e*e/(4.0*pi*pi)*(thetaxmax-thetaxmin)/(double)thetaxfine*(thetaymax-thetaymin)/(double)thetayfine/tprim;
			}
		}
	}

	cudaFree(result_dev);
	cudaFree(trajec_dev);
	cudaFree(const_dev);
	delete[] result;

Error:
	return cudaStatus;

} // end calc2
 
cudaError_t calculator3(double *consts, double *trajec, double *out1, double *out2, double *out3, double *out4, double *out5, double *out6, double *out7, double *out8) {
	cudaError_t cudaStatus;
	//double *consts = new double[9];
	
	double thetaxmin = consts[3];
	double thetaxmax = consts[4];
	double thetaymin = consts[5];
	double thetaymax = consts[6];

	
	
	



	double *const_dev;
	double *result_dev;
	double *trajec_dev;

	cudaStatus = cudaMalloc((void**)&const_dev, 9 * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&result_dev, 8*thetaxfine*thetayfine*nfreq * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&trajec_dev, 8 * pointsprim * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(const_dev, consts, 9 * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 1 failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(trajec_dev, trajec, 8 * pointsprim * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 2 failed!");
		goto Error;
	}



	myradKernelspin << <BLOCKS, THREADS >> > (const_dev, trajec_dev, result_dev);

	cudaStatus = cudaPeekAtLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "myKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}


	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching myKernel!\n", cudaStatus);
		goto Error;
	}


	double *result = new double[thetaxfine*thetayfine*nfreq*8]();

	cudaStatus = cudaMemcpy(result, result_dev, 8*thetaxfine*thetayfine*nfreq * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy 1 failed!");
		goto Error;
	}

	
		for (int i = 0; i < nfreq; i++) {
			for (int j = 0; j < thetayfine; j++) {
				for (int l = 0; l < thetaxfine; l++) {
					int id = i * thetaxfine*thetayfine + j * thetaxfine + l;
					out1[i] += result[id] *(thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out2[i] += result[id+ 1*thetaxfine * thetayfine*nfreq] *(thetaxmax - thetaxmin) / ((double)thetaxfine-1.0)*(thetaymax - thetaymin) / ((double)thetayfine-1.0);
					out3[i] += result[id + 2*thetaxfine * thetayfine*nfreq] *(thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out4[i] += result[id + 3 * thetaxfine * thetayfine*nfreq] *(thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out5[i] += result[id + 4 * thetaxfine * thetayfine*nfreq] * (thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out6[i] += result[id + 5 * thetaxfine * thetayfine*nfreq] * (thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out7[i] += result[id + 6 * thetaxfine * thetayfine*nfreq] * (thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
					out8[i] += result[id + 7 * thetaxfine * thetayfine*nfreq] * (thetaxmax - thetaxmin) / ((double)thetaxfine - 1.0)*(thetaymax - thetaymin) / ((double)thetayfine - 1.0);
				}
			}
		}
	

	cudaFree(result_dev);
	cudaFree(trajec_dev);
	cudaFree(const_dev);
	delete[] result;

Error:
	return cudaStatus;

} // end calc3

//
//cudaError_t calculator4(double *trajec, double vorthosqavg, double *out, double phimin, double phimax) {
//	cudaError_t cudaStatus;
//	double *consts = new double[10];
//
//	consts[0] = m;
//	consts[1] = wmin;
//	consts[2] = wmax;
//	consts[3] = phimin;
//	consts[4] = phimax;
//	consts[5] = thetamin;
//	consts[6] = thetamax;
//	consts[7] = E0 / m;
//	consts[8] = (double)points1period;
//	consts[9] = vorthosqavg;
//
//
//
//
//	int nDevices;
//	cudaGetDeviceCount(&nDevices);
//	for (int i = 0; i < nDevices; i++) {
//		cudaDeviceProp prop;
//		cudaGetDeviceProperties(&prop, i);
//		printf("Device Number: %d\n", i);
//		printf("  Device name: %s\n", prop.name);
//		printf("  Memory Clock Rate (KHz): %d\n",
//			prop.memoryClockRate);
//		printf("  Memory Bus Width (bits): %d\n",
//			prop.memoryBusWidth);
//		printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
//			2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);
//	}
//
//
//	cudaStatus = cudaSetDevice(0);
//	//cudaDeviceReset();
//
////int devID = 0;
////int argc;
////char argv[1000];
////devID =findCudaDevice(argc, (const char **)argv);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
//		goto Error;
//	}
//
//
//
//	double *const_dev;
//	double *result_dev;
//	double *trajec_dev;
//
//	cudaStatus = cudaMalloc((void**)&const_dev, 10 * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMalloc((void**)&result_dev, harmonicsfine*phifine*nfreq * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMalloc((void**)&trajec_dev, 8 * points1period * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMemcpy(const_dev, consts, 10 * sizeof(double), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy 1 failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMemcpy(trajec_dev, trajec, 8 * points1period * sizeof(double), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy 2 failed!");
//		goto Error;
//	}
//
//
//
//	myradKernelspin << <BLOCKS2, THREADS2 >> > (const_dev, trajec_dev, result_dev);
//
//	cudaStatus = cudaPeekAtLastError();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "myKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//		goto Error;
//	}
//
//
//	cudaStatus = cudaDeviceSynchronize();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching myKernel!\n", cudaStatus);
//		goto Error;
//	}
//
//
//	double *result = new double[harmonicsfine*phifine*nfreq]();
//
//	cudaStatus = cudaMemcpy(result, result_dev, harmonicsfine*phifine*nfreq * sizeof(double), cudaMemcpyDeviceToHost);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy 1 failed!");
//		goto Error;
//	}
//
//	for (int i = 0; i < nfreq; i++) {
//		for (int j = 0; j < harmonicsfine; j++) {
//			for (int l = 0; l < phifine; l++) {
//				int id = i * harmonicsfine*phifine + j * phifine + l;
//				out[i] += result[id] * e*e*(phimax - phimin) / (double)phifine;
//			}
//		}
//	}
//
//	cudaFree(result_dev);
//	cudaFree(trajec_dev);
//	cudaFree(const_dev);
//
//	delete[] result;
//
//Error:
//	return cudaStatus;
//
//} // end calc2