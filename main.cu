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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "constants.h"
#include <random>

using namespace std;


int main() {
	
	char a[150];
	char b[150];
	ofstream enh1,enh2,enh3,file;
	strcpy_s(a, "spectrum.txt");
	enh1.open(a);
	
	int div = 500;
	file.open("potential.txt");
	for (int i = 0; i < div; i++) {
		double x = (double)i / (double)div*dp;
		file << x << " " << UDT(x) << endl;
	}
	file.close();

	cudaError_t cudaStatus;
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		/*printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n",
			prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",
			prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
			2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);*/
	}


	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	

	

	
		double y[7];
		
		double xpos, ypos, z, vx, vy, vz, gam, t, ok;

		t = 0.0;
		vx = eta/1.0/gam0;
		vy = 0.0;
		xpos = 0.0;

		y[0] = vx;
		y[1] = vy;
		y[2] = (-pow(vx, 2.0) - pow(vy, 2.0)) / 2.0;
		y[3] = xpos;
		y[4] = ypos;
		y[5] = 0.0;
		y[6] = 0.0;
		gam = gam0;

		double *omega = new double[freqiters]();
		double *dEdtspinupup1 = new double[wfine]();
		double *dEdtspindowndown1 = new double[wfine]();
		double *dEdtspinupdown1 = new double[wfine]();
		double *dEdtspindownup1 = new double[wfine]();
		double *dEdtspinupup2 = new double[wfine]();
		double *dEdtspindowndown2 = new double[wfine]();
		double *dEdtspinupdown2 = new double[wfine]();
		double *dEdtspindownup2 = new double[wfine]();
		double *dWdw = new double[wfine]();
		double *dWdwspin = new double[wfine]();
		double *yrad = new double[7];
		double *radresult = new double[thetaxfine*thetayfine*wfine]();
		double *thetas = new double[4]();

		yrad[0] = y[0];
		yrad[1] = y[1];
		yrad[2] = y[2];
		yrad[3] = y[3];
		yrad[4] = y[4];
		yrad[5] = y[5];
		yrad[6] = y[6];

		int curN = 0;
		int j = 0;

		double l_f = pi*N/w0;
		//double l_f = 1.0*L0;


		W(yrad, gam, omega, dEdtspinupup1, dEdtspindowndown1, dEdtspinupdown1, dEdtspindownup1, dEdtspinupup2, dEdtspindowndown2, dEdtspinupdown2, dEdtspindownup2, radresult, thetas, l_f, curN, j, t);
		cudaDeviceSynchronize();
		double *dWdwavg1 = new double[wfine]();
		double *dWdwupup1 = new double[wfine]();
		double *dWdwdowndown1 = new double[wfine]();
		double *dWdwupdown1 = new double[wfine]();
		double *dWdwdownup1 = new double[wfine]();

		double *dWdwupup2 = new double[wfine]();
		double *dWdwdowndown2 = new double[wfine]();
		double *dWdwupdown2 = new double[wfine]();
		double *dWdwdownup2 = new double[wfine]();
		


		for (int l = 0; l < wfine; l++) {
			dWdwupup1[l] = dEdtspinupup1[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwdowndown1[l] = dEdtspindowndown1[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwupdown1[l] = dEdtspinupdown1[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwdownup1[l] = dEdtspindownup1[l] * alfa / (4.0*pi*pi) / omega[l];

			dWdwupup2[l] = dEdtspinupup2[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwdowndown2[l] = dEdtspindowndown2[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwupdown2[l] = dEdtspinupdown2[l] * alfa / (4.0*pi*pi) / omega[l];
			dWdwdownup2[l] = dEdtspindownup2[l] * alfa / (4.0*pi*pi) / omega[l];

			dWdwavg1[l] = 0.5*(dWdwupup1[l] + dWdwdowndown1[l] + dWdwupdown1[l] + dWdwdownup1[l]);
			enh1 << omega[l] << " " << dWdwupup1[l] * omega[l] << " " << dWdwdowndown1[l] * omega[l] << " " << dWdwupdown1[l] * omega[l] << " " << dWdwdownup1[l] * omega[l] << " " << dWdwupup2[l] * omega[l] << " " << dWdwdowndown2[l] * omega[l] << " " << dWdwupdown2[l] * omega[l] << " " << dWdwdownup2[l] * omega[l] << endl;

		}

		
				
		
		
		delete[] omega;
		delete[] radresult;
		delete[] thetas;
		delete[] dEdtspinupup1;
		delete[] dEdtspindowndown1;
		delete[] dEdtspinupdown1;
		delete[] dEdtspindownup1;
		delete[] dEdtspinupup2;
		delete[] dEdtspindowndown2;
		delete[] dEdtspinupdown2;
		delete[] dEdtspindownup2;
		delete[] yrad;

		
		delete[] dWdwupup2;
		delete[] dWdwdowndown2;
		delete[] dWdwupdown2;
		delete[] dWdwdownup2;
		
	enh1.close();



Error:
	return cudaStatus;

}