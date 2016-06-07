/* Simulation of chemotactic cells
 * In this file, simulations are managed and output controlled
 *
 * F.Tostevin,M.Flores 2010,2011
 */

/******************************************************************
constants loaded from file constants via load_params() function.
modified 12 July 2010 by marlo flores.
******************************************************************/
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "load_params.h"
#include "cell.h"

//Definitions
const int RUN=1, TUMBLE=0;
const int X=0, Y=1, Z=2;

using namespace std;
void load_program_params(double *);

//program constants
double dt;			//time resolution
int nMax;				//number of iterations
int nCells;			//number of cells
int aveRsmpl;		//intervals between sampling ensemble averages (in time steps)
double ds;

int main (
		int argc,
		char *argv[]
) {
	ofstream outfile("output.dat");
	outfile << "t "
	        << "<x> "
	        << "<y> "
	        << "<z> "
	        << "<[L]> "
	        << "<e> "
	        << "<a> "
	        << "<m> "
	        << "<Yp> "
	        << "bias "
	        << endl;
	
	//read constants from file
	const int numparams = 21;
	double params[numparams];
	load_params("constants",params,numparams);
	
	//load program constants
	load_program_params(params);

	cell C[nCells];
	ligandfield L;
	L.load_params(params);
	
	for(int i=0; i<nCells; ++i) {
		C[i]=cell(&L);
		C[i].load_params(params);
		C[i].adapt();
	}
	
	//temporary variables
	unsigned int n=0, i=0;
	
	//ds and dt
	ds = cell::v0*dt;	
	
	double av_x, av_y, av_z, av_e, av_L, av_a, av_m, av_yp, av_r2, bias;
	ofstream paramsout("params.txt");
	paramsout << "dt" << '\t' << dt << endl
	       << "tmax" << '\t' << nMax*dt << endl
	       << "nCells" << '\t' << nCells << endl << endl
	       << "receptor KI" << '\t' << C[0].KI << endl
	       << "receptor KA" << '\t' << C[0].KA << endl << endl
	       << "meth alpha" << '\t' << C[0].alpha << endl
	       << "meth m0" << '\t' << C[0].m0 << endl << endl
	       << "meth kR" << '\t' << C[0].kR << endl
	       << "meth kB" << '\t' << C[0].kB << endl
	       << "adapted a" << '\t' << C[0].a0 << endl
	       << "adapted y" << '\t' << C[0].y0 << endl
	       << "CheY kZ" << '\t' << C[0].kZ << endl << endl
	       << "CheY k" << '\t' << C[0].k << endl
	       << "motor t_tumble" << '\t' << C[0].t_tumble << endl
	       << "adabpted bias" << '\t' << C[0].b0 << endl
	       << "motor k_run" << '\t' << C[0].k_run << endl
	       << "speed" << '\t' << C[0].v0 << endl
	       << "rot diff" << '\t' << C[0].Dtheta << endl << endl
	       << "noise m" << '\t' << C[0].mNoise << endl
	       << "receptor N" << '\t' << C[0].N << endl
	       << "motor H" << '\t' << C[0].H << endl << endl
	       << "ligand L0" << '\t' << L.baselevel << endl
	       << "ligand l" << '\t' << L.lengthscale << endl << endl;
				
	//Propagate the simulation
	while(n<nMax){
		for(int i=0; i<nCells; ++i) {
			int j=C[i].propagate(dt, ds);
		}

		cell::t = n*dt;

		if(n%aveRsmpl==0) {
			av_x=0.;
			av_y=0.;
			av_z=0.;
			av_L=0.;
			av_e=0.;
			av_a=0.;
			av_m=0.;
			av_yp=0.;
			bias=0.;
			//calculate ensemble averages
			for(int i=0; i<nCells; ++i) {
				av_L+=C[i].localligand()/nCells;
				av_e+=C[i].e/nCells;
				av_a+=C[i].a/nCells;
				av_x+=C[i].pos[X]/nCells;
				av_y+=C[i].pos[Y]/nCells;
				av_z+=C[i].pos[Z]/nCells;
				av_m+=C[i].m/nCells;
				av_yp+=C[i].y/nCells;
				bias+=C[i].state/(1.*nCells);
			}
		
			outfile << n*dt << " "
			        << av_x << " " 
			        << av_y << " "
			        << av_z << " "
			        << av_L << " "
			        << av_e << " "
			        << av_a << " "
			        << av_m << " "
			        << av_yp << " "
			        << bias << " "
			        << endl;
		}
		
		++n;
	}

	ofstream finalout; finalout.open("final.dat");
	for(int i=0; i<nCells; ++i) {
		finalout << n << " ";
		C[i].print_internal_state(finalout);
	}
	finalout.close();

	return EXIT_SUCCESS;	
}

void load_program_params(double *params) {
	dt = params[15];
	nMax = static_cast<int>(params[16]);
	nCells = static_cast<int>(params[17]);
	aveRsmpl = static_cast<int>(params[18]);
}
