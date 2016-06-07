/* Simulation of chemotactic cells
 * Definitions of paramters and functions of the model
 *
 * F.Tostevin,M.Flores 2010,2011
 */

#ifndef CELL_H
#define CELL_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "load_params.h"

class ligandfield{
public:
	double lengthscale, baselevel, amplitude, frequency;

	ligandfield();
	void load_params(double *p);
	double at(const double *, const double);
};

class cell{
public:
	//If you want different parameters in each cell, remove the static  definition
	//and add a way to assign parameters to each cell individually.
	static double KI, KA, alpha, m0,	//parameters of a(L,m)
	              kR, kB,							//kinetic (de-)methylation rates
	              a0, y0,							//adapted activity
	              kZ, k,							//kinetic yp-modification rates
	              t_tumble,						//mean tumble duration
	              b0,									//adapted run-bias
	              k_run,							//tumble rate parameter
	              v0,									//swimming speed
	              Dtheta,							//rotational diffusion constant
	              mNoise;							//noise strength
	static int N,											//receptor cluster factor
	           H;											//yp-cooperativity
	static double t;
	static gsl_rng * r;
	static ligandfield * L;

	//Cell state variables
	double pos[3];
	double e;
	double a;
	double m;
	double y;
	double Ll;
	double theta;
	double phi;
	bool state;
	cell();
	cell(ligandfield*);
	
	void init();
	void derive_params();
	void pick_random_direction();
	static void load_params(double *p);
	int propagate(const double, const double);
	double localligand(); 
	void adapt();
	void print_internal_state(std::ofstream&);
	
private:
	static int population;
};

#endif
