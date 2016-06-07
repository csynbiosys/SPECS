/* Simulation of chemotactic cells
 * This file contains the dynamic model
 *
 * F.Tostevin,M.Flores 2010,2011
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "cell.h"
#include <iostream>
#include <cfloat>

//definitions
const int RUN=1, TUMBLE=0;
const int X=0, Y=1, Z=2;

void cell::init(
		void
){
	//setup random number generator
	gsl_rng_env_setup();
	const gsl_rng_type * T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(cell::r,time(NULL));
	t=0;
}
cell::cell(
	void
){
	if(population == 0) {
		init();
	}
	pick_random_direction();
	pos[X] = pos[Y] = pos[Z] = 0.0;
	++population;
}
cell::cell(
	ligandfield * L_in
){
	if(population == 0) {
		init();
	}
	pick_random_direction();
	L=L_in;
	pos[X] = pos[Y] = pos[Z] = 0.0;

	++population;
}

void cell::pick_random_direction(
  void
){
	phi=2*M_PI*gsl_rng_uniform(r);
	double z=2.*gsl_rng_uniform(r)-1.;
	theta=acos(z);
}

void cell::derive_params(
	void
){
	kB=kR*(1.-a0)/a0;
	k=y0*kZ/(a0*(1-y0));
	k_run=(1./b0-1.)/(t_tumble*pow(y0,H));
}
void cell::adapt(
	void
){
	derive_params();
	Ll=localligand();
	a=a0;
	double eps0 =log(1./a0 - 1.)/N;
	m=m0 - (1./alpha)*(eps0 + log(( 1. +  Ll/KA )/( 1. + Ll/KI )));
	y=y0;
	state=(gsl_rng_uniform(r)<b0)?RUN:TUMBLE;
}

double cell::localligand(
	void
){
	return (*L).at(pos,t);
}
int cell::propagate(
	const double dt, const double ds
){
	Ll=localligand();

	//Iterate the dynamical intracellular model
	e=alpha*(m0-m)-log((1.+Ll/KA)/(1.+Ll/KI));
	a=1./(1.+exp(N*e));
	m+=dt*(kR*(1.0-a)-kB*a)+gsl_ran_gaussian(r,sqrt(dt*(kR*(1.0-a)+kB*a)*mNoise));
	if(m<0.) {
		m=0.;
	}
	y+=dt*(k*a*(1.0-y)-kZ*y);
	if(y<0.) {
		y=0.;
	}
	if(y>1.) {
		y=1.;
	}

	//swimming/tumbling dynamics
	if(state==RUN) {
		double costh=cos(theta),
	         cosphi=cos(phi),
	         sinth=sin(theta),
	         sinphi=sin(phi);
		pos[X]+=ds*sinth*cosphi;
		pos[Y]+=ds*sinth*sinphi;
		pos[Z]+=ds*costh;

		double dtheta=fabs(gsl_ran_gaussian(r,sqrt(2.*dt*Dtheta)));
		double cosdth=cos(dtheta),
		       sindth=sin(dtheta);
		double dphi=2.*M_PI*gsl_rng_uniform(r);
		double cosdphi=cos(dphi),
		       sindphi=sin(dphi);
		theta=acos(costh*cosdth-sinth*sindth*cosdphi);
		phi=atan2(sinphi*(costh*sindth*cosdphi+sinth*cosdth)+cosphi*sindth*sindphi,
		          costh*cosphi*sindth*cosdphi-sinphi*sindth*sindphi+sinth*cosphi*cosdth);

		if(gsl_rng_uniform(r)<dt*k_run*std::pow(y,H)) {
			state=TUMBLE;
			return -1;
		}
	} else {
		if(gsl_rng_uniform(r)<dt/t_tumble) {
			state=RUN;
			pick_random_direction();
			return 1;
		}
	}
	return 0;
}
void cell::print_internal_state(
	std::ofstream& outfile
){
	outfile << pos[X] << " "
	        << pos[Y] << " "
	        << pos[Z] << " "
	        << theta << " "
	        << phi << " "
	        << a << " "
	        << m << " "
	        << y << " "
	        << state << " "
	        << std::endl;
}

void cell::load_params(double *params) {
	KI = params[0];
	KA = params[1];
	alpha = params[2];
	m0 = params[3];
	kR = params[4];
	a0 = params[5];
	y0 = params[6];
	kZ = params[7];
	t_tumble = params[8];
	b0 = params[9];
	v0 = params[10];
	Dtheta = params[11];
	mNoise = 1./params[12];
	N = static_cast<int>(params[13]);
	H = static_cast<int>(params[14]);
}

ligandfield::ligandfield(
	void
){
	lengthscale=0.;
	baselevel=0.;
	amplitude=0.;
	frequency=0.;
}

void ligandfield::load_params(double *params) {
	baselevel = params[19];
	lengthscale = params[20];
}
double ligandfield::at(
	const double * pos,
	const double time
){
	double val;
	val=baselevel*exp(pos[X]/(1.*lengthscale));
	if(val > DBL_MAX) {
		std::cerr << "Error: maximum ligand reached at t=" << time << std::endl;
		exit(EXIT_FAILURE);
	}
	if(val<0.) {
		std::cerr << "Error: negative ligand at " << pos[0] << "," << pos[1] << " t=" << time << std::endl;
		exit(EXIT_FAILURE);
	} else {
		return val;
	}
}

gsl_rng * cell::r;
ligandfield * cell::L;
double cell::t;
double cell::KI, cell::KA, cell::alpha, cell::m0,
       cell::kR, cell::kB, 
       cell::a0, cell::y0,
       cell::kZ, cell::k,
       cell::t_tumble,
       cell::b0,
       cell::k_run,
       cell::v0,
       cell::Dtheta,
       cell::mNoise;
int cell::N,
    cell::H;
int cell::population=0;
