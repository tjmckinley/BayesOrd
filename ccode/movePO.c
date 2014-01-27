#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include"functions.h"

//function to move PO parameter
void movePO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double propsdb, double *loglikeorig, gsl_rng *rand_gen, double mnb, double *sdb)
{
	int m;
	double betastore,loglikeprop;
	double accorig,accprop,acc;
	
	betastore=beta[index2(k,0,nbetagroup)];
	beta[index2(k,0,nbetagroup)]+=gsl_ran_gaussian(rand_gen,propsdb);
	for(m=1;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=beta[index2(k,0,nbetagroup)];
	
	//calculate log-likelihood
	loglikeprop=loglikelihood(beta,theta,psi,variables,n,ntheta,nbetagroup);
	//work out acceptance probability  				
	if(isfinite(loglikeprop)!=0)
	{
		accorig=(*loglikeorig);
		accprop=loglikeprop;
		//adjust for beta priors
		accorig+=log(gsl_ran_gaussian_pdf((betastore-mnb),sdb[index2(k,0,nbetagroup)]));
		accprop+=log(gsl_ran_gaussian_pdf((beta[index2(k,0,nbetagroup)]-mnb),sdb[index2(k,0,nbetagroup)]));
		//proposals cancel
		acc=accprop-accorig;
		acc=exp(acc);
		if(gsl_rng_uniform_pos(rand_gen)<acc) (*loglikeorig)=loglikeprop;
		else
		{
			for(m=0;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=betastore;
		}
	}
	else
	{
		for(m=0;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=betastore;
	}
	return;
}
