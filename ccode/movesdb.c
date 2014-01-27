#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include"functions.h"

//function to move SD for beta parameter
void movesdb(int k, int m, int nbetagroup, double *beta, double propsdb, gsl_rng *rand_gen, double mnb, double *sdb,double maxsdb)
{
	double sdbprop;
	double accorig,accprop,acc;
	
	sdbprop=sdb[index2(k,m,nbetagroup)]+gsl_ran_flat(rand_gen,-propsdb,propsdb);
	
	if(sdbprop>0.0 && sdbprop < maxsdb)
	{
		accorig=0.0; accprop=0.0;
	
		//adjust for beta priors
		accorig+=log(gsl_ran_gaussian_pdf((beta[index2(k,m,nbetagroup)]-mnb),sdb[index2(k,m,nbetagroup)]));
		accprop+=log(gsl_ran_gaussian_pdf((beta[index2(k,m,nbetagroup)]-mnb),sdbprop));
		//priors and proposals cancel for SD terms
			
		//calculate acceptance
		acc=accprop-accorig;
		acc=exp(acc);
		if(gsl_rng_uniform_pos(rand_gen)<acc) sdb[index2(k,m,nbetagroup)]=sdbprop;
	}	
	return;
}
