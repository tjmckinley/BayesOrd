#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include"functions.h"

//function to move PO parameter to NPO parameter
void movePOtoNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, gsl_rng *rand_gen, double mnb, double *sdb, double sdt, double maxsdb, int fixed)
{
	int m,valid;
	double *betacondstore = (double *)malloc((ntheta-1)*sizeof(double));
	double *u = (double *)malloc((ntheta-1)*sizeof(double));
	double betastore,loglikeprop,temp,sdbstore,lower=0.0,upper=0.0;
	double accorig,accprop,acc;
		
	//propose move from PO to NPO model
	for(m=0;m<(ntheta-1);m++) betacondstore[m]=betacond[m];
	
	//first level is current value
	betastore=beta[index2(k,0,nbetagroup)];
	
	//update beta parameters
	beta[index2(k,0,nbetagroup)]=0.0;
	for(m=0;m<(ntheta-1);m++)
	{
		u[m]=gsl_ran_flat(rand_gen,-propsdb,propsdb);
		beta[index2(k,0,nbetagroup)]+=u[m];
	}
	beta[index2(k,0,nbetagroup)]+=betastore*(4.0-(ntheta+1));
	beta[index2(k,0,nbetagroup)]=(ntheta/2.0)*beta[index2(k,0,nbetagroup)];
	for(m=1;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=(ntheta/2.0)*(betastore-u[m-1]);
	
	//update SD parameters
	sdbstore=sdb[index2(k,0,nbetagroup)];
	if(fixed==0)
	{
		lower=((sdbstore-propsdb)<0.0 ? 0.0:(sdbstore-propsdb));		
		upper=((sdbstore+propsdb)>maxsdb ? maxsdb:(sdbstore+propsdb));
		for(m=1;m<ntheta;m++) sdb[index2(k,m,nbetagroup)]=gsl_ran_flat(rand_gen,lower,upper);
	}
	
	//update validity conditions
	m=0; valid=1;
	while(m<(ntheta-1) && valid==1)
	{
		temp=beta[index2(k,m,nbetagroup)]-beta[index2(k,m+1,nbetagroup)];
		betacond[m]+=((minvar[k]*temp)<(maxvar[k]*temp) ? (minvar[k]*temp):(maxvar[k]*temp));
		if((theta[m]-theta[m+1])>betacond[m]) valid=0;
		m++;
	}
	if(valid==1)
	{
		//calculate log-likelihood
		loglikeprop=loglikelihood(beta,theta,psi,variables,n,ntheta,nbetagroup);
		//work out acceptance probability  				
		if(isfinite(loglikeprop)!=0)
		{
			accorig=(*loglikeorig);
			accprop=loglikeprop;
		
			//adjust for beta priors
			accorig+=log(gsl_ran_gaussian_pdf((betastore-mnb),sdbstore));
			for(m=0;m<ntheta;m++) accprop+=log(gsl_ran_gaussian_pdf((beta[index2(k,m,nbetagroup)]-mnb),sdb[index2(k,m,nbetagroup)]));
			
			if(fixed==0)
			{
				//adjust for SD priors
				accorig-=log(maxsdb);
				accprop-=ntheta*log(maxsdb);
			}
		
			//adjust for theta priors conditional on beta
			for(m=1;m<ntheta;m++)
			{
				accorig-=log(1.0-gsl_cdf_gaussian_P(theta[m-1]-betacondstore[m-1],sdt));
				accprop-=log(1.0-gsl_cdf_gaussian_P(theta[m-1]-betacond[m-1],sdt));
			}
		
			//adjust for proposals
			accorig-=(ntheta-1)*log(2*propsdb);
			if(fixed==0) accorig-=(ntheta-1)*log(upper-lower);
		
			//equal probs of move proposal cancel
			
			//adjust for Jacobian
			accprop+=log(ntheta)+(ntheta-1.0)*log(ntheta/2.0);
		
			//calculate acceptance
			acc=accprop-accorig;
			acc=exp(acc);
			if(gsl_rng_uniform_pos(rand_gen)<acc)
			{
				(*loglikeorig)=loglikeprop;
				//set beta status
				betastatus[k]=1;
			}
			else
			{
				for(m=0;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=betastore;
				if(fixed==0){for(m=0;m<ntheta;m++) sdb[index2(k,m,nbetagroup)]=sdbstore;}
				for(m=0;m<(ntheta-1);m++) betacond[m]=betacondstore[m];
			}
		}
		else
		{
			for(m=0;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=betastore;
			if(fixed==0){for(m=0;m<ntheta;m++) sdb[index2(k,m,nbetagroup)]=sdbstore;}
			for(m=0;m<(ntheta-1);m++) betacond[m]=betacondstore[m];
		}
	}
	else
	{
		for(m=0;m<ntheta;m++) beta[index2(k,m,nbetagroup)]=betastore;
		if(fixed==0){for(m=0;m<ntheta;m++) sdb[index2(k,m,nbetagroup)]=sdbstore;}
		for(m=0;m<(ntheta-1);m++) betacond[m]=betacondstore[m];
	}
	//free memory from the heap
	free(betacondstore);free(u);
	return;
}
