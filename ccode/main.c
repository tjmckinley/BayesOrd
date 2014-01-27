/*Metropolis-Hastings routine for fitting ordinal regression model*/
/*Written for running multiple chains with multicore in R*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//random number generators from GSL library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include"functions.h"
 
//main MCMC function
int main(int argc, char * argv[])
{
	//check number of arguments
	if(argc<3)
	{
		printf("Wrong number of arguments being passed to function!\n");
		return 0;
	}
	// initialize random seed 		
	const gsl_rng_type * rand_type;
	gsl_rng_env_setup();
	rand_type = gsl_rng_default;
	gsl_rng * rand_gen = gsl_rng_alloc (rand_type);
	//gsl_rng_set(rand_gen,(int)(time(NULL)));
	gsl_rng_set(rand_gen,atoi(argv[2]));
	       
	if(atoi(argv[1])==1) printf ("\n#### GSL_generator ####\ngenerator type: %s\n", gsl_rng_name (rand_gen));
	//printf ("seed = %lu\n", gsl_rng_default_seed);
	//printf ("first value = %lu\n", gsl_rng_get (rand_gen));
		
	int i,j,k,m;
	double u;
  		  	
	//import initial conditions and priors into function
	//and set up parameters to be estimated
	int n,niter,nvariables,nbeta,nbetagroup,ntheta,nsavecoda,npsi,RE,movetype,varselect,fixed,maxinteraction,runtraining,nitertrain;
	double mnb,varb,maxsdb,vart,propsdb,propsdt,mnpsi,shvarp,rtvarp,propsdp,propsdvarp;
	
	//read in priors and initial values
	char dummy[40];

	FILE *PriorFile;
	PriorFile = fopen("priors.txt","r");
	if(!PriorFile)
	{
		printf("Could not open priors file\n");
		return 1;
	}
	j=fscanf(PriorFile, "%d", &n);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &niter);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &nsavecoda);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &nvariables);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &nbeta);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &nbetagroup);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &ntheta);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &mnb);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &varb);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &maxsdb);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &fixed);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &vart);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &propsdb);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &propsdt);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &mnpsi);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &shvarp);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &rtvarp);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &propsdp);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%lf", &propsdvarp);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &npsi);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &RE);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &movetype);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &varselect);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &maxinteraction);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &runtraining);j=fscanf(PriorFile, "%s", dummy);
	j=fscanf(PriorFile, "%d", &nitertrain);j=fscanf(PriorFile, "%s", dummy);
	fclose(PriorFile);
	
	//set up parameter to control printing in multiple chains run in parallel
	int printall=1;
	if(atoi(argv[1])==1) printall=1;
	else printall=0;
	
	if(printall==1)
	{	
		//visual check of prior info
		printf("\n#### Data ####\n");
		printf("number of data points=%d\n",n);
		printf("number of iterations=%d\n",niter);
		printf("Save coda after every %d iterations\n",nsavecoda);
	
		printf("\n#### Parameters ####\n");
		printf("number of variables=%d\n",nvariables);
		printf("number of beta parameters overall=%d\n",nbeta);
		printf("number of beta parameters in each group=%d\n",nbetagroup);
		printf("number of theta parameters=%d\n",ntheta);
	
		printf("\n#### Priors ####\n");
		printf("mean beta=%f\n",mnb);
		if(fixed==1) printf("var beta=%f\n",varb);
		else printf("SD beta ~ N(0,sigma^2), where sigma ~ U(0,%f)\n",maxsdb);
		printf("var theta=%f\n",vart);
		if(RE==1)
		{
			printf("mn psi=%f\n",mnpsi);
			printf("shape varp=%f\n",shvarp);
			printf("rate varp=%f\n",rtvarp);
			printf("mean varp=%f\tvariance varp=%f\n",shvarp/rtvarp,shvarp/(rtvarp*rtvarp));
		}
	
		printf("\n#### Proposal ####\n");
		printf("proposal SD for beta=%f\n",propsdb);
		printf("proposal SD for theta=%f\n",propsdt);
		if(RE==1)
		{
			printf("proposal SD for psi=%f\n",propsdp);
			printf("proposal SD for varp=%f\n",propsdvarp);
		}
	
		if(movetype==0) printf("\n#### PROPORTIONAL ODDS MODEL ####\n");
		else
		{
			if(movetype==1) printf("\n#### NON-PROPORTIONAL ODDS MODEL ####\n");
			else printf("\n#### PROPORTIONAL ODDS OR NON PROPORTIONAL ODDS MODEL ####\n");
		}
		if(RE==0) printf("Model includes NO random effect terms\n");
		else printf("Model includes RANDOM EFFECT terms\n");
		if(varselect==1) printf("Model includes VARIABLE SELECTION\n");
		if(runtraining==1) printf("\nTraining run being used for %d iterations\n",nitertrain);
	}
	
	//import data
	double *variables = (double *) malloc(n*(nbetagroup+3)*sizeof(double));
	FILE *DataFile;
	DataFile = fopen("data.dat","r");
	if(!DataFile)
	{
		printf("Could not open data file\n");
		return 1;
	}
	for(i=0;i<n;i++) for(j=0;j<(nbetagroup+3);j++) k=fscanf(DataFile,"%lf",&variables[index2(i,j,n)]);
	fclose(DataFile);
	
/*	//visual check of data*/
/*	printf("\n#### Visual check of data ####\n");*/
/*	for(i=0;i<5;i++)*/
/*	{*/
/*		for(j=0;j<(nbetagroup+2);j++) printf("%f\t",variables[index2(i,j,n)]);*/
/*		if(RE==1) printf("%f",variables[index2(i,j,n)]);*/
/*		printf("\n");*/
/*	}*/
		
	//declare arrays and variables for MCMC
	int npars=nbeta+ntheta+npsi+1+nbetagroup;
	if(fixed==0) npars+=nbeta;
	
	int kcoda=0;
	double *posterior = (double *) malloc((npars+1)*nsavecoda*sizeof(double));
	for(j=0;j<((npars+1)*nsavecoda);j++) posterior[j]=0;
	double *beta = (double *) malloc(nbeta*sizeof(double));
	double *theta = (double *) malloc(ntheta*sizeof(double));
	double *psi = (double *) malloc(npsi*sizeof(double));
	double varp=1.0;
	//set up conditions for beta parameters
	double *betacond = (double *)malloc((ntheta-1)*sizeof(double));
	double *betacondstore_v = (double *)malloc((ntheta-1)*sizeof(double));
	for(i=0;i<(ntheta-1);i++)
	{
		betacond[i]=0.0;
		betacondstore_v[i]=0.0;
	}
	//set minimum and maximum values for each variable based on the data
	double *minvar = (double *) malloc(nbetagroup*sizeof(double));
	double *maxvar = (double *) malloc(nbetagroup*sizeof(double));
	for(j=0;j<nbetagroup;j++)
	{
		minvar[j]=variables[index2(0,j,n)];
		maxvar[j]=variables[index2(0,j,n)];
		for(i=1;i<n;i++)
		{
			minvar[j]=(minvar[j]<variables[index2(i,j,n)] ? minvar[j]:variables[index2(i,j,n)]);
			maxvar[j]=(maxvar[j]>variables[index2(i,j,n)] ? maxvar[j]:variables[index2(i,j,n)]);
		}
	}
	for(i=0;i<nbetagroup;i++)
	{
		if(minvar[i]==maxvar[i])
		{
			printf("No variation in %dth variable\n",i);
			return 1;
		}
	}
/*	for(i=0;i<nbetagroup;i++) printf("min[%d]=%f\tmax[%d]=%f\n",i,minvar[i],i,maxvar[i]);*/
/*	printf("\n");*/
	double temp=0.0;
	
	//set up temporary storage variables
	double psistore,varpstore;
	int valid;
	//set up vector for shuffling when updating beta terms
	int *betashuffle = (int *)malloc(nvariables*sizeof(int));
	for(i=0;i<nvariables;i++) betashuffle[i]=i;
	int *betashufflepar = (int *)malloc(nbetagroup*sizeof(int));
	for(i=0;i<nbetagroup;i++) betashufflepar[i]=i;
	int *betashuffletheta = (int *)malloc(ntheta*sizeof(int));
	for(i=0;i<ntheta;i++) betashuffletheta[i]=i;
	//set up vector for recording parameter status (0=PO, 1=NPO)
	int *betastatus = (int *)malloc(nbetagroup*sizeof(int));
	for(i=0;i<nbetagroup;i++) betastatus[i]=1;
	//set up vector for recording variable status (0=exc, 1=inc)
	int *betastatusvar = (int *)malloc(nvariables*sizeof(int));
	for(i=0;i<nvariables;i++) betastatusvar[i]=1;
	//import assignment information
	int *xassign = (int *)malloc((nvariables+1)*sizeof(int));
	DataFile = fopen("xassign.txt","r");
	if(!DataFile)
	{
		printf("Could not open xassign file\n");
		return 1;
	}
	for(i=0;i<(nvariables+1);i++) k=fscanf(DataFile,"%d",&xassign[i]);
	fclose(DataFile);
	//visual check
/*	printf("xassign\n");*/
/*	for(i =0;i<(nvariables+1);i++) printf("%d\t",xassign[i]);*/
/*	printf("\n");*/
		
	int *intfactor = (int *)malloc(nvariables*nvariables*sizeof(int));
	int *interaction = (int *)malloc(nvariables*sizeof(int));
	int *intstart = (int *)malloc((maxinteraction+2)*sizeof(int));
	for(i=0;i<nvariables;i++)
	{
		interaction[i]=0.0;
		for(j=0;j<nvariables;j++) intfactor[index2(i,j,nvariables)]=0.0;
	}
	if(maxinteraction>0)
	{
		//read in interaction information if necessary
		DataFile = fopen("intfactor.txt","r");
		if(!DataFile)
		{
			printf("Could not open intfactor file\n");
			return 1;
		}
		for(i=0;i<nvariables;i++) for(j=0;j<nvariables;j++) k=fscanf(DataFile,"%d",&intfactor[index2(i,j,nvariables)]);
		fclose(DataFile);
		
		DataFile = fopen("interaction.txt","r");
		if(!DataFile)
		{
			printf("Could not open interaction file\n");
			return 1;
		}
		for(i=0;i<nvariables;i++) k=fscanf(DataFile,"%d",&interaction[i]);
		fclose(DataFile);
		
		DataFile = fopen("intstart.txt","r");
		if(!DataFile)
		{
			printf("Could not open intstart file\n");
			return 1;
		}
		for(i=0;i<(maxinteraction+2);i++) k=fscanf(DataFile,"%d",&intstart[i]);
		fclose(DataFile);
		
		//visual check of data inputs
		
		if(printall==1)
		{
			printf("intfactor:\n");
			for(i=0;i<nvariables;i++)
			{
				for(j=0;j<nvariables;j++) printf("%d ",intfactor[index2(i,j,nvariables)]);
				printf("\n");
			}
			printf("\ninteraction:\n");
			for(i=0;i<nvariables;i++) printf("%d ",interaction[i]);
			printf("\n");
			printf("\nintstart:\n");
			for(i=0;i<(maxinteraction+2);i++) printf("%d ",intstart[i]);
			printf("\n");
		}
	}	
			
	//set up loglikelihood stuff
	double loglikeorig,loglikeprop,accorig,accprop,acc;
	
	double *sdb=(double *)malloc(nbeta*sizeof(double));
	for(i=0;i<nbeta;i++) sdb[i]=sqrt(varb);
	
	double sdt=sqrt(vart);
	double sdp;
	//sample initial values
	m=0;
	loglikeorig=0.0/0.0;
	while(isfinite(loglikeorig)==0 && m<1000)
	{
		//reset parameters
		for(i=0;i<nbetagroup;i++) betastatus[i]=0;
		for(i=0;i<nbeta;i++) beta[i]=0.0;
		//choose initial variables for inclusion if variable selection being implemented
/*		if(varselect==1)*/
/*		{*/
/*			u=gsl_ran_flat(rand_gen,0.0,1.0);*/
/*			i=(int) gsl_ran_binomial(rand_gen,u,(unsigned int) nbetagroup);*/
/*			i=(i<1 ? 1:i);*/
/*			//choose variables for inclusion at random without replacement*/
/*	 		gsl_ran_shuffle(rand_gen,betashuffle,nbetagroup,sizeof(int));*/
/*			for(j=0;j<i;j++) betastatus[betashuffle[j]]=2;*/
/*		}*/
		//sample values for beta	
		for(i=0;i<nbetagroup;i++)
		{
			if(betastatus[i]!=2) beta[index2(i,0,nbetagroup)]=mnb+gsl_ran_gaussian(rand_gen,1.0);
			else beta[index2(i,0,nbetagroup)]=0.0;
		}
		if(movetype==0 || runtraining == 1)
		{
			for(i=0;i<nbetagroup;i++)
			{
				if(betastatus[i]!=2)
				{
					for(j=1;j<ntheta;j++) beta[index2(i,j,nbetagroup)]=beta[index2(i,j-1,nbetagroup)];
					betastatus[i]=0;
				}
			}
		}
		else
		{
			if(movetype==1)
			{
				for(i=0;i<nbetagroup;i++)
				{
					if(betastatus[i]!=2)
					{
						for(j=1;j<ntheta;j++) beta[index2(i,j,nbetagroup)]=beta[index2(i,j-1,nbetagroup)]-gsl_ran_flat(rand_gen,0.0,1.0);
						betastatus[i]=1;
					}
				}
			}
			else
			{
				for(i=0;i<nbetagroup;i++)
				{
					if(betastatus[i]!=2)
					{
						betastatus[i]=0;//gsl_ran_bernoulli(rand_gen,0.5);
						if(betastatus[i]==0)
						{
							for(j=1;j<ntheta;j++)  beta[index2(i,j,nbetagroup)]=beta[index2(i,j-1,nbetagroup)];
						}
						else
						{
							for(j=1;j<ntheta;j++)  beta[index2(i,j,nbetagroup)]=beta[index2(i,j-1,nbetagroup)]-gsl_ran_flat(rand_gen,0.0,1.0);
						}
					}
				}
			}
		}	
		//set and check prior conditions
		for(i=0;i<(ntheta-1);i++)
		{
			betacond[i]=0.0;
			for(j=0;j<nbetagroup;j++)
			{
				temp=beta[index2(j,i,nbetagroup)]-beta[index2(j,i+1,nbetagroup)];
				betacond[i]+=((minvar[j]*temp)<(maxvar[j]*temp) ? (minvar[j]*temp):(maxvar[j]*temp));
			}
		}
		//randomly sample cut-points
		theta[0]=gsl_ran_gaussian(rand_gen,1.0);
		for(j=1;j<ntheta;j++)
		{
			theta[j]=gsl_ran_gaussian(rand_gen,1.0);
			while(theta[j]<theta[j-1]) theta[j]=gsl_ran_gaussian(rand_gen,1.0);
		}
		//check validity
		valid=1;
		for(j=0;j<(ntheta-1);j++)
		{
			if(theta[j]>=(theta[j+1]+betacond[j])) valid=0;
		}
		//sample hyperparameters
		if(fixed==0)
		{
			for(j=0;j<nbeta;j++) sdb[j]=gsl_ran_flat(rand_gen,0.0,maxsdb);
		}
		if(valid==1)
		{
			//sample random effects
			if(RE==1)
			{
				varp=gsl_ran_gamma(rand_gen,1.0,1.0);
				for(j=0;j<npsi;j++) psi[j]=mnpsi+gsl_ran_gaussian(rand_gen,sqrt(varp));
			}
			else
			{
				varp=1.0;
				for(j=0;j<npsi;j++) psi[j]=0.0;
			}
			//calculate log-likelihood
			loglikeorig=loglikelihood(beta,theta,psi,variables,n,ntheta,nbetagroup);
		}
		m++;
	}
		
	if(m==1000){printf("\nCrap initialisation: chain %s!\n",argv[1]);return 0;}
	printf("\nSystem initialised OK: chain %s!\n",argv[1]);
	
/*	//print theta's to the screen as check*/
/*	printf("thetas:\n");*/
/*	for(j=0;j<ntheta;j++) printf("theta[%d]=%f\n",j,theta[j]);*/
/*	printf("\n");*/
	
	//set up and read in cumulative count vector for efficient psi updates
	int *psicount = (int *)malloc((npsi+1)*sizeof(int));
	if(printall==1) printf("RE=%d\n",RE);
	if(RE==1)
	{
		DataFile = fopen("psicount.txt","r");
		if(!DataFile)
		{
			printf("Could not open psicount file\n");
			return 1;
		}
		psicount[0]=0;
		for(j=1;j<(npsi+1);j++) k=fscanf(DataFile,"%d",&psicount[j]);
		fclose(DataFile);
	
		if(printall==1)
		{
			//visual check 
			for(j=0;j<5;j++) printf("psicount[%d]=%d\n",j,psicount[j]);
			printf("...\n");
			for(j=npsi-5;j<(npsi+1);j++) printf("psicount[%d]=%d\n",j,psicount[j]);
			printf("\n");
		}
	}
	
	//function to check validity of proposals
	if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
	
	//initialise time
	time_t start,end;
	time (&start);
	double runtime;
/*	printf("\n#### Timing started: chain %s ####\n",argv[1]);*/

	double pinteraction;
	//set move probabilities (0: keep, 1: add/remove)
	double pvecvar[2];
	if(varselect==1) pvecvar[0]=0.5;
	else pvecvar[0]=1.0;
	pvecvar[1]=1.0;
	double pvecvarmove[2];
	pvecvarmove[0]=0.0;
	pvecvarmove[1]=1.0;
	
	//set move probabilities (0: move, 1: PO/NPO change)
	double pvec[2];
	if(movetype==0 || movetype==1)
	{
		pvec[0]=1.0;
		pvec[1]=1.0;
	}
	else
	{
		pvec[0]=0.5;
		pvec[1]=1.0;
	}
	//set zero to PO and NPO move probabilities
	double pzero[2];
	if(varselect==1)
	{
		if(movetype==0) pzero[0]=1.0;
		else
		{
			if(movetype==1)	pzero[0]=0.0;
			else pzero[0]=0.5;
		}
		pzero[1]=1.0;
	}
/*	printf("\nSelect probs:\tpvar1=%f\tpvar2=%f\n",pvecvar[0],pvecvar[1]-pvecvar[0]);*/
/*	printf("\nSelect probs:\tp1=%f\tp2=%f\n",pvec[0],pvec[1]-pvec[0]);*/
/*	printf("\nSelect probs:\tpzero1=%f\tpzero2=%f\n",pzero[0],pzero[1]-pzero[0]);*/
/*	*/
/*	printf("betastatus:\n");*/
/*	for(i=0;i<nbetagroup;i++) printf("%d ",betastatus[i]);*/
/*	printf("\n");*/
/*	*/
/*	printf("maxinteraction=%d\n",maxinteraction);*/
		
	//start MCMC chain
	if(runtraining==1)
	{
	 	for(i=1;i<nitertrain;i++)
	 	{
			//update linear terms in random order
	 		gsl_ran_shuffle(rand_gen,betashufflepar,nbetagroup,sizeof(int));
	 		for(j=0;j<nbetagroup;j++)
	 		{
	 			k=betashufflepar[j];
				movePO(k,n,nbetagroup,ntheta,beta,theta,psi,variables,propsdb,&loglikeorig,rand_gen,mnb,sdb);
			}			
			//update theta terms in random order
	 		gsl_ran_shuffle(rand_gen,betashuffletheta,ntheta,sizeof(int));
			for(m=0;m<ntheta;m++) movetheta(betashuffletheta[m],n,nbetagroup,ntheta,beta,theta,psi,variables,betacond,propsdt,&loglikeorig,rand_gen,sdt);
			//now check validity conditions
			if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
			//move SD for beta parameter
			if(fixed==0)
			{
				//update SD terms in random order
	 			gsl_ran_shuffle(rand_gen,betashufflepar,nbetagroup,sizeof(int));
				for(j=0;j<nbetagroup;j++)
				{
					k=betashufflepar[j];
					movesdb(k,0,nbetagroup,beta,propsdb,rand_gen,mnb,sdb,maxsdb);
				}
			}
		}
		if(movetype == 1) for(j = 0; j < nbetagroup; j++) betastatus[j] = 1;
		printf("Finished training run, restarting chains...\n");
	}
	
	FILE *fileres;
	char filename[20];
	sprintf(filename,"%s%s%s","codaMCMC",argv[1],".txt");
	
	fileres=fopen(filename,"w");
	for(j=0;j<nbeta;j++) fprintf(fileres,"beta%d\t",j);
	for(j=0;j<ntheta;j++) fprintf(fileres,"theta%d\t",j);
	for(j=0;j<nbetagroup;j++) fprintf(fileres,"betastatus%d\t",j);
	for(j=0;j<npsi;j++) fprintf(fileres,"psi%d\t",j);
	if(fixed==0)
	{
		for(j=0;j<nbeta;j++) fprintf(fileres,"sdb%d\t",j);
	}
	fprintf(fileres,"varp\t");
	fprintf(fileres,"loglikelihood\n");
	fclose(fileres);
	
	//save posterior for initial samples
	k=0;
	for(j=0;j<nbeta;j++){posterior[index2(kcoda,j,nsavecoda)]=beta[j];} k=nbeta;
	for(j=0;j<ntheta;j++){posterior[index2(kcoda,j+k,nsavecoda)]=theta[j];} k+=ntheta;
	for(j=0;j<nbetagroup;j++){posterior[index2(kcoda,j+k,nsavecoda)]=betastatus[j];} k+=nbetagroup;
	for(j=0;j<npsi;j++){posterior[index2(kcoda,j+k,nsavecoda)]=psi[j];} k+=npsi;
	if(fixed==0)
	{
		for(j=0;j<nbeta;j++){posterior[index2(kcoda,j+k,nsavecoda)]=sdb[j];} k+=nbeta;
	}
	posterior[index2(kcoda,k,nsavecoda)]=varp;
	posterior[index2(kcoda,k+1,nsavecoda)]=loglikeorig;
	kcoda++;
	
  	for(i=1;i<niter;i++)
 	{
/* 		printf("i=%d\n",i);*/
 		//update linear terms in random order
 		gsl_ran_shuffle(rand_gen,betashuffle,nvariables,sizeof(int));
		if(varselect==1)
		{
			for(j=0;j<nvariables;j++)
			{
				//select variable to update from shuffled set
				k=betashuffle[j];
				//if variable is currently excluded, then propose to keep or add
				if(betastatusvar[k]==0)
				{
					//check whether variable selection can occur and if so check and amend probabilities for interaction effects
					if(interaction[k]>0)
					{
						//check that earlier main and interaction effects are included if we wish to add this variable
						pinteraction=pvecvar[0];
						for(m=0;m<nvariables;m++)
						{
							if(intfactor[index2(m,k,nbetagroup)]==1 && betastatusvar[m]==2) pinteraction=1.0;
						}
						pvecvarmove[0]=pinteraction;
						//if variable is currently excluded, then propose to keep or add
						if(gsl_ran_flat(rand_gen,0.0,1.0)>pvecvar[0])
						{
							betastatusvar[k]=1;
							//add variable to model choosing PO or NPO status for each parameter randomly
							for(m=xassign[k];m<xassign[k+1];m++)
							{
								if(gsl_ran_flat(rand_gen,0.0,1.0)<pzero[0]) betastatus[m]=0;
								else betastatus[m]=1;
							}
							moveexctoinc(k,xassign,betastatusvar,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,pvecvarmove,pzero,maxsdb,fixed);
						}
						//now check validity conditions
						if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
					}
					else
					{
						if(gsl_ran_flat(rand_gen,0.0,1.0)>pvecvar[0])
						{
							betastatusvar[k]=1;
							//add variable to model choosing PO or NPO status for each parameter randomly
							for(m=xassign[k];m<xassign[k+1];m++)
							{
								if(gsl_ran_flat(rand_gen,0.0,1.0)<pzero[0]) betastatus[m]=0;
								else betastatus[m]=1;
							}
							moveexctoinc(k,xassign,betastatusvar,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,pvecvar,pzero,maxsdb,fixed);
						}
						//now check validity conditions
						if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
					}
				}
				else
				{
					//check whether variable selection can occur and if so check and amend probabilities for interaction effects
					if(maxinteraction>0)
					{
						//check that later interaction effects are not currently included if we wish to drop this variable
						pinteraction=pvecvar[0];
						if(maxinteraction>interaction[k])
						{
							for(m=intstart[interaction[k]+1];m<intstart[maxinteraction+1];m++)
							{
								if(intfactor[index2(k,m,nvariables)]==1 && betastatusvar[m]!=2) pinteraction=1.0;
							}
						}
						pvecvarmove[0]=pinteraction;
						//if variable is currently included, then propose to keep or drop
						if(gsl_ran_flat(rand_gen,0.0,1.0)>pvecvarmove[0])
						{
							betastatusvar[k]=0;
							//drop variable from model
							moveinctoexc(k,xassign,betastatusvar,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,pvecvarmove,pzero,maxsdb,fixed);
						}
						//now check validity conditions
						if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
					}
					else
					{
						//if variable is currently included, then propose to keep or drop
						if(gsl_ran_flat(rand_gen,0.0,1.0)>pvecvar[0])
						{
							betastatusvar[k]=0;
							//drop variable from model
							moveinctoexc(k,xassign,betastatusvar,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,pvecvar,pzero,maxsdb,fixed);
						}
						//now check validity conditions
						if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
					}
				}
			}
		}
		//update linear terms in random order
 		gsl_ran_shuffle(rand_gen,betashufflepar,nbetagroup,sizeof(int));
 		for(j=0;j<nbetagroup;j++)
 		{
 			k=betashufflepar[j];
			//if variable currently included in model then update all parameters associated with it
			if(betastatus[k]!=2)
			{
				//if parameter is currently PO, then choose to move parameter or shift to NPO
				if(betastatus[k]==0)
				{
					u=gsl_ran_flat(rand_gen,0.0,1.0);
					if(u<pvec[0]) movePO(k,n,nbetagroup,ntheta,beta,theta,psi,variables,propsdb,&loglikeorig,rand_gen,mnb,sdb);
					else movePOtoNPO(k,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,maxsdb,fixed);
				}
				else
				{
					//if parameter is currently NPO, then choose to move parameter or shift to PO
					u=gsl_ran_flat(rand_gen,0.0,1.0);
					if(u<pvec[0]) moveNPO(k,n,nbetagroup,ntheta,beta,theta,psi,variables,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt);
					else moveNPOtoPO(k,n,nbetagroup,ntheta,beta,theta,psi,variables,betastatus,betacond,minvar,maxvar,propsdb,&loglikeorig,rand_gen,mnb,sdb,sdt,maxsdb,fixed);
				}
				//now check validity conditions
				if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
			}
		}			
		//update theta terms in random order
 		gsl_ran_shuffle(rand_gen,betashuffletheta,ntheta,sizeof(int));
		for(m=0;m<ntheta;m++) movetheta(betashuffletheta[m],n,nbetagroup,ntheta,beta,theta,psi,variables,betacond,propsdt,&loglikeorig,rand_gen,sdt);
		//now check validity conditions
		if(validitycheck(ntheta,nbetagroup,betastatus,beta,theta,betacond,minvar,maxvar)==1) return 1;
	
		//move SD for beta parameter
		if(fixed==0)
		{
			//update SD terms in random order
 			gsl_ran_shuffle(rand_gen,betashufflepar,nbetagroup,sizeof(int));
			for(j=0;j<nbetagroup;j++)
			{
				k=betashufflepar[j];
				if(betastatus[k]==0) movesdb(k,0,nbetagroup,beta,propsdb,rand_gen,mnb,sdb,maxsdb);
				else
				{
					if(betastatus[k]==1)
					{
						gsl_ran_shuffle(rand_gen,betashuffletheta,ntheta,sizeof(int));
						for(m=0;m<ntheta;m++) movesdb(k,betashuffletheta[m],nbetagroup,beta,propsdb,rand_gen,mnb,sdb,maxsdb);
					}
				}
			}
		}
		if(RE==1)
		{
			//randomly update individual effect terms
			sdp=sqrt(varp);
			for(j=0;j<((int) floor(0.3*npsi));j++)
			{
				k=(int) floor(npsi*gsl_ran_flat(rand_gen,0.0,1.0));
				psistore=psi[k];
				loglikeorig=loglikelihood_psi(beta,theta,psi,variables,n,ntheta,nbetagroup,psicount,k);
				psi[k]=psistore+gsl_ran_flat(rand_gen,-propsdp,propsdp);
				loglikeprop=loglikelihood_psi(beta,theta,psi,variables,n,ntheta,nbetagroup,psicount,k);
				//work out acceptance probability  				
				if(isfinite(loglikeprop)!=0)
				{
					accorig=loglikeorig;
					accprop=loglikeprop;
					//adjust for priors
					accorig+=log(gsl_ran_gaussian_pdf((psistore-mnpsi),sdp));
					accprop+=log(gsl_ran_gaussian_pdf((psi[k]-mnpsi),sdp));
					//proposals cancel
					acc=accprop-accorig;
					acc=exp(acc);
					if(gsl_rng_uniform_pos(rand_gen)>=acc) psi[k]=psistore;
				}
				else psi[k]=psistore;
			}
		
			//update full likelihood
			loglikeorig=loglikelihood(beta,theta,psi,variables,n,ntheta,nbetagroup);
	
			//update individual effect variance term
			varpstore=varp;
			varp+=gsl_ran_flat(rand_gen,-propsdvarp,propsdvarp);
			if(varp>0)
			{
				accorig=0.0; accprop=0.0;
				//adjust for priors
				sdp=sqrt(varpstore);
				for(j=0;j<npsi;j++) accorig+=log(gsl_ran_gaussian_pdf((psi[j]-mnpsi),sdp));
				sdp=sqrt(varp);
				for(j=0;j<npsi;j++) accprop+=log(gsl_ran_gaussian_pdf((psi[j]-mnpsi),sdp));
				accorig+=log(gsl_ran_gamma_pdf(varpstore,shvarp,1.0/rtvarp));
				accprop+=log(gsl_ran_gamma_pdf(varp,shvarp,1.0/rtvarp));
				//proposals cancel
				acc=accprop-accorig;
				acc=exp(acc);
				if(gsl_rng_uniform_pos(rand_gen)>=acc) varp=varpstore;
			}
			else varp=varpstore;
		}
/*		//record posterior*/
/*		temp=loglikeorig;*/
/*		//add beta and SD priors*/
/*		for(k=0;k<nbetagroup;k++)*/
/*		{*/
/*			if(betastatus[k]<2)*/
/*			{*/
/*				if(betastatus[k]==0) temp+=log(gsl_ran_gaussian_pdf((beta[index2(k,0,nbetagroup)]-mnb),sdb[index2(k,0,nbetagroup)]));*/
/*				else*/
/*				{*/
/*					for(j=0;j<ntheta;j++) temp+=log(gsl_ran_gaussian_pdf((beta[index2(k,j,nbetagroup)]-mnb),sdb[index2(k,j,nbetagroup)]));*/
/*				}*/
/*				if(fixed==0) temp-=log(maxsdb);*/
/*			}*/
/*		}*/
/*		//add theta priors*/
/*		temp+=log(gsl_ran_gaussian_pdf(theta[0],sdt));*/
/*		for(k=1;k<ntheta;k++) temp+=log(gsl_ran_gaussian_pdf(theta[k],sdt))-log(1.0-gsl_cdf_gaussian_P(theta[k-1]-betacond[k-1],sdt));*/
/*		//add RE priors*/
/*		if(RE==1)*/
/*		{*/
/*			for(k=0;k<npsi;k++) temp+=log(gsl_ran_gaussian_pdf((psi[j]-mnpsi),sqrt(varp)));*/
/*			temp+=log(gsl_ran_gamma_pdf(varp,shvarp,1.0/rtvarp));*/
/*		}	*/
	
		//record output		
		k=0;
		for(j=0;j<nbeta;j++){posterior[index2(kcoda,j,nsavecoda)]=beta[j];} k=nbeta;
		for(j=0;j<ntheta;j++){posterior[index2(kcoda,j+k,nsavecoda)]=theta[j];} k+=ntheta;
		for(j=0;j<nbetagroup;j++){posterior[index2(kcoda,j+k,nsavecoda)]=betastatus[j];} k+=nbetagroup;
		for(j=0;j<npsi;j++){posterior[index2(kcoda,j+k,nsavecoda)]=psi[j];} k+=npsi;
		if(fixed==0)
		{
			for(j=0;j<nbeta;j++){posterior[index2(kcoda,j+k,nsavecoda)]=sdb[j];} k+=nbeta;
		}
		posterior[index2(kcoda,k,nsavecoda)]=varp;
		posterior[index2(kcoda,k+1,nsavecoda)]=loglikeorig;
		m=k+1;
		kcoda++;
	
		//save output into text files for examination/editing in R
		if((i+1)%nsavecoda==0)
		{
			time (&end);
	  		runtime=difftime(end,start);
			printf("Chain %s: current iterations=%d\tt=%f\n",argv[1],i+1,runtime);
			//append current results to output		
			fileres=fopen(filename,"a");
			for(k=0;k<kcoda;k++)
			{
				for(j=0;j<(m+1);j++) fprintf(fileres,"%f\t",posterior[index2(k,j,nsavecoda)]);
				fprintf(fileres,"\n");
			}
			kcoda=0;
			fclose(fileres);
 		}
	}		
  	// terminate the program:
  	return 0;
}
