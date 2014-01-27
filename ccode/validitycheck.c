#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"

//function to check validity of proposals
int validitycheck(int ntheta, int nbetagroup, int *betastatus, double *beta, double *theta, double *betacond, double *minvar, double *maxvar)
{
	int j,m;
	double betastore,temp;
	
	//now check validity conditions
	for(j=0;j<(ntheta-1);j++)
	{
		betastore=0.0;
		for(m=0;m<nbetagroup;m++)
		{
			if(betastatus[m]==1)
			{
				temp=beta[index2(m,j,nbetagroup)]-beta[index2(m,j+1,nbetagroup)];
				betastore+=(minvar[m]*temp<maxvar[m]*temp ? (minvar[m]*temp):(maxvar[m]*temp));
			}
		}
		if(fabs(betastore-betacond[j])>1e-10)
		{
			printf("Something invalid in proposals for beta1\n");
			return 1;
		}
		if((theta[j]-theta[j+1])>betastore)
		{
			printf("Something invalid in proposals for beta2\n");
			return 1;
		}
	}
	return 0;
}
