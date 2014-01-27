#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include"functions.h"

//function to drop parameters for variable currently present in model
void moveinctoexc(int k, int *xassign, int *betastatusvar, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, gsl_rng *rand_gen, double mnb, double *sdb, double sdt, double *pvec, double *pzero, double maxsdb, int fixed)
{
    int m, l;
    double temp;
    double loglikeprop;
    double accorig, accprop, acc;

    int npars = xassign[k + 1] - xassign[k];
    //declare variables for storing original values
    int * betastatusstore = (int *)malloc(npars * sizeof(int));
    double *betastore = (double *)malloc(npars * ntheta * sizeof(double));
    double *sdstore = (double *)malloc(npars * ntheta * sizeof(double));
    double *betacondstore = (double *)malloc((ntheta - 1) * sizeof(double));

    //store conditions
    for(m = 0; m < (ntheta - 1); m++) betacondstore[m] = betacond[m];

    //drop parameters out of the model
    for(l = 0; l < npars; l++)
    {
        betastatusstore[l] = betastatus[xassign[k] + l];
        betastatus[xassign[k] + l] = 2;
        if(betastatusstore[l] == 1)
        {
            //update conditions
            for(m = 0; m < (ntheta - 1); m++)
            {
                temp = beta[index2(xassign[k] + l, m, nbetagroup)] - beta[index2(xassign[k] + l, m + 1, nbetagroup)];
                betacond[m] -= ((minvar[xassign[k] + l] * temp) < (maxvar[xassign[k] + l] * temp) ? (minvar[xassign[k] + l] * temp) : (maxvar[xassign[k] + l] * temp));
            }
        }
        //set parameters to zero for variable
        for(m = 0; m < ntheta; m++)
        {
            betastore[index2(m, l, ntheta)] = beta[index2(xassign[k] + l, m, nbetagroup)];
            beta[index2(xassign[k] + l, m, nbetagroup)] = 0.0;
        }
        if(fixed == 0)
        {
            for(m = 0; m < ntheta; m++)
            {
                sdstore[index2(m, l, ntheta)] = sdb[index2(xassign[k] + l, m, nbetagroup)];
                sdb[index2(xassign[k] + l, m, nbetagroup)] = 0.0;
            }
        }
    }

    //calculate log-likelihood
    loglikeprop = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
    //work out acceptance probability
    if(isfinite(loglikeprop) != 0)
    {
        accorig = (*loglikeorig);
        accprop = loglikeprop;
        //adjust for beta priors
        for(l = 0; l < npars; l++)
        {
            if(betastatusstore[l] == 0) accorig += log(gsl_ran_gaussian_pdf((betastore[index2(0, l, ntheta)] - mnb), sdstore[index2(0, l, ntheta)]));
            else
            {
                for(m = 0; m < ntheta; m++) accorig += log(gsl_ran_gaussian_pdf((betastore[index2(m, l, ntheta)] - mnb), sdstore[index2(m, l, ntheta)]));
                //adjust for theta priors conditional on beta
                for(m = 1; m < ntheta; m++)
                {
                    accorig -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore[m - 1], sdt));
                    accprop -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
                }
            }
        }

        //adjust for SD prior
        if(fixed == 0)
        {
            for(l = 0; l < npars; l++)
            {
                if(betastatusstore[l] == 0) accorig -= log(maxsdb);
                else accorig -= ntheta * log(maxsdb);
            }
        }

        //adjust for proposals
        for(l = 0; l < npars; l++)
        {
            if(betastatusstore[l] == 0) accprop += log(gsl_ran_gaussian_pdf(betastore[index2(0, l, ntheta)], propsdb));
            else
            {
                for(m = 0; m < ntheta; m++) accprop += log(gsl_ran_gaussian_pdf(betastore[index2(m, l, ntheta)], propsdb));
            }
            if(fixed == 0)
            {
                if(betastatusstore[l] == 0) accprop -= log(maxsdb);
                else accprop -= ntheta * log(maxsdb);
            }
        }

        //adjust for move probabilities

        //prob. of dropping variable
        accorig += log(pvec[1] - pvec[0]);
        //prob. of adding variable in reverse move
        accprop += log(pvec[1] - pvec[0]);
        for(l = 0; l < npars; l++)
        {
            //prob. of zero -> PO/NPO
            if(betastatusstore[l] == 0) accprop += log(pzero[0]);
            else accprop += log(pzero[1] - pzero[0]);
        }

        //calculate acceptance
        acc = accprop - accorig;
        acc = exp(acc);
        if(gsl_rng_uniform_pos(rand_gen) < acc) (*loglikeorig) = loglikeprop;
        else
        {
            //if invalid proposal then reset everything to previous values
            for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
            for(l = 0; l < npars; l++)
            {
                betastatus[xassign[k] + l] = betastatusstore[l];
                for(m = 0; m < ntheta; m++) beta[index2(xassign[k] + l, m, nbetagroup)] = betastore[index2(m, l, ntheta)];
                if(fixed == 0)
                {
                    for(m = 0; m < ntheta; m++) sdb[index2(xassign[k] + l, m, nbetagroup)] = sdstore[index2(m, l, ntheta)];
                }
            }
            betastatusvar[k] = 1;
        }
    }
    else
    {
        //if invalid proposal then reset everything to previous values
        for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
        for(l = 0; l < npars; l++)
        {
            betastatus[xassign[k] + l] = betastatusstore[l];
            for(m = 0; m < ntheta; m++) beta[index2(xassign[k] + l, m, nbetagroup)] = betastore[index2(m, l, ntheta)];
            if(fixed == 0)
            {
                for(m = 0; m < ntheta; m++) sdb[index2(xassign[k] + l, m, nbetagroup)] = sdstore[index2(m, l, ntheta)];
            }
        }
        betastatusvar[k] = 1;
    }
    //free memory from the heap
    free(betastatusstore);
    free(betastore);
    free(betacondstore);
    free(sdstore);
    return;
}
