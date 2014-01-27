#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include"functions.h"

//function to add parameters for variable not currently present in model
void moveexctoinc(int k, int *xassign, int *betastatusvar, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, gsl_rng *rand_gen, double mnb, double *sdb, double sdt, double *pvec, double *pzero, double maxsdb, int fixed)
{
    int m, l, valid;
    double temp;
    double loglikeprop;
    double accorig, accprop, acc;
    double *betacondstore = (double *)malloc((ntheta - 1) * sizeof(double));
    for(m = 0; m < (ntheta - 1); m++) betacondstore[m] = betacond[m];

    //set parameters
    valid = 1;
    for(l = xassign[k]; l < xassign[k + 1]; l++)
    {
        if(betastatus[l] == 0)
        {
            //if parameter is to have PO structure
            beta[index2(l, 0, nbetagroup)] = gsl_ran_gaussian(rand_gen, propsdb);
            for(m = 1; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = beta[index2(l, 0, nbetagroup)];
        }
        else
        {
            //if parameter is to have NPO structure
            for(m = 0; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = gsl_ran_gaussian(rand_gen, propsdb);

            //update conditions
            for(m = 0; m < (ntheta - 1); m++)
            {
                temp = beta[index2(l, m, nbetagroup)] - beta[index2(l, m + 1, nbetagroup)];
                betacond[m] += ((minvar[l] * temp) < (maxvar[l] * temp) ? (minvar[l] * temp) : (maxvar[l] * temp));
            }

            //check validity
            m = 0;
            while(m < (ntheta - 1) && valid == 1)
            {
                if((theta[m] - theta[m + 1]) > betacond[m]) valid = 0;
                m++;
            }
        }
    }

    if(valid == 1)
    {
        //if we don't have a fixed prior variance for the regression parameters
        //then update SDs
        if(fixed == 0)
        {
            for(l = xassign[k]; l < xassign[k + 1]; l++)
            {
                if(betastatus[l] == 0)
                {
                    //sample new parameters for PO
                    sdb[index2(l, 0, nbetagroup)] = gsl_ran_flat(rand_gen, 0.0, maxsdb);
                    for(m = 1; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = sdb[index2(l, 0, nbetagroup)];
                }
                else
                {
                    //sample new parameters for NPO
                    for(m = 0; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = gsl_ran_flat(rand_gen, 0.0, maxsdb);
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
            for(l = xassign[k]; l < xassign[k + 1]; l++)
            {
                if(betastatus[l] == 0) accprop += log(gsl_ran_gaussian_pdf((beta[index2(l, 0, nbetagroup)] - mnb), sdb[index2(l, 0, nbetagroup)]));
                else
                {
                    for(m = 0; m < ntheta; m++) accprop += log(gsl_ran_gaussian_pdf((beta[index2(l, m, nbetagroup)] - mnb), sdb[index2(l, m, nbetagroup)]));
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
                for(l = xassign[k]; l < xassign[k + 1]; l++)
                {
                    if(betastatus[l] == 0) accprop -= log(maxsdb);
                    else accprop -= ntheta * log(maxsdb);
                }
            }

            //adjust for proposals
            for(l = xassign[k]; l < xassign[k + 1]; l++)
            {
                if(betastatus[l] == 0) accorig += log(gsl_ran_gaussian_pdf(beta[index2(l, 0, nbetagroup)], propsdb));
                else
                {
                    for(m = 0; m < ntheta; m++) accorig += log(gsl_ran_gaussian_pdf(beta[index2(l, m, nbetagroup)], propsdb));
                }
                if(fixed == 0)
                {
                    if(betastatus[l] == 0) accorig -= log(maxsdb);
                    else accorig -= ntheta * log(maxsdb);
                }
            }

            //adjust for move probabilities

            //prob. of adding variable
            accorig += log(pvec[1] - pvec[0]);
            for(l = xassign[k]; l < xassign[k + 1]; l++)
            {
                //prob. of zero -> PO/NPO
                if(betastatus[l] == 0) accorig += log(pzero[0]);
                else accorig += log(pzero[1] - pzero[0]);
            }
            //prob. of dropping variable in reverse move
            accprop += log(pvec[1] - pvec[0]);

            //calculate acceptance
            acc = accprop - accorig;
            acc = exp(acc);
            if(gsl_rng_uniform_pos(rand_gen) < acc) (*loglikeorig) = loglikeprop;
            else
            {
                //if invalid proposal then reset everything to zero
                for(l = xassign[k]; l < xassign[k + 1]; l++)
                {
                    for(m = 0; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = 0.0;
                    if(fixed == 0) {
                        for(m = 0; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = 0.0;
                    }
                    betastatus[l] = 2;
                }
                for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
                betastatusvar[k] = 0;
            }
        }
        else
        {
            //if invalid proposal then reset everything to zero
            for(l = xassign[k]; l < xassign[k + 1]; l++)
            {
                for(m = 0; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = 0.0;
                if(fixed == 0) {
                    for(m = 0; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = 0.0;
                }
                betastatus[l] = 2;
            }
            for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
            betastatusvar[k] = 0;
        }
    }
    else
    {
        //if invalid proposal then reset everything to zero
        for(l = xassign[k]; l < xassign[k + 1]; l++)
        {
            for(m = 0; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = 0.0;
            /*			if(fixed==0){for(m=0;m<ntheta;m++) sdb[index2(l,m,nbetagroup)]=0.0;}*/
            betastatus[l] = 2;
        }
        for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
        betastatusvar[k] = 0;
    }
    //free memory from the heap
    free(betacondstore);
    return;
}
