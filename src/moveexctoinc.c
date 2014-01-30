#include <R.h>
#include <Rmath.h>

#include"functions.h"

//function to add parameters for variable not currently present in model
void moveexctoinc(int k, int *xassign, int *betastatusvar, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double *pvec, double *pzero, double maxsdb, int fixed)
{
    int m, l, valid;
    double temp;
    double loglikeprop;
    double accorig, accprop, acc;
    double *betacondstore = (double *) Calloc((ntheta - 1), double);
    for(m = 0; m < (ntheta - 1); m++) betacondstore[m] = betacond[m];

    //set parameters
    valid = 1;
    for(l = xassign[k]; l < xassign[k + 1]; l++)
    {
        if(betastatus[l] == 0)
        {
            //if parameter is to have PO structure
            beta[index2(l, 0, nbetagroup)] = rnorm(0.0, propsdb);
            for(m = 1; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = beta[index2(l, 0, nbetagroup)];
        }
        else
        {
            //if parameter is to have NPO structure
            for(m = 0; m < ntheta; m++) beta[index2(l, m, nbetagroup)] = rnorm(0.0, propsdb);

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
                    sdb[index2(l, 0, nbetagroup)] = runif(0.0, maxsdb);
                    for(m = 1; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = sdb[index2(l, 0, nbetagroup)];
                }
                else
                {
                    //sample new parameters for NPO
                    for(m = 0; m < ntheta; m++) sdb[index2(l, m, nbetagroup)] = runif(0.0, maxsdb);
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
                if(betastatus[l] == 0) accprop += dnorm((beta[index2(l, 0, nbetagroup)] - mnb), 0.0, sdb[index2(l, 0, nbetagroup)], 1);
                else
                {
                    for(m = 0; m < ntheta; m++) accprop += dnorm((beta[index2(l, m, nbetagroup)] - mnb), 0.0, sdb[index2(l, m, nbetagroup)], 1);
                    //adjust for theta priors conditional on beta
                    for(m = 1; m < ntheta; m++)
                    {
                        accorig -= pnorm(theta[m - 1] - betacondstore[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore[m - 1], sdt));
                        accprop -= pnorm(theta[m - 1] - betacond[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
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
                if(betastatus[l] == 0) accorig += dnorm(beta[index2(l, 0, nbetagroup)], 0.0, propsdb, 1);
                else
                {
                    for(m = 0; m < ntheta; m++) accorig += dnorm(beta[index2(l, m, nbetagroup)], 0.0, propsdb, 1);
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
            if(runif(0.0, 1.0) < acc) (*loglikeorig) = loglikeprop;
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
    Free(betacondstore);
    return;
}
