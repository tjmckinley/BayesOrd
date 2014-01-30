#include <R.h>
#include<Rmath.h>
#include"functions.h"

//function to move NPO parameter to PO parameter
void moveNPOtoPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double maxsdb, int fixed)
{
    int m;
    double *betacondstore = (double *) Calloc(ntheta - 1, double);
    double *betastore = (double *) Calloc(ntheta, double);
    double *sdbstore = (double *) Calloc(ntheta, double);
    double loglikeprop, temp, lower = 0.0, upper = 0.0;
    double accorig, accprop, acc;

    //propose move from NPO to PO model

    //update conditions
    for(m = 0; m < (ntheta - 1); m++)
    {
        betacondstore[m] = betacond[m];
        temp = beta[index2(k, m, nbetagroup)] - beta[index2(k, m + 1, nbetagroup)];
        betacond[m] -= ((minvar[k] * temp) < (maxvar[k] * temp) ? (minvar[k] * temp) : (maxvar[k] * temp));
    }

    for(m = 0; m < ntheta; m++) betastore[m] = beta[index2(k, m, nbetagroup)];
    if(fixed == 0) for(m = 0; m < ntheta; m++) sdbstore[m] = sdb[index2(k, m, nbetagroup)];

    //calculate ucand variables to calculate reverse move
    temp = 0.0;
    for(m = 0; m < ntheta; m++) temp += beta[index2(k, m, nbetagroup)];
    temp = temp / ((double) ntheta);
    //set new values for beta
    for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = temp;

    if(fixed == 0)
    {
        lower = ((sdbstore[0] - propsdb) < 0.0 ? 0.0 : (sdbstore[0] - propsdb));
        upper = ((sdbstore[0] + propsdb) > maxsdb ? maxsdb : (sdbstore[0] + propsdb));
        for(m = 1; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore[0];
    }

    //calculate log-likelihood
    loglikeprop = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
    //work out acceptance probability
    if(isfinite(loglikeprop) != 0)
    {
        accorig = (*loglikeorig);
        accprop = loglikeprop;
        //adjust for beta priors
        accprop += dnorm((beta[index2(k, 0, nbetagroup)] - mnb), 0.0, sdb[index2(k, 0, nbetagroup)], 1);
        for(m = 0; m < ntheta; m++) accorig += dnorm((betastore[m] - mnb), 0.0, sdbstore[m], 1);
        //adjust for theta priors conditional on beta
        for(m = 1; m < ntheta; m++)
        {
            accorig -= pnorm(theta[m - 1] - betacondstore[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore[m - 1], sdt));
            accprop -= pnorm(theta[m - 1] - betacond[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
        }
        //adjust for SD priors
        if(fixed == 0)
        {
            accorig -= ntheta * log(maxsdb);
            accprop -= log(maxsdb);
        }

        //adjust for proposals for ucand
        accprop -= (ntheta - 1) * log(2 * propsdb);
        if(fixed == 0) accprop -= (ntheta - 1) * log(upper - lower);

        //equal probs of move proposal cancel

        //adjust for Jacobian
        accprop += (-log(ntheta) + (ntheta - 1.0) * log(2.0 / ((double) ntheta)));

        //calculate acceptance
        acc = accprop - accorig;
        acc = exp(acc);
        if(runif(0.0, 1.0) < acc)
        {
            (*loglikeorig) = loglikeprop;
            //set beta status
            betastatus[k] = 0;
        }
        else
        {
            for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore[m];
            if(fixed == 0) {
                for(m = 0; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore[m];
            }
            for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
        }
    }
    else
    {
        for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore[m];
        if(fixed == 0) {
            for(m = 0; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore[m];
        }
        for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
    }
    //free memory from the heap
    Free(betacondstore);
    Free(betastore);
    Free(sdbstore);
    return;
}
