#include <R.h>
#include <Rmath.h>
#include"functions.h"

//function to move PO parameter to NPO parameter
void movePOtoNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double maxsdb, int fixed)
{
    int m, valid;
    double *betacondstore = (double *) Calloc(ntheta - 1, double);
    double *u = (double *) Calloc(ntheta - 1, double);
    double betastore, loglikeprop, temp, sdbstore, lower = 0.0, upper = 0.0;
    double accorig, accprop, acc;

    //propose move from PO to NPO model
    for(m = 0; m < (ntheta - 1); m++) betacondstore[m] = betacond[m];

    //first level is current value
    betastore = beta[index2(k, 0, nbetagroup)];

    //update beta parameters
    beta[index2(k, 0, nbetagroup)] = 0.0;
    for(m = 0; m < (ntheta - 1); m++)
    {
        u[m] = runif(-propsdb, propsdb);
        beta[index2(k, 0, nbetagroup)] += u[m];
    }
    beta[index2(k, 0, nbetagroup)] += betastore * (4.0 - (ntheta + 1));
    beta[index2(k, 0, nbetagroup)] = (ntheta / 2.0) * beta[index2(k, 0, nbetagroup)];
    for(m = 1; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = (ntheta / 2.0) * (betastore - u[m - 1]);

    //update SD parameters
    sdbstore = sdb[index2(k, 0, nbetagroup)];
    if(fixed == 0)
    {
        lower = ((sdbstore - propsdb) < 0.0 ? 0.0 : (sdbstore - propsdb));
        upper = ((sdbstore + propsdb) > maxsdb ? maxsdb : (sdbstore + propsdb));
        for(m = 1; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = runif(lower, upper);
    }

    //update validity conditions
    m = 0;
    valid = 1;
    while(m < (ntheta - 1) && valid == 1)
    {
        temp = beta[index2(k, m, nbetagroup)] - beta[index2(k, m + 1, nbetagroup)];
        betacond[m] += ((minvar[k] * temp) < (maxvar[k] * temp) ? (minvar[k] * temp) : (maxvar[k] * temp));
        if((theta[m] - theta[m + 1]) > betacond[m]) valid = 0;
        m++;
    }
    if(valid == 1)
    {
        //calculate log-likelihood
        loglikeprop = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
        //work out acceptance probability
        if(isfinite(loglikeprop) != 0)
        {
            accorig = (*loglikeorig);
            accprop = loglikeprop;

            //adjust for beta priors
            accorig += dnorm((betastore - mnb), 0.0, sdbstore, 1);
            for(m = 0; m < ntheta; m++) accprop += dnorm((beta[index2(k, m, nbetagroup)] - mnb), 0.0, sdb[index2(k, m, nbetagroup)], 1);

            if(fixed == 0)
            {
                //adjust for SD priors
                accorig -= log(maxsdb);
                accprop -= ntheta * log(maxsdb);
            }

            //adjust for theta priors conditional on beta
            for(m = 1; m < ntheta; m++)
            {
                accorig -= pnorm(theta[m - 1] - betacondstore[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore[m - 1], sdt));
                accprop -= pnorm(theta[m - 1] - betacond[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
            }

            //adjust for proposals
            accorig -= (ntheta - 1) * log(2 * propsdb);
            if(fixed == 0) accorig -= (ntheta - 1) * log(upper - lower);

            //equal probs of move proposal cancel

            //adjust for Jacobian
            accprop += log(ntheta) + (ntheta - 1.0) * log(ntheta / 2.0);

            //calculate acceptance
            acc = accprop - accorig;
            acc = exp(acc);
            if(runif(0.0, 1.0) < acc)
            {
                (*loglikeorig) = loglikeprop;
                //set beta status
                betastatus[k] = 1;
            }
            else
            {
                for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore;
                if(fixed == 0) {
                    for(m = 0; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore;
                }
                for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
            }
        }
        else
        {
            for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore;
            if(fixed == 0) {
                for(m = 0; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore;
            }
            for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
        }
    }
    else
    {
        for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore;
        if(fixed == 0) {
            for(m = 0; m < ntheta; m++) sdb[index2(k, m, nbetagroup)] = sdbstore;
        }
        for(m = 0; m < (ntheta - 1); m++) betacond[m] = betacondstore[m];
    }
    //free memory from the heap
    Free(betacondstore);
    Free(u);
    return;
}
