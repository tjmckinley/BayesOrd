#include <R.h>
#include <Rmath.h>
#include"functions.h"

//function to move PO parameter
void movePO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double propsdb, double *loglikeorig, double mnb, double *sdb)
{
    int m;
    double betastore, loglikeprop;
    double accorig, accprop, acc;

    betastore = beta[index2(k, 0, nbetagroup)];
    beta[index2(k, 0, nbetagroup)] += rnorm(0.0, propsdb);
    for(m = 1; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = beta[index2(k, 0, nbetagroup)];

    //calculate log-likelihood
    loglikeprop = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
    //work out acceptance probability
    if(isfinite(loglikeprop) != 0)
    {
        accorig = (*loglikeorig);
        accprop = loglikeprop;
        //adjust for beta priors
        accorig += dnorm((betastore - mnb), 0.0, sdb[index2(k, 0, nbetagroup)], 1);
        accprop += dnorm((beta[index2(k, 0, nbetagroup)] - mnb), 0.0, sdb[index2(k, 0, nbetagroup)], 1);
        //proposals cancel
        acc = accprop - accorig;
        acc = exp(acc);
        if(runif(0.0, 1.0) < acc) (*loglikeorig) = loglikeprop;
        else
        {
            for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore;
        }
    }
    else
    {
        for(m = 0; m < ntheta; m++) beta[index2(k, m, nbetagroup)] = betastore;
    }
    return;
}
