#include <R.h>
#include <Rmath.h>
#include"functions.h"

//function to move theta parameter
void movetheta(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double *betacond, double propsdt, double *loglikeorig, double sdt)
{
    double thetastore, loglikeprop;
    double accorig, accprop, acc, low, upp, low1, upp1;

    thetastore = theta[k];
    //set bounds for forward move
    if(k == 0)
    {
        low = thetastore - propsdt;
        upp = ((thetastore + propsdt) > (theta[k + 1] + betacond[k]) ? (theta[k + 1] + betacond[k]) : (thetastore + propsdt));
    }
    else
    {
        if(k < (ntheta - 1))
        {
            low = ((thetastore - propsdt) < (theta[k - 1] - betacond[k - 1]) ? (theta[k - 1] - betacond[k - 1]) : (thetastore - propsdt));
            upp = ((thetastore + propsdt) > (theta[k + 1] + betacond[k]) ? (theta[k + 1] + betacond[k]) : (thetastore + propsdt));
        }
        else
        {
            low = ((thetastore - propsdt) < (theta[k - 1] - betacond[k - 1]) ? (theta[k - 1] - betacond[k - 1]) : (thetastore - propsdt));
            upp = thetastore + propsdt;
        }
    }
    //sample new parameter value
    theta[k] = runif(low, upp);
    //set bounds for reverse move
    if(k == 0)
    {
        low1 = theta[k] - propsdt;
        upp1 = ((theta[k] + propsdt) > (theta[k + 1] + betacond[k]) ? (theta[k + 1] + betacond[k]) : (theta[k] + propsdt));
    }
    else
    {
        if(k < (ntheta - 1))
        {
            low1 = ((theta[k] - propsdt) < (theta[k - 1] - betacond[k - 1]) ? (theta[k - 1] - betacond[k - 1]) : (theta[k] - propsdt));
            upp1 = ((theta[k] + propsdt) > (theta[k + 1] + betacond[k]) ? (theta[k + 1] + betacond[k]) : (theta[k] + propsdt));
        }
        else
        {
            low1 = ((theta[k] - propsdt) < (theta[k - 1] - betacond[k - 1]) ? (theta[k - 1] - betacond[k - 1]) : (theta[k] - propsdt));
            upp1 = theta[k] + propsdt;
        }
    }

    //calculate likelihood
    loglikeprop = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
    //work out acceptance probability
    if(isfinite(loglikeprop) != 0)
    {
        accorig = (*loglikeorig);
        accprop = loglikeprop;
        //adjust for priors
        if(k == 0)
        {
            accorig += dnorm(thetastore, 0.0, sdt, 1);
            accprop += dnorm(theta[k], 0.0, sdt, 1);
            accorig -= pnorm(thetastore - betacond[k], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(thetastore - betacond[k], sdt));
            accprop -= pnorm(theta[k] - betacond[k], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[k] - betacond[k], sdt));
        }
        else
        {
            if(k < (ntheta - 1))
            {
                accorig += dnorm(thetastore, 0.0, sdt, 1);
                accprop += dnorm(theta[k], 0.0, sdt, 1);
                accorig -= pnorm(thetastore - betacond[k], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(thetastore - betacond[k], sdt));
                accprop -= pnorm(theta[k] - betacond[k], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[k] - betacond[k], sdt));
            }
            else
            {
                accorig += dnorm(thetastore, 0.0, sdt, 1);
                accprop += dnorm(theta[k], 0.0, sdt, 1);
            }
        }
        //adjust for proposals
        accorig -= log(upp - low);
        accprop -= log(upp1 - low1);
        //calculate acceptance
        acc = accprop - accorig;
        acc = exp(acc);
        if(runif(0.0, 1.0) < acc) (*loglikeorig) = loglikeprop;
        else theta[k] = thetastore;
    }
    else theta[k] = thetastore;

    return;
}
