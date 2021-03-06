#include <R.h>
#include <Rmath.h>
#include "functions.h"

//function to move NPO parameters
void moveNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt)
{
    int m, valid;
    double betastore, betacondstore = 0.0, betacondstore1 = 0.0, loglikeprop, temp;
    double accorig, accprop, acc;

    //update each regression parameter in turn for the given covariate
    for(m = 0; m < ntheta; m++)
    {
        betastore = beta[index2(k, m, nbetagroup)];
        //sample new parameter value
        beta[index2(k, m, nbetagroup)] = betastore + runif(-propsdb, propsdb);

        //update and check prior constraints
        if(m == 0)
        {
            betacondstore1 = betacond[m];
            //upper condition
            temp = betastore - beta[index2(k, m + 1, nbetagroup)];
            betacond[m] -= (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
            temp = beta[index2(k, m, nbetagroup)] - beta[index2(k, m + 1, nbetagroup)];
            betacond[m] += (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
            //check prior conditions
            valid = (theta[m] < (theta[m + 1] + betacond[m]) ? 1 : 0);
        }
        else
        {
            if(m == (ntheta - 1))
            {
                betacondstore = betacond[m - 1];
                //lower condition
                temp = beta[index2(k, m - 1, nbetagroup)] - betastore;
                betacond[m - 1] -= (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                temp = beta[index2(k, m - 1, nbetagroup)] - beta[index2(k, m, nbetagroup)];
                betacond[m - 1] += (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                //check prior conditions
                valid = (theta[m - 1] < (theta[m] + betacond[m - 1]) ? 1 : 0);
            }
            else
            {
                betacondstore = betacond[m - 1];
                betacondstore1 = betacond[m];
                //upper condition
                temp = betastore - beta[index2(k, m + 1, nbetagroup)];
                betacond[m] -= (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                temp = beta[index2(k, m, nbetagroup)] - beta[index2(k, m + 1, nbetagroup)];
                betacond[m] += (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                //lower condition
                temp = beta[index2(k, m - 1, nbetagroup)] - betastore;
                betacond[m - 1] -= (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                temp = beta[index2(k, m - 1, nbetagroup)] - beta[index2(k, m, nbetagroup)];
                betacond[m - 1] += (minvar[k] * temp < maxvar[k] * temp ? (minvar[k] * temp) : (maxvar[k] * temp));
                //check both prior conditions
                valid = ((theta[m] < (theta[m + 1] + betacond[m])) && (theta[m - 1] < (theta[m] + betacond[m - 1])) ? 1 : 0);
            }
        }
        /*		if(valid==0) printf("Invalid proposal\n");*/
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
                accorig += dnorm((betastore - mnb), 0.0, sdb[index2(k, m, nbetagroup)], 1);
                accprop += dnorm((beta[index2(k, m, nbetagroup)] - mnb), 0.0, sdb[index2(k, m, nbetagroup)], 1);

                //adjust for conditional theta priors
                if(m == 0)
                {
                    accorig -= pnorm(theta[m] - betacondstore1, 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacondstore1, sdt));
                    accprop -= pnorm(theta[m] - betacond[m], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacond[m], sdt));
                }
                else
                {
                    if(m == (ntheta - 1))
                    {
                        accorig -= pnorm(theta[m - 1] - betacondstore, 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore, sdt));
                        accprop -= pnorm(theta[m - 1] - betacond[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
                    }
                    else
                    {
                        accorig -= pnorm(theta[m] - betacondstore1, 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacondstore1, sdt));
                        accprop -= pnorm(theta[m] - betacond[m], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacond[m], sdt));
                        accorig -= pnorm(theta[m - 1] - betacondstore, 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore, sdt));
                        accprop -= pnorm(theta[m - 1] - betacond[m - 1], 0.0, sdt, 0, 1);//log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
                    }
                }
                //proposals cancel
                //calculate acceptance
                acc = accprop - accorig;
                acc = exp(acc);
                if(runif(0.0, 1.0) < acc) (*loglikeorig) = loglikeprop;
                else
                {
                    //reset update to previous value (as well as validity conditions)
                    beta[index2(k, m, nbetagroup)] = betastore;
                    if(m == 0) betacond[m] = betacondstore1;
                    else
                    {
                        if(m == (ntheta - 1)) betacond[m - 1] = betacondstore;
                        else
                        {
                            betacond[m] = betacondstore1;
                            betacond[m - 1] = betacondstore;
                        }
                    }
                }
            }
            else
            {
                //reset update to previous value (as well as validity conditions)
                beta[index2(k, m, nbetagroup)] = betastore;
                if(m == 0) betacond[m] = betacondstore1;
                else
                {
                    if(m == (ntheta - 1)) betacond[m - 1] = betacondstore;
                    else
                    {
                        betacond[m] = betacondstore1;
                        betacond[m - 1] = betacondstore;
                    }
                }
            }
        }
        else
        {
            //reset update to previous value (as well as validity conditions)
            beta[index2(k, m, nbetagroup)] = betastore;
            if(m == 0) betacond[m] = betacondstore1;
            else
            {
                if(m == (ntheta - 1)) betacond[m - 1] = betacondstore;
                else
                {
                    betacond[m] = betacondstore1;
                    betacond[m - 1] = betacondstore;
                }
            }
        }
    }
    return;
}
