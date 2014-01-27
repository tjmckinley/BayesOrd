#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "functions.h"

//function to move NPO parameters
void moveNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, gsl_rng *rand_gen, double mnb, double *sdb, double sdt)
{
    int m, valid;
    double betastore, betacondstore = 0.0, betacondstore1 = 0.0, loglikeprop, temp;
    double accorig, accprop, acc;

    //update each regression parameter in turn for the given covariate
    for(m = 0; m < ntheta; m++)
    {
        betastore = beta[index2(k, m, nbetagroup)];
        //sample new parameter value
        beta[index2(k, m, nbetagroup)] = betastore + gsl_ran_flat(rand_gen, -propsdb, propsdb);

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
                accorig += log(gsl_ran_gaussian_pdf((betastore - mnb), sdb[index2(k, m, nbetagroup)]));
                accprop += log(gsl_ran_gaussian_pdf((beta[index2(k, m, nbetagroup)] - mnb), sdb[index2(k, m, nbetagroup)]));

                //adjust for conditional theta priors
                if(m == 0)
                {
                    accorig -= log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacondstore1, sdt));
                    accprop -= log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacond[m], sdt));
                }
                else
                {
                    if(m == (ntheta - 1))
                    {
                        accorig -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore, sdt));
                        accprop -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
                    }
                    else
                    {
                        accorig -= log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacondstore1, sdt));
                        accprop -= log(1.0 - gsl_cdf_gaussian_P(theta[m] - betacond[m], sdt));
                        accorig -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacondstore, sdt));
                        accprop -= log(1.0 - gsl_cdf_gaussian_P(theta[m - 1] - betacond[m - 1], sdt));
                    }
                }
                //proposals cancel
                //calculate acceptance
                acc = accprop - accorig;
                acc = exp(acc);
                if(gsl_rng_uniform_pos(rand_gen) < acc) (*loglikeorig) = loglikeprop;
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
