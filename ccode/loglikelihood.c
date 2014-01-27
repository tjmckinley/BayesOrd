#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include"functions.h"

//function to calculate the log-likelihood
double loglikelihood(double *beta, double *theta, double *psi, double *variables, int n, int ntheta, int nbetagroup)
{
    int i, j, k;
    double *gamma = (double *) malloc(ntheta * sizeof(double));
    double *p = (double *) malloc((ntheta + 1) * sizeof(double));
    double loglike, mu;

    loglike = 0.0;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < ntheta; j++)
        {
            mu = 0.0;
            for(k = 0; k < nbetagroup; k++) mu += beta[index2(k, j, nbetagroup)] * variables[index2(i, k, n)];
            mu += psi[(int) (variables[index2(i, nbetagroup + 2, n)] - 1)];
            gamma[j] = theta[j] - mu;
            gamma[j] = exp(gamma[j]) / (1.0 + exp(gamma[j]));
        }
        p[0] = gamma[0];
        for(j = 1; j < ntheta; j++) p[j] = gamma[j] - gamma[j - 1];
        p[ntheta] = 1.0 - gamma[ntheta - 1];
        loglike += variables[index2(i, nbetagroup, n)] * log(p[(int) variables[index2(i, nbetagroup + 1, n)]]);
    }
    free(gamma);
    free(p);
    return loglike;
}
