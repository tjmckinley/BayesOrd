#include <R.h>
#include <Rmath.h>
#include"functions.h"

//function to calculate the log-likelihood component relating to a single RE (psi) term
double loglikelihood_psi(double *beta, double *theta, double *psi, double *variables, int n, int ntheta, int nbetagroup, int *psicount, int currpsi)
{
    int i, j, k;
    double *gamma = (double *) Calloc(ntheta, double);
    double *p = (double *) Calloc((ntheta + 1), double);
    double loglike, mu;

    loglike = 0.0;
    for(i = psicount[currpsi]; i < psicount[currpsi + 1]; i++)
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
    Free(gamma);
    Free(p);
    return loglike;
}
