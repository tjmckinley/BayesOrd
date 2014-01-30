/*Metropolis-Hastings routine for fitting ordinal regression model*/
/*Written for running multiple chains with multicore in R*/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "functions.h"

//main MCMC function
SEXP bayesord(SEXP intinputs, SEXP doubleinputs, SEXP Rvariables, SEXP Rxassign, SEXP Rintfactor, SEXP Rinteraction, SEXP Rintstart, SEXP Rpsicount)
{
	//declare counters etc.
    int i, j, k, m;
    double u;

    //import initial conditions and priors into function
    //and set up parameters to be estimated
    int chain, n, niter, nvariables, nbeta, nbetagroup, ntheta, noutputsum, npsi, RE, movetype, varselect, fixed, maxinteraction, runtraining, nitertrain, printall, multi;
    double mnb, varb, maxsdb, vart, propsdb, propsdt, mnpsi, shvarp, rtvarp, propsdp, propsdvarp;
    
    //check input lengths
    if(LENGTH(intinputs) != 18 && LENGTH(doubleinputs) != 11)
    {
        Rprintf("Input information of incorrect length\n");
        return R_NilValue;
    }

    //set variables from inputs
    intinputs = PROTECT(coerceVector(intinputs, INTSXP));
    doubleinputs = PROTECT(coerceVector(doubleinputs, REALSXP));
    
    chain = INTEGER(intinputs)[0];
    n = INTEGER(intinputs)[1];
    niter = INTEGER(intinputs)[2];
    noutputsum = INTEGER(intinputs)[3];
    nvariables = INTEGER(intinputs)[4];
    nbeta = INTEGER(intinputs)[5];
    nbetagroup = INTEGER(intinputs)[6];
    ntheta = INTEGER(intinputs)[7];
    fixed = INTEGER(intinputs)[8];
    npsi = INTEGER(intinputs)[9];
    RE = INTEGER(intinputs)[10];
    movetype = INTEGER(intinputs)[11];
    varselect = INTEGER(intinputs)[12];
    maxinteraction = INTEGER(intinputs)[13];
    runtraining = INTEGER(intinputs)[14];
    nitertrain = INTEGER(intinputs)[15];
    printall = INTEGER(intinputs)[16];
    multi = INTEGER(intinputs)[17];
    
    mnb = REAL(doubleinputs)[0];
    varb = REAL(doubleinputs)[1];
    vart = REAL(doubleinputs)[2];
    mnpsi = REAL(doubleinputs)[3];
    shvarp = REAL(doubleinputs)[4];
    rtvarp = REAL(doubleinputs)[5];
    propsdb = REAL(doubleinputs)[6];
    propsdt = REAL(doubleinputs)[7];
    propsdp = REAL(doubleinputs)[8];
    propsdvarp = REAL(doubleinputs)[9];
    maxsdb = REAL(doubleinputs)[10];
    
    //set up data structures
    Rvariables = PROTECT(coerceVector(Rvariables, REALSXP));
    Rxassign = PROTECT(coerceVector(Rxassign, INTSXP));
    Rintfactor = PROTECT(coerceVector(Rintfactor, INTSXP));
    Rinteraction = PROTECT(coerceVector(Rinteraction, INTSXP));
    Rintstart = PROTECT(coerceVector(Rintstart, INTSXP));
    Rpsicount = PROTECT(coerceVector(Rpsicount, INTSXP));
    //set pointers to data structures for succinct coding in C
    double *variables = REAL(Rvariables);
    int *xassign = INTEGER(Rxassign);
    int *intfactor = INTEGER(Rintfactor);
    int *interaction = INTEGER(Rinteraction);
    int *intstart = INTEGER(Rintstart);
    int *psicount = INTEGER(Rpsicount);

    if(printall == 1)
    {
        //visual check of input information
        Rprintf("\n#### Data ####\n");
        Rprintf("number of data points = %d\n", n);
        Rprintf("number of iterations = %d\n", niter);
        Rprintf("output chain summaries after every %d iterations\n", noutputsum);

        Rprintf("\n#### Parameters ####\n");
        Rprintf("number of variables = %d\n", nvariables);
        Rprintf("number of beta parameters overall = %d\n", nbeta);
        Rprintf("number of beta parameters in each group = %d\n", nbetagroup);
        Rprintf("number of theta parameters = %d\n", ntheta);

        Rprintf("\n#### Priors ####\n");
        Rprintf("mean beta = %f\n", mnb);
        if(fixed == 1) Rprintf("var beta = %f\n", varb);
        else Rprintf("SD beta ~ N(0,sigma^2), where sigma ~ U(0,%f)\n", maxsdb);
        Rprintf("var theta = %f\n", vart);
        if(RE == 1)
        {
            Rprintf("mn psi = %f\n", mnpsi);
            Rprintf("shape varp = %f\n", shvarp);
            Rprintf("rate varp = %f\n", rtvarp);
            Rprintf("mean varp = %f\tvariance varp = %f\n", shvarp / rtvarp, shvarp / (rtvarp * rtvarp));
        }

        Rprintf("\n#### Proposal ####\n");
        Rprintf("proposal SD for beta = %f\n", propsdb);
        Rprintf("max. jump for theta = %f\n", propsdt);
        if(RE == 1)
        {
            Rprintf("max. jump for psi=%f\n", propsdp);
            Rprintf("max. jump for varp=%f\n", propsdvarp);
        }

        if(movetype == 0) Rprintf("\n#### PROPORTIONAL ODDS MODEL ####\n");
        else
        {
            if(movetype == 1) Rprintf("\n#### NON-PROPORTIONAL ODDS MODEL ####\n");
            else Rprintf("\n#### PROPORTIONAL ODDS OR NON PROPORTIONAL ODDS MODEL ####\n");
        }
        if(RE == 0) Rprintf("Model includes NO random effect terms\n");
        else Rprintf("Model includes RANDOM EFFECT terms\n");
        if(varselect == 1) Rprintf("Model includes VARIABLE SELECTION\n");
        if(runtraining == 1) Rprintf("\nTraining run being used for %d iterations\n", nitertrain);
    }

    //visual check of data
    if(printall == 1)
    {
        Rprintf("\n#### Visual check of data ####\n");
        for(i = 0; i < 5; i++)
        {
            for(j = 0; j < (nbetagroup + 2); j++) Rprintf("%f\t", variables[index2(i, j, n)]);
            if(RE == 1) Rprintf("%f", variables[index2(i, j, n)]);
	        Rprintf("\n");
	    }
	}

    //set number of parameters
    int npars = nbeta + ntheta + npsi + 1 + nbetagroup;
    if(fixed == 0) npars += nbeta;
    
    //declare output vector and pointer for saving out MCMC output
    SEXP Rposterior;
    PROTECT(Rposterior = allocVector(REALSXP, (npars + 1) * niter));
    double *posterior = REAL(Rposterior);
    for(j = 0; j < ((npars + 1) * niter); j++) posterior[j] = 0;
    //declare vectors relating to MCMC
    int kcoda = 0;
    double *beta = (double *) Calloc(nbeta, double);
    double *theta = (double *) Calloc(ntheta, double);
    double *psi = (double *) Calloc(npsi, double);
    double varp = 1.0;
    //set up conditions for beta parameters
    double *betacond = (double *) Calloc((ntheta - 1), double);
    double *betacondstore_v = (double *) Calloc((ntheta - 1), double);
    for(i = 0; i < (ntheta - 1); i++)
    {
        betacond[i] = 0.0;
        betacondstore_v[i] = 0.0;
    }
    //set minimum and maximum values for each variable based on the data
    double *minvar = (double *) Calloc(nbetagroup, double);
    double *maxvar = (double *) Calloc(nbetagroup, double);
    for(j = 0; j < nbetagroup; j++)
    {
        minvar[j] = variables[index2(0, j, n)];
        maxvar[j] = variables[index2(0, j, n)];
        for(i = 1; i < n; i++)
        {
            minvar[j] = (minvar[j] < variables[index2(i, j, n)] ? minvar[j] : variables[index2(i, j, n)]);
            maxvar[j] = (maxvar[j] > variables[index2(i, j, n)] ? maxvar[j] : variables[index2(i, j, n)]);
        }
    }
    for(i = 0; i < nbetagroup; i++)
    {
        if(minvar[i] == maxvar[i])
        {
            Rprintf("No variation in %dth variable\n", i);
            return R_NilValue;
        }
    }
    double temp = 0.0;

    //set up temporary storage variables
    double psistore, varpstore;
    int valid;
    //set up vector for shuffling when updating beta terms
    int *betashuffle = (int *) Calloc(nvariables, int);
    for(i = 0; i < nvariables; i++) betashuffle[i] = i;
    int *betashufflepar = (int *) Calloc(nbetagroup, int);
    for(i = 0; i < nbetagroup; i++) betashufflepar[i] = i;
    int *betashuffletheta = (int *) Calloc(ntheta, int);
    for(i = 0; i < ntheta; i++) betashuffletheta[i] = i;
    //set up vector for recording parameter status (0=PO, 1=NPO)
    int *betastatus = (int *) Calloc(nbetagroup, int);
    for(i = 0; i < nbetagroup; i++) betastatus[i] = 1;
    //set up vector for recording variable status (0=exc, 1=inc)
    int *betastatusvar = (int *) Calloc(nvariables, int);
    for(i = 0; i < nvariables; i++) betastatusvar[i] = 1;
    
    if(maxinteraction > 0)
    {
        //visual check of data inputs
        if(printall == 1)
        {
            Rprintf("intfactor:\n");
            for(i = 0; i < nvariables; i++)
            {
                for(j = 0; j < nvariables; j++) Rprintf("%d ", intfactor[index2(i, j, nvariables)]);
                Rprintf("\n");
            }
            Rprintf("\ninteraction:\n");
            for(i = 0; i < nvariables; i++) Rprintf("%d ", interaction[i]);
            Rprintf("\n");
            Rprintf("\nintstart:\n");
            for(i = 0; i < (maxinteraction + 2); i++) Rprintf("%d ", intstart[i]);
            Rprintf("\n");
        }
    }

    //set up loglikelihood stuff
    double loglikeorig, loglikeprop, accorig, accprop, acc;

    double *sdb = (double *) Calloc(nbeta, double);
    for(i = 0; i < nbeta; i++) sdb[i] = sqrt(varb);
    double sdt = sqrt(vart);
    double sdp;
    //sample initial values
    m = 0;
    loglikeorig = 0.0 / 0.0;
    while(isfinite(loglikeorig) == 0 && m < 1000)
    {
        //reset parameters
        for(i = 0; i < nbetagroup; i++) betastatus[i] = 0;
        for(i = 0; i < nbeta; i++) beta[i] = 0.0;
        //sample values for beta
        for(i = 0; i < nbetagroup; i++)
        {
            if(betastatus[i] != 2) beta[index2(i, 0, nbetagroup)] = mnb + rnorm(0.0, 1.0);
            else beta[index2(i, 0, nbetagroup)] = 0.0;
        }
        if(movetype == 0 || runtraining == 1)
        {
            for(i = 0; i < nbetagroup; i++)
            {
                if(betastatus[i] != 2)
                {
                    for(j = 1; j < ntheta; j++) beta[index2(i, j, nbetagroup)] = beta[index2(i, j - 1, nbetagroup)];
                    betastatus[i] = 0;
                }
            }
        }
        else
        {
            if(movetype == 1)
            {
                for(i = 0; i < nbetagroup; i++)
                {
                    if(betastatus[i] != 2)
                    {
                        for(j = 1; j < ntheta; j++) beta[index2(i, j, nbetagroup)] = beta[index2(i, j - 1, nbetagroup)] - runif(0.0, 1.0);
                        betastatus[i] = 1;
                    }
                }
            }
            else
            {
                for(i = 0; i < nbetagroup; i++)
                {
                    if(betastatus[i] != 2)
                    {
                        betastatus[i] = 0;
                        if(betastatus[i] == 0)
                        {
                            for(j = 1; j < ntheta; j++)  beta[index2(i, j, nbetagroup)] = beta[index2(i, j - 1, nbetagroup)];
                        }
                        else
                        {
                            for(j = 1; j < ntheta; j++)  beta[index2(i, j, nbetagroup)] = beta[index2(i, j - 1, nbetagroup)] - runif(0.0, 1.0);
                        }
                    }
                }
            }
        }
        //set and check prior conditions
        for(i = 0; i < (ntheta - 1); i++)
        {
            betacond[i] = 0.0;
            for(j = 0; j < nbetagroup; j++)
            {
                temp = beta[index2(j, i, nbetagroup)] - beta[index2(j, i + 1, nbetagroup)];
                betacond[i] += ((minvar[j] * temp) < (maxvar[j] * temp) ? (minvar[j] * temp) : (maxvar[j] * temp));
            }
        }
        //randomly sample cut-points
        theta[0] = rnorm(0.0, 1.0);
        for(j = 1; j < ntheta; j++)
        {
            theta[j] = rnorm(0.0, 1.0);
            while(theta[j] < theta[j - 1]) theta[j] = rnorm(0.0, 1.0);
        }
        //check validity
        valid = 1;
        for(j = 0; j < (ntheta - 1); j++)
        {
            if(theta[j] >= (theta[j + 1] + betacond[j])) valid = 0;
        }
        //sample hyperparameters
        if(fixed == 0)
        {
            for(j = 0; j < nbeta; j++) sdb[j] = runif(0.0, maxsdb);
        }
        if(valid == 1)
        {
            //sample random effects
            if(RE == 1)
            {
                varp = rgamma(1.0, 1.0);
                for(j = 0; j < npsi; j++) psi[j] = mnpsi + rnorm(0.0, sqrt(varp));
            }
            else
            {
                varp = 1.0;
                for(j = 0; j < npsi; j++) psi[j] = 0.0;
            }
            //calculate log-likelihood
            loglikeorig = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);
        }
        m++;
    }

    if(m == 1000) {
        if(multi == 1 && chain > 1) Rprintf("Can't initialise: chain %d!\n", chain);
        else Rprintf("\nCan't initialise: chain %d!\n", chain);
        return R_NilValue;
    }
    if(multi == 1 && chain > 1) Rprintf("System initialised OK: chain %d!\n", chain);
    else Rprintf("\nSystem initialised OK: chain %d!\n", chain);

    //set up and read in cumulative count vector for efficient psi updates
/*    if(printall == 1) Rprintf("RE = %d\n", RE);*/
    if(RE == 1)
    {
        if(printall == 1)
        {
            //visual check
            for(j = 0; j < 5; j++) Rprintf("psicount[%d]=%d\n", j, psicount[j]);
            Rprintf("...\n");
            for(j = npsi - 5; j < (npsi + 1); j++) Rprintf("psicount[%d]=%d\n", j, psicount[j]);
            Rprintf("\n");
        }
    }

    //function to check validity of proposals
    if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;

    double pinteraction;
    //set move probabilities (0: keep, 1: add/remove)
    double pvecvar[2];
    if(varselect == 1) pvecvar[0] = 0.5;
    else pvecvar[0] = 1.0;
    pvecvar[1] = 1.0;
    double pvecvarmove[2];
    pvecvarmove[0] = 0.0;
    pvecvarmove[1] = 1.0;

    //set move probabilities (0: move, 1: PO/NPO change)
    double pvec[2];
    if(movetype == 0 || movetype == 1)
    {
        pvec[0] = 1.0;
        pvec[1] = 1.0;
    }
    else
    {
        pvec[0] = 0.5;
        pvec[1] = 1.0;
    }
    //set zero to PO and NPO move probabilities
    double pzero[2];
    if(varselect == 1)
    {
        if(movetype == 0) pzero[0] = 1.0;
        else
        {
            if(movetype == 1)	pzero[0] = 0.0;
            else pzero[0] = 0.5;
        }
        pzero[1] = 1.0;
    }

    //start MCMC chain
    if(runtraining == 1)
    {
        for(i = 1; i < nitertrain; i++)
        {
            //update linear terms in random order
            ran_shuffle_int(betashufflepar, nbetagroup);
            for(j = 0; j < nbetagroup; j++)
            {
                k = betashufflepar[j];
                movePO(k, n, nbetagroup, ntheta, beta, theta, psi, variables, propsdb, &loglikeorig, mnb, sdb);
            }
            //update theta terms in random order
            ran_shuffle_int(betashuffletheta, ntheta);
            for(m = 0; m < ntheta; m++) movetheta(betashuffletheta[m], n, nbetagroup, ntheta, beta, theta, psi, variables, betacond, propsdt, &loglikeorig, sdt);
            //now check validity conditions
            if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
            //move SD for beta parameter
            if(fixed == 0)
            {
                //update SD terms in random order
                ran_shuffle_int(betashufflepar, nbetagroup);
                for(j = 0; j < nbetagroup; j++)
                {
                    k = betashufflepar[j];
                    movesdb(k, 0, nbetagroup, beta, propsdb, mnb, sdb, maxsdb);
                }
            }
        }
        if(movetype == 1) for(j = 0; j < nbetagroup; j++) betastatus[j] = 1;
        Rprintf("Finished training run, restarting chains...\n");
    }

    //save posterior for initial samples
    k = 0;
    for(j = 0; j < nbeta; j++) {
        posterior[index2(kcoda, j, niter)] = beta[j];
    }
    k = nbeta;
    for(j = 0; j < ntheta; j++) {
        posterior[index2(kcoda, j + k, niter)] = theta[j];
    }
    k += ntheta;
    for(j = 0; j < nbetagroup; j++) {
        posterior[index2(kcoda, j + k, niter)] = betastatus[j];
    }
    k += nbetagroup;
    for(j = 0; j < npsi; j++) {
        posterior[index2(kcoda, j + k, niter)] = psi[j];
    }
    k += npsi;
    if(fixed == 0)
    {
        for(j = 0; j < nbeta; j++) {
            posterior[index2(kcoda, j + k, niter)] = sdb[j];
        }
        k += nbeta;
    }
    posterior[index2(kcoda, k, niter)] = varp;
    posterior[index2(kcoda, k + 1, niter)] = loglikeorig;
    kcoda++;

    for(i = 1; i < niter; i++)
    {
        //update linear terms in random order
        ran_shuffle_int(betashuffle, nvariables);
        if(varselect == 1)
        {
            for(j = 0; j < nvariables; j++)
            {
                //select variable to update from shuffled set
                k = betashuffle[j];
                //if variable is currently excluded, then propose to keep or add
                if(betastatusvar[k] == 0)
                {
                    //check whether variable selection can occur and if so check and amend probabilities for interaction effects
                    if(interaction[k] > 0)
                    {
                        //check that earlier main and interaction effects are included if we wish to add this variable
                        pinteraction = pvecvar[0];
                        for(m = 0; m < nvariables; m++)
                        {
                            if(intfactor[index2(m, k, nbetagroup)] == 1 && betastatusvar[m] == 2) pinteraction = 1.0;
                        }
                        pvecvarmove[0] = pinteraction;
                        //if variable is currently excluded, then propose to keep or add
                        if(runif(0.0, 1.0) > pvecvar[0])
                        {
                            betastatusvar[k] = 1;
                            //add variable to model choosing PO or NPO status for each parameter randomly
                            for(m = xassign[k]; m < xassign[k + 1]; m++)
                            {
                                if(runif(0.0, 1.0) < pzero[0]) betastatus[m] = 0;
                                else betastatus[m] = 1;
                            }
                            moveexctoinc(k, xassign, betastatusvar, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, pvecvarmove, pzero, maxsdb, fixed);
                        }
                        //now check validity conditions
                        if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
                    }
                    else
                    {
                        if(runif(0.0, 1.0) > pvecvar[0])
                        {
                            betastatusvar[k] = 1;
                            //add variable to model choosing PO or NPO status for each parameter randomly
                            for(m = xassign[k]; m < xassign[k + 1]; m++)
                            {
                                if(runif(0.0, 1.0) < pzero[0]) betastatus[m] = 0;
                                else betastatus[m] = 1;
                            }
                            moveexctoinc(k, xassign, betastatusvar, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, pvecvar, pzero, maxsdb, fixed);
                        }
                        //now check validity conditions
                        if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
                    }
                }
                else
                {
                    //check whether variable selection can occur and if so check and amend probabilities for interaction effects
                    if(maxinteraction > 0)
                    {
                        //check that later interaction effects are not currently included if we wish to drop this variable
                        pinteraction = pvecvar[0];
                        if(maxinteraction > interaction[k])
                        {
                            for(m = intstart[interaction[k] + 1]; m < intstart[maxinteraction + 1]; m++)
                            {
                                if(intfactor[index2(k, m, nvariables)] == 1 && betastatusvar[m] != 2) pinteraction = 1.0;
                            }
                        }
                        pvecvarmove[0] = pinteraction;
                        //if variable is currently included, then propose to keep or drop
                        if(runif(0.0, 1.0) > pvecvarmove[0])
                        {
                            betastatusvar[k] = 0;
                            //drop variable from model
                            moveinctoexc(k, xassign, betastatusvar, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, pvecvarmove, pzero, maxsdb, fixed);
                        }
                        //now check validity conditions
                        if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
                    }
                    else
                    {
                        //if variable is currently included, then propose to keep or drop
                        if(runif(0.0, 1.0) > pvecvar[0])
                        {
                            betastatusvar[k] = 0;
                            //drop variable from model
                            moveinctoexc(k, xassign, betastatusvar, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, pvecvar, pzero, maxsdb, fixed);
                        }
                        //now check validity conditions
                        if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
                    }
                }
            }
        }
        //update linear terms in random order
        ran_shuffle_int(betashufflepar, nbetagroup);
        for(j = 0; j < nbetagroup; j++)
        {
            k = betashufflepar[j];
            //if variable currently included in model then update all parameters associated with it
            if(betastatus[k] != 2)
            {
                //if parameter is currently PO, then choose to move parameter or shift to NPO
                if(betastatus[k] == 0)
                {
                    u = runif(0.0, 1.0);
                    if(u < pvec[0]) movePO(k, n, nbetagroup, ntheta, beta, theta, psi, variables, propsdb, &loglikeorig, mnb, sdb);
                    else movePOtoNPO(k, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, maxsdb, fixed);
                }
                else
                {
                    //if parameter is currently NPO, then choose to move parameter or shift to PO
                    u = runif(0.0, 1.0);
                    if(u < pvec[0]) moveNPO(k, n, nbetagroup, ntheta, beta, theta, psi, variables, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt);
                    else moveNPOtoPO(k, n, nbetagroup, ntheta, beta, theta, psi, variables, betastatus, betacond, minvar, maxvar, propsdb, &loglikeorig, mnb, sdb, sdt, maxsdb, fixed);
                }
                //now check validity conditions
                if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;
            }
        }
        //update theta terms in random order
        ran_shuffle_int(betashuffletheta, ntheta);
        for(m = 0; m < ntheta; m++) movetheta(betashuffletheta[m], n, nbetagroup, ntheta, beta, theta, psi, variables, betacond, propsdt, &loglikeorig, sdt);
        //now check validity conditions
        if(validitycheck(ntheta, nbetagroup, betastatus, beta, theta, betacond, minvar, maxvar) == 1) return R_NilValue;

        //move SD for beta parameter
        if(fixed == 0)
        {
            //update SD terms in random order
            ran_shuffle_int(betashufflepar, nbetagroup);
            for(j = 0; j < nbetagroup; j++)
            {
                k = betashufflepar[j];
                if(betastatus[k] == 0) movesdb(k, 0, nbetagroup, beta, propsdb, mnb, sdb, maxsdb);
                else
                {
                    if(betastatus[k] == 1)
                    {
                        ran_shuffle_int(betashuffletheta, ntheta);
                        for(m = 0; m < ntheta; m++) movesdb(k, betashuffletheta[m], nbetagroup, beta, propsdb, mnb, sdb, maxsdb);
                    }
                }
            }
        }
        if(RE == 1)
        {
            //randomly update individual effect terms
            sdp = sqrt(varp);
            for(j = 0; j < ((int) floor(0.3 * npsi)); j++)
            {
                k = (int) floor(npsi * runif(0.0, 1.0));
                psistore = psi[k];
                loglikeorig = loglikelihood_psi(beta, theta, psi, variables, n, ntheta, nbetagroup, psicount, k);
                psi[k] = psistore + runif(-propsdp, propsdp);
                loglikeprop = loglikelihood_psi(beta, theta, psi, variables, n, ntheta, nbetagroup, psicount, k);
                //work out acceptance probability
                if(isfinite(loglikeprop) != 0)
                {
                    accorig = loglikeorig;
                    accprop = loglikeprop;
                    //adjust for priors
                    accorig += dnorm((psistore - mnpsi), 0.0, sdp, 1);
                    accprop += dnorm((psi[k] - mnpsi), 0.0, sdp, 1);
                    //proposals cancel
                    acc = accprop - accorig;
                    acc = exp(acc);
                    if(runif(0.0, 1.0) >= acc) psi[k] = psistore;
                }
                else psi[k] = psistore;
            }

            //update full likelihood
            loglikeorig = loglikelihood(beta, theta, psi, variables, n, ntheta, nbetagroup);

            //update individual effect variance term
            varpstore = varp;
            varp += runif(-propsdvarp, propsdvarp);
            if(varp > 0)
            {
                accorig = 0.0;
                accprop = 0.0;
                //adjust for priors
                sdp = sqrt(varpstore);
                for(j = 0; j < npsi; j++) accorig += dnorm((psi[j] - mnpsi), 0.0, sdp, 1);
                sdp = sqrt(varp);
                for(j = 0; j < npsi; j++) accprop += dnorm((psi[j] - mnpsi), 0.0, sdp, 1);
                accorig += dgamma(varpstore, shvarp, 1.0 / rtvarp, 1);
                accprop += dgamma(varp, shvarp, 1.0 / rtvarp, 1);
                //proposals cancel
                acc = accprop - accorig;
                acc = exp(acc);
                if(runif(0.0, 1.0) >= acc) varp = varpstore;
            }
            else varp = varpstore;
        }

        //record output
        k = 0;
        for(j = 0; j < nbeta; j++) {
            posterior[index2(kcoda, j, niter)] = beta[j];
        }
        k = nbeta;
        for(j = 0; j < ntheta; j++) {
            posterior[index2(kcoda, j + k, niter)] = theta[j];
        }
        k += ntheta;
        for(j = 0; j < nbetagroup; j++) {
            posterior[index2(kcoda, j + k, niter)] = betastatus[j];
        }
        k += nbetagroup;
        for(j = 0; j < npsi; j++) {
            posterior[index2(kcoda, j + k, niter)] = psi[j];
        }
        k += npsi;
        if(fixed == 0)
        {
            for(j = 0; j < nbeta; j++) {
                posterior[index2(kcoda, j + k, niter)] = sdb[j];
            }
            k += nbeta;
        }
        posterior[index2(kcoda, k, niter)] = varp;
        posterior[index2(kcoda, k + 1, niter)] = loglikeorig;
        m = k + 1;
        kcoda++;

        //save output into text files for examination/editing in R
        if((i + 1) % noutputsum == 0) Rprintf("Chain %d: current iterations=%d\n", chain, i + 1);
    }
    
    //free memory from the heap
    Free(beta); Free(theta); Free(psi);
    Free(betacond); Free(betacondstore_v); Free(minvar); Free(maxvar);
    Free(betashuffle); Free(betashufflepar); Free(betashuffletheta);
    Free(betastatus); Free(betastatusvar); Free(sdb);
    
    UNPROTECT(9);
    
    // terminate the program:
    return Rposterior;
}
