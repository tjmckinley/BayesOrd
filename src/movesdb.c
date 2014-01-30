#include <R.h>
#include <Rmath.h>
#include"functions.h"

//function to move SD for beta parameter
void movesdb(int k, int m, int nbetagroup, double *beta, double propsdb, double mnb, double *sdb, double maxsdb)
{
    double sdbprop;
    double accorig, accprop, acc;

    sdbprop = sdb[index2(k, m, nbetagroup)] + runif(-propsdb, propsdb);

    if(sdbprop > 0.0 && sdbprop < maxsdb)
    {
        accorig = 0.0;
        accprop = 0.0;

        //adjust for beta priors
        accorig += dnorm((beta[index2(k, m, nbetagroup)] - mnb), 0.0, sdb[index2(k, m, nbetagroup)], 1);
        accprop += dnorm((beta[index2(k, m, nbetagroup)] - mnb), 0.0, sdbprop, 1);
        //priors and proposals cancel for SD terms

        //calculate acceptance
        acc = accprop - accorig;
        acc = exp(acc);
        if(runif(0.0, 1.0) < acc) sdb[index2(k, m, nbetagroup)] = sdbprop;
    }
    return;
}
