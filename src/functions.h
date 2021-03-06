#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

//indexing function for matrices
int index2(int i, int j, int nrow);

//main MCMC function
SEXP bayesord(SEXP intinputs, SEXP doubleinputs, SEXP Rvariables, SEXP Rxassign, SEXP Rintfactor, SEXP Rinteraction, SEXP Rintstart, SEXP Rpsicount);

//function to shuffle integer vector
void ran_shuffle_int(int *vec, int length);

//function to calculate the log-likelihood
double loglikelihood(double *beta, double *theta, double *psi, double *variables, int n, int ntheta, int nbetagroup);

//function to calculate the log-likelihood component relating to a single RE (psi) term
double loglikelihood_psi(double *beta, double *theta, double *psi, double *variables, int n, int ntheta, int nbetagroup, int *psicount, int currpsi);

//function to check validity of proposals
int validitycheck(int ntheta, int nbetagroup, int *betastatus, double *beta, double *theta, double *betacond, double *minvar, double *maxvar);

//function to move PO parameter
void movePO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double propsdb, double *loglikeorig, double mnb, double *sdb);

//function to move PO parameter to NPO parameter
void movePOtoNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double maxsdb, int fixed);

//function to move NPO parameter to PO parameter
void moveNPOtoPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double maxsdb, int fixed);

//function to move NPO parameters
void moveNPO(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt);

//function to drop parameters for variable currently present in model
void moveinctoexc(int k, int *xassign, int *betastatusvar, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double *pvec, double *pzero, double maxsdb, int fixed);

//function to add parameters for variable not currently present in model
void moveexctoinc(int k, int *xassign, int *betastatusvar, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, int *betastatus, double *betacond, double *minvar, double *maxvar, double propsdb, double *loglikeorig, double mnb, double *sdb, double sdt, double *pvec, double *pzero, double maxsdb, int fixed);

//function to move theta parameter
void movetheta(int k, int n, int nbetagroup, int ntheta, double *beta, double *theta, double *psi, double *variables, double *betacond, double propsdt, double *loglikeorig, double sdt);

//function to move SD for beta parameter
void movesdb(int k, int m, int nbetagroup, double *beta, double propsdb, double mnb, double *sdb, double maxsdb);
