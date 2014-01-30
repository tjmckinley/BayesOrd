#include <R.h>
#include <Rmath.h>
#include "functions.h"

//function to shuffle integer vector
void ran_shuffle_int(int *vec, int length)
{
	int i, j, k, length1;
	int *prop = (int *) Calloc(length, int);
	
	length1 = length;
	for(j = 0; j < length; j++)
	{
		//sample random value in 0-(length1 - 1)
		i = (int) floor(runif(0.0, (double) length1));
		prop[j] = vec[i];
		for(k = i; k < length1; k++) vec[k] = vec[k + 1];
		length1--;
	}
	for(i = 0; i < length; i++) vec[i] = prop[i];
	//free memory from heap
	Free(prop);
	return;
}
