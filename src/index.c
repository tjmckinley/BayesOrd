#include <R.h>
#include "functions.h"

//indexing function for matrices
int index2(int i, int j, int nrow)
{
    return (nrow * j) + i;
}
