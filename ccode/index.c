#include <stdlib.h>
#include"functions.h"

//indexing function for matrices
inline int index2(int i, int j, int nrow)
{
	return (nrow*j)+i;
}
