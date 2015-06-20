#include <fftw3.h>
#include "sse_mathfun.h"
#include <stdlib.h>
#include <stdbool.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SQR(x)    ((x)*(x))
#define ASSERT(x) if(!x){fprintf(stderr,"Error, assertion failed: " __FILE__ " line %d.\n", __LINE__); exit(EXIT_FAILURE);}
