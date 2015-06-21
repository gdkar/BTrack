#include <fftw3.h>
#include "sse_mathfun.h"
#include <stdlib.h>
#include <stdbool.h>

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef SQR
#define SQR(x)    ((x)*(x))
#endif
#ifndef ASSERT
#define ASSERT(x) if(!x){fprintf(stderr,"Error, assertion failed: " __FILE__ " line %d.\n", __LINE__); exit(EXIT_FAILURE);}
#endif
