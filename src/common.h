#ifndef __BTRACK_COMMON_H
#define __BTRACK_COMMON_H

#include <cmath>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <functional>
#include <memory>
#include <cstdlib>
#include <cstring>

#include "SysUtils.hpp"
#include "Math.hpp"
#include "Plan.hpp"

typedef float btrack_chunk_t;
#define BTRACK_STRINGIFY(x) BTRACK_STRINGIFY2(x)
#define BTRACK_STRINGIFY2(x) #x
#define BTRACK_ASSERT(x) if(!(x)){fprintf(stderr,"Error, assertion failed: " __FILE__ " line %d: " BTRACK_STRINGIFY(x) "\n", __LINE__); exit(EXIT_FAILURE);}

#endif
