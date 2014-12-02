#ifndef __CISI__
#define __CISI__

#include "common.h"
#include "math.h"
#include "float.h"
#include "complex.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int cisi (
        double x,
        double * ci,
        double * si,
        ErrorMsg errmsg
        );

#ifdef __cplusplus
}
#endif

#endif