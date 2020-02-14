/*
 * From Guillaume Guénard - Université de Montréal
 * August 2019
 * License: GPL version 3
 * C header
 */

#ifndef __geodesics_h__

#define __geodesics_h__
#define INTERRUPT_CHECK 10000

#include<R.h>
#include<math.h>

void dist_geo_hvs(double* from, double* to, int* n, int* tri,
                  double* d, double* r);

void dist_geo_vif(double* from, double* to, int* n, int* tri,
                  double* d, int* niter,
                  double* a, double* f, int* maxiter, double* tol);

#endif
