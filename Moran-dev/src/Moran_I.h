/*
 * From Guillaume Guénard - Université de Montréal
 * August 2019
 * License: GPL version 3
 * C header
 */

#ifndef __Moran_I_h__

#define __Moran_I_h__
#define INTNUM 25000

#include<R.h>
#include<math.h>

void moran(double* x, double* w, int* n, double* I,
           int* nperm, int* pv, int* retIp, double* Ip_ret);
void permute_idx(int* idx, int n);

#endif
