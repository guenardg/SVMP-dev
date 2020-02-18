/*************************************************************************
 
 (c) 2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Geodesic distances**
 
 This file is part of SVMP
 
 SVMP is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SVMP is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SVMP.  If not, see <https://www.gnu.org/licenses/>.
 
 C header
 
 *************************************************************************/

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
