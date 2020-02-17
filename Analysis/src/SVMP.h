/*************************************************************************
 
 (c) 2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada

 **Core functions**
 
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

 C header files

*************************************************************************/

#ifndef __spectR_h__

#define __spectR_h__

#include<R.h>
#include<math.h>

// MLE-friendly MEM weighting functions

void scf_spher(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res);
void scf_expon(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res);
void scf_power(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res);
void scf_hyper(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res);
void scf_super(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res);

/* Ideally: n should be length 3 and contain size for d, alpha, and beta in
 * order to implement recyling more efficiently than by using a "int* recycle"
 */

// Accessory functions
void dist_Euclid(double* a, double* b, int* na, int* nb,
		int* m, double* d, int* sq);
void get_colmeans(double* x, int n, int m, double* mn_col);
void get_rowmeans(double* x, int n, int m, double* mn_row, double* mn_mat);
void mat_center(double* x, int* n, int* m, double* mn_row,
		        double* mn_col, double* mn_mat, int* rowToo);
void get_center(double* x, int* n, int* m, double* mn_row,
                double* mn_col, double* mn_mat, int* rowToo);
void mat_recenter(double* x, int* n, int* m, double* mn_row,
		          double* mn_col, double* mn_mat, int* rowToo);

#endif
