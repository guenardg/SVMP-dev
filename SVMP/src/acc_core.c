/*************************************************************************
 
 (c) 2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Accessory and core functions**
 
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

 C fonctions definitions

*************************************************************************/

#include"acc_core.h"
#include"matrix.h"

// Accessory functions

/* Computes the Euclidean distance between two matrices.
   The only implementation I could find in R takes a
   single matrix and renders a lower triangular
   distance matrix. That function takes two matrices with
   identical numbers of columns and output a rectangular
   distance matrix.
*/
void dist_Euclid(double* a, double* b, int* na, int* nb,
		 int* m, double* d, int* sq) {
  int i, j, k, os_a, os_b, os_d;
  double acc;
  os_d = 0;
  for(i = 0; i < *na; i++) {
    for(j = 0; j < *nb; j++, os_d++) {
      os_a = i;
      os_b = j;
      for(k = 0; k < *m; k++, os_a += *na, os_b += *nb) {
        acc = a[os_a] - b[os_b];
        acc *= acc;
        d[os_d] += acc;
      }
      if(!*sq)
        d[os_d] = sqrt(d[os_d]);
    }
  }
  return;
}

/*
 * Allows one to center a matrix on its columns (or on both
 * its columns and rows) without prior calculation or
 * the means
 * Expects mn_row, mn_col, and mn_mat to have been
 * initialized to 0.0 by the calling function.
 */

// Columns means calculation, assuming mn_col initialized to 0.0
void get_colmeans(double* x, int n, int m, double* mn_col) {
  int i, j, os;
  os = 0;
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++, os++)
      mn_col[i] += x[os];
    mn_col[i] /= n;
  }
  return;
}

// Row means and overall mean calculation, assuming both are initialized to 0.0
void get_rowmeans(double* x, int n, int m, double* mn_row, double* mn_mat) {
  int i, j, os;
  for(i = 0; i < n; i++) {
    os = i;
    for(j = 0; j < m; j++, os += n)
      mn_row[i] += x[os];
    *mn_mat += mn_row[i];
    mn_row[i] /= m;
  }
  *mn_mat /= (n * m);
  return;
}

void mat_center(double* x, int* n, int* m, double* mn_row,
		            double* mn_col, double* mn_mat, int* rowToo) {
  int i, j, os;
  get_colmeans(x,*n,*m,mn_col); // Calculating column means
  if(*rowToo) {                  // When required: row means and overall mean
    get_rowmeans(x,*n,*m,mn_row,mn_mat);
    // Double centering:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
	      x[os] += (*mn_mat - mn_col[i] - mn_row[j]);
  } else {
    // Columns-only centering:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
	      x[os] -= mn_col[i];
  }
  return;
}

void get_center(double* x, int* n, int* m, double* mn_row,
                double* mn_col, double* mn_mat, int* rowToo) {
  int i, j, os;
  get_colmeans(x,*n,*m,mn_col); // Calculating column means
  if(*rowToo) {                  // When required: row means and overall mean
    get_rowmeans(x,*n,*m,mn_row,mn_mat);
    // Double centering:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
        x[os] = mn_col[i] + mn_row[j] - *mn_mat;
  } else {
    // Columns-only centering:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
        x[os] = mn_col[i];
  }
  return;
}

void mat_recenter(double*x, int* n, int* m, double* mn_row,
		              double* mn_col, double* mn_mat, int* rowToo) {
  int i, j, os;
  if(*rowToo) {
    // When necessary row means and the overall mean
    for(i = 0; i < *n; i++) {
      os = i;
      for(j = 0; j < *m; j++, os += *n)
        mn_row[i] += x[os];
      mn_row[i] /= *m;
    }
    // Double centering here:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
        x[os] += (*mn_mat - mn_col[i] - mn_row[j]);
  } else {
    // Columns-only centering there:
    os = 0;
    for(i = 0; i < *m; i++)
      for(j = 0; j < *n; j++, os++)
        x[os] -= mn_col[i];
  }
  return;
}
