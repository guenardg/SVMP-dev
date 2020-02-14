/*************************************************************************

    (c) 2019 Guillaume Guénard
        Univeersité de Montréal, Montreal, Quebec, Canada

    This file is part of spectR

    spectR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    spectR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <https://www.gnu.org/licenses/>.

    C header files

*************************************************************************/

#ifndef __spectR_h__

#define __spectR_h__
#include<R.h>
#include<math.h>

// Simple structure and functions for linear algebra

// A structure to wrap arrays into matrices.
typedef struct matrix {
    char* id;           // Matrix idenfification string.
    unsigned int nr;    // Number of row(s) of the matrix.
    unsigned int nc;    // Number of column(s) of the matrix.
    double* v;          // Values attached to the matrix, ordered by column(s).
} matrix;

matrix initmatrix(char* id, unsigned int nr, unsigned int nc);
matrix assignmatrix(char* id, unsigned int nr, unsigned int nc, double* v);
void freematrix(matrix* mat);
void deassignmatrix(matrix* mat);
matrix copymatrix(matrix* a);
void rowsums(matrix* a, double* s);
void colsums(matrix* a, double* s);
void rowcentering(matrix* a, matrix* b, double* c);
void colcentering(matrix* a, matrix* b, double* c);
void rowcentermeans(matrix* a, matrix* b, double* m);
void colcentermeans(matrix* a, matrix* b, double* m);
void rowweighting(matrix* a, matrix* b, double* w);
void colweighting(matrix* a, matrix* b, double* w);
void addmatrix(matrix* a, matrix* b, matrix* c);
void subtractmatrix(matrix* a, matrix* b, matrix* c);
void matrixscalar(matrix* a, double b, matrix* c);
void matrixdotproduct(matrix* a, matrix* b, matrix* c);
void matrixproduct(matrix* a, matrix* b, matrix* c);
void matrixweightedproduct(matrix* a, double*d, matrix* b, matrix* c);
void matrixtransproduct(matrix* a, matrix* b, matrix* c);
void matrixproducttrans(matrix* a, matrix* b, matrix* c);
void matrixproductweightedtrans(matrix* a, double* d, matrix* b, matrix* c);
void getdiagonal(matrix* mat, double* a);
void getrow(matrix* mat, unsigned int i, double* a);
void getcolumn(matrix* mat, unsigned int j, double* a);

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
