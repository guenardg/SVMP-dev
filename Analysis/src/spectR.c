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

    C fonctions definitions

*************************************************************************/

#include"spectR.h"

/* Simple structure and functions for linear algebra (all double precision
 * floating points (IEEE 754-2008) */

/* Initialize a matrix structure and allocate necessary ressources.
 * Arguments
 * id   : identifier of the matrix
 * nr   : number of row on the matrix
 * nc   : number of columns of the matrix */
matrix initmatrix(char* id, unsigned int nr, unsigned int nc) {
  matrix mat;
  mat.id = id;
  mat.nr = nr;
  mat.nc = nc;
  mat.v = (double*)Calloc(nr*nc,double);
  if(mat.v == NULL)
    error("Unable to allocate ressources for matrix %s",mat.id);
  return mat;
}

/* Create a matrix structure and assign it existing data.
 * Arguments
 * id   : identifier of the matrix
 * nr   : number of row on the matrix
 * nc   : number of columns of the matrix
 * v    : pointer to an existing array (to be treated as a matrix) */
matrix assignmatrix(char* id, unsigned int nr, unsigned int nc, double* v) {
  matrix mat;
  mat.id = id;
  mat.nr = nr;
  mat.nc = nc;
  mat.v = v;
  return mat;
}

/* Clear a matrix that has been assigned in the calling code and free its data
 * pointer. Warning: to be used with matrices created using initmatrix(). */
void freematrix(matrix* mat) {
  Free(mat->v);
  if(mat->v != NULL)
    warning("Data from matrix %s could not be freed from memory.",mat->id);
  else {
    mat->nr = 0;
    mat->nc = 0;
  }
  return;
}

/* Clear a matrix and de-assign its data pointer. Warning: to be used for
 * matrices created with assignmatrix(). */
void deassignmatrix(matrix* mat) {
  mat->v = NULL;
  mat->nr = 0;
  mat->nc = 0;
  return;
}

// Make an element-wise copy of a matrix.
matrix copymatrix(matrix* a) {
  unsigned int i, n = a->nr*a->nc;
  matrix b = initmatrix(a->id,a->nr,a->nc);
  for(i = 0; i < n; i++)
    b.v[i] = a->v[i];
  return b;
}

// Sums of rows of matrix 'a' stored in array 's'
void rowsums(matrix* a, double* s) {
  unsigned int i, j, offset;
  double acc;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    acc = 0.0;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      acc += a->v[offset];
    s[i] = acc;
  }
  return;
}

// Sums of columns of matrix 'a' stored in array 's'
void colsums(matrix* a, double* s) {
  unsigned int i, j, offset;
  double acc;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr) {
    acc = 0.0;
    for (i = 0; i < a->nr; i++)
      acc += a->v[i + offset];
    s[j] = acc;
  }
  return;
}

/* Center the rows of matrix 'a' on values in array 'c' and write the result in 
 * matrix 'b'. To overwrite 'a', simply pass the same pointer for both arguments
 * 'a' and 'b'. */
void rowcentering(matrix* a, matrix* b, double* c) {
  unsigned int i, j, offset;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] - c[i];
  }
  return;
}

/* Center the columns of matrix 'a' on values in array 'c' and write the result
 * in matrix 'b'. To overwrite 'a', simply pass the same pointer for both
 * arguments 'a' and 'b'. */
void colcentering(matrix* a, matrix* b, double* c) {
  unsigned int i, j, offset;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr)
    for (i = 0; i < a->nr; i++)
      b->v[i+offset] = a->v[i+offset] - c[j];
  return;
}

/* Calculate the row means of matrix 'a' and store then as array 'm', then
 * center the rows of matrix 'a' on values in array 'c' and write the result in 
 * matrix 'b'. To overwrite 'a', simply pass the same pointer for both arguments
 * 'a' and 'b'. */
void rowcentermeans(matrix* a, matrix* b, double* m) {
  unsigned int i, j, offset;
  double acc;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    acc = 0.0;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      acc += a->v[offset];
    m[i] = acc/(a->nc);
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] - m[i];
  }
  return;
}

/* Calculate the column means of matrix 'a' and store then as array 'm', then
 * center the columns of matrix 'a' on values in array 'c' and write the result
 * in matrix 'b'. To overwrite 'a', simply pass the same pointer for both
 * arguments 'a' and 'b'. */
void colcentermeans(matrix* a, matrix* b, double* m) {
  unsigned int i, j, offset;
  double acc;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr) {
    acc = 0.0;
    for (i = 0; i < a->nr; i++)
      acc += a->v[i + offset];
    m[j] = acc/(a->nr);
    for (i = 0; i < a->nr; i++)
      b->v[i+offset] = a->v[i+offset] - m[j];
  }
  return;
}

/* Weight the rows of matrix 'a' with values in array 'w' and write the result
 * in matrix 'b'. To overwrite 'a', simply pass the same pointer for both
 * arguments 'a' and 'b'. */
void rowweighting(matrix* a, matrix* b, double* w) {
  unsigned int i, j, offset;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] * w[i];
  }
  return;
}

/* Weight the columns of matrix 'a' with values in array 'w' and write the
 * result in matrix 'b'. To overwrite 'a', simply pass the same pointer for both
 * arguments 'a' and 'b'. */
void colweighting(matrix* a, matrix* b, double* w) {
  unsigned int i, j, offset;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr)
    for (i = 0; i < a->nr; i++)
      b->v[i+offset] = a->v[i+offset] * w[j];
  return;
}

/* Add matrix 'a' to matrix 'b' and send the result in matrix 'c'. */
void addmatrix(matrix* a, matrix* b, matrix* c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] + b->v[i];
  return;
}

/* Subtract matrix 'b' from matrix 'a' and send the result in matrix 'c'. */
void subtractmatrix(matrix* a, matrix* b, matrix* c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] - b->v[i];
  return;
}

/* Multiply matrix 'a' by a scalar 'b' and send the result in matrix 'c'. */
void matrixscalar(matrix* a, double b, matrix* c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] * b;
  return;
}

/* Dot multiply (Haddamar product) matrices 'a' and 'b', and store the result
 * in matrix 'c' */
void matrixdotproduct(matrix* a, matrix* b, matrix* c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] * b->v[i];
  return;
}

/* Multiply (matrix product) matrices 'a' and 'b', and store the result in
 * matrix 'c' */
void matrixproduct(matrix* a, matrix* b, matrix* c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    offset2 = 0;
    for (j = 0; j < c->nc; j++, offset2 += b->nr, offset1 += c->nr) {
      offset3 = 0;
      acc = 0.0;
      for (k = 0; k < a->nc; k++, offset3 += a->nc)
        acc += a->v[i+offset3] * b->v[k+offset2];
      c->v[i+offset1] = acc;
    }
  }
  return;
}

/* Calculates C = ADB, where:
 * C is a n x p matrix passed as argument 'c',
 * A is a n x m matrix passed as argument 'a',
 * D is a m x m square matrix whose main diagonal values are passed as argument
 * 'd' and whose off-diagnal values are 0 (i.e., a diagonal matrix), and
 * B is a m x p matrix passed as argument 'b'. */
void matrixweightedproduct(matrix* a, double* d, matrix* b, matrix* c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    offset2 = 0;
    for (j = 0; j < c->nc; j++, offset2 += b->nr, offset1 += c->nr) {
      offset3 = 0;
      acc = 0.0;
      for (k = 0; k < a->nc; k++, offset3 += a->nc)
        acc += a->v[i+offset3] * d[k] * b->v[k+offset2];
      c->v[i+offset1] = acc;
    }
  }
  return;
}

/* Calculates C = A'B, where:
 * C is a n x p matrix passed as argument 'c',
 * A is a m x n matrix passed as argument 'a', and
 * B is a m x p matrix passed as argument 'b'. */
void matrixtransproduct(matrix* a, matrix* b, matrix* c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  offset1 = 0;
  for (i = 0; i < c->nr; i++, offset1 += a->nr) {
    offset2 = 0;
    offset3 = 0;
    for (j = 0; j < c->nc; j++, offset3 += b->nr, offset2 += c->nr) {
      acc = 0.0;
      for (k = 0; k < a->nr; k++)
        acc += a->v[k+offset1] * b->v[k+offset3];
      c->v[i+offset2] = acc;
    }
  }
  return;
}

/* Calculates C = AB', where:
 * C is a n x p matrix passed as argument 'c',
 * A is a n x m matrix passed as argument 'a', and
 * B is a p x m matrix passed as argument 'b'. */
void matrixproducttrans(matrix* a, matrix* b, matrix* c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    for (j = 0; j < c->nc; j++, offset1 += c->nr) {
      acc = 0.0;
      offset2 = 0;
      offset3 = 0;
      for(k = 0; k < a->nc; k++, offset2 += a->nr, offset3 += b->nr)
        acc += a->v[i+offset2] * b->v[j+offset3];
      c->v[i+offset1] = acc;
    }
  }
  return;
}

/* Calculates C = ADB', where:
 * C is a n x p matrix passed as argument 'c',
 * A is a n x m matrix passed as argument 'a',
 * D is a m x m square matrix whose main diagonal values are passed as argument
 * 'd' and whose off-diagnal values are 0 (i.e., a diagonal matrix), and
 * B is a p x m matrix passed as argument 'b'. */
void matrixproductweightedtrans(matrix* a, double* d, matrix* b, matrix* c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    for (j = 0; j < c->nc; j++, offset1 += c->nr) {
      acc = 0.0;
      offset2 = 0;
      offset3 = 0;
      for(k = 0; k < a->nc; k++, offset2 += a->nr, offset3 += b->nr)
        acc += a->v[i+offset2] * d[k] * b->v[j+offset3];
      c->v[i+offset1] = acc;
    }
  }
  return;
}

/* Extract the main diagonal of matrix 'a' as store it in array 'd'. */
void getdiagonal(matrix* a, double* d) {
  unsigned int i, order, offset = 0;
  order = (a->nr < a->nc)?a->nr:a->nc;
  for (i = 0; i < order; i++, offset += a->nr)
    d[i] = a->v[i+offset];
  return;
}

/* Extract row 'i' of matrix 'a' and store it in array 'r'
 * ATTENTION: indices begin with 0 and < a->nr) */
void getrow(matrix* a, unsigned int i, double* r) {
  unsigned int j, offset = i;
  for (j = 0; j < a->nc; j++, offset += a->nr)
    r[j] = a->v[offset];
  return;
}

/* Extract column 'j' of matrix 'a' and store it in array 'c'
 * ATTENTION: indices begin with 0 and < a->nc) */
void getcolumn(matrix* a, unsigned int j, double* c) {
  unsigned int i = 0, offset = a->nr * i;
  for (; i < a->nc; i++, offset++)
    c[i] = a->v[offset];
  return;
}

// MLE-friendly MEM weighting functions

void scf_spher(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res) {
  int i;
  double a, am, ep, dd;
  if(*recycle) {
    if((*alpha) == 0.5) {
      for(i = 0; i < *n; i++) {
        if(d[i] == 0.0)
          res[i] = 1.0;
        else if(d[i] < (*beta)) {
          dd = d[i]/(*beta);
          res[i] = dd*log(dd) - dd + 1.0;
        } else
          res[i] = 0.0;
      }
    } else if((*alpha) > 0) {
      a = 0.5 - 0.25/(0.5 - (*alpha));
      am = a - 1.0;
      ep = a / am;
      for(i = 0; i < *n; i++) {
        if(d[i] < (*beta)) {
          dd = d[i]/(*beta);
          res[i] = 1.0 - a*dd + am*pow(dd,ep);
        } else
          res[i] = 0.0;
      }
    } else
      for(i = 0; i < *n; i++) {
        if(d[i] == 0.0)
          res[i] = 1.0;
        else
          res[i] = 0.0;
      }
  } else {
    for(i = 0; i < *n; i++)
      if(d[i] < beta[i]) {
        dd = d[i]/beta[i];
        if(alpha[i] == 0.5)
          res[i] = dd*log(dd) - dd + 1.0;
        else if(alpha[i] > 0.0) {
          a = 0.5 - 0.25/(0.5 - alpha[i]);
          am = a - 1.0;
          ep = a / am;
          res[i] = 1.0 - a*dd + am*pow(dd,ep);
        } else
          if(d[i] == 0.0)
            res[i] = 1.0;
          else
            res[i] = 0.0;
      } else
        res[i] = 0.0;
  }
  return;
}

void scf_expon(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res) {
  int i;
  double bi, a, os, sc;
  if(*recycle) {
    if(((*alpha)<1.0)&&((*alpha)>0.0)) {
      bi = 1.0/(*beta);
      a = ((*alpha) - 1.0)/(*alpha);
      os = exp(a);
      sc = 1.0/(1.0 - os);
      for(i = 0; i < *n; i++)
        if(d[i] < (*beta))
          res[i] = sc*(exp(a*bi*d[i]) - os);
        else
          res[i] = 0.0;
    } else if(!((*alpha)<1.0)) {
      bi = 1.0/(*beta);
      for(i = 0; i < *n; i++)
        if(d[i] < (*beta))
          res[i] = 1.0 - bi*d[i];
        else
          res[i] = 0.0;
    } else
        for(i = 0; i < *n; i++)
          res[i] = (d[i]==0.0)?1.0:0.0;
  } else
    for(i = 0; i < *n; i++)
      if(d[i] < beta[i]) {
        bi = 1.0/beta[i];
        if((alpha[i]<1.0)&&(alpha[i]>0.0)) {
          a = (alpha[i] - 1.0)/alpha[i];
          os = exp(a);
          sc = 1.0/(1.0 - os);
          res[i] = sc*(exp(a*bi*d[i]) - os);
        } else if(!((*alpha)<1.0))
          res[i] = 1.0 - bi*d[i];
        else
          res[i] = (d[i]==0.0)?1.0:0.0;
      } else
        res[i] = 0.0;
  return;
}

void scf_power(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res) {
  int i;
  double dd, aa;
  if(*recycle) {
    aa = 1.0/(*alpha);
    for(i = 0; i < *n; i++) {
      dd = d[i]/(*beta);
      res[i] = (dd < 1.0)?pow(1.0 - dd,aa):0.0;
    }
  } else {
    for(i = 0; i < *n; i++) {
      dd = d[i]/beta[i];
      res[i] = (dd < 1.0)?pow(1.0 - dd,1.0/alpha[i]):0.0;
    }
  }
  return;
}

void scf_hyper(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res) {
  int i;
  double a, os, sc;
  if(*recycle) {
    if((*alpha) == 0.5) {
      a = 1.0/log((*beta) + 1.0);
      for(i = 0; i < *n; i++)
        if(d[i] < (*beta))
          res[i] = 1.0 - a*log(d[i] + 1.0);
        else
          res[i] = 0.0;
    } else {
      // a = 1.0/(0.5 - 0.25/(0.5 - (*alpha)));  // simplifies to:
      a = 2.0 - 1.0/(*alpha);
      os = pow((*beta) + 1.0,a);
      sc = 1.0/(1.0 - os);
      for(i = 0; i < *n; i++)
        if(d[i] < (*beta))
          res[i] = sc*(pow(d[i] + 1.0,a) - os);
        else
          res[i] = 0.0;
    }
  } else {
    for(i = 0; i < *n; i++) {
      if(alpha[i] == 0.5) {
        if(d[i] < beta[i]) {
          a = 1.0/log(beta[i] + 1.0);
          res[i] = 1.0 - a*log(d[i] + 1.0);
        } else {
          res[i] = 0.0;
        }
      } else {
        if(d[i] < beta[i]) {
          // a = 1.0/(0.5 - 0.25/(0.5 - alpha[i])); // simplifies to:
          a = 2.0 - 1.0/alpha[i];
          os = pow(beta[i] + 1.0,a);
          sc = 1.0/(1.0 - os);
          res[i] = sc*(pow(d[i] + 1.0,a) - os);
        } else
          res[i] = 0.0;
      }
    }
  }
  return;
}

void scf_super(double* d, double* alpha, double* beta,
               int* n, int* recycle, double* res) {
  int i;
  double ep, de, di;
  if(*recycle) {
    if((*alpha) > 0) {
      ep = 1.0/(*alpha);
      de = pow(*beta,ep);
      di = 1.0/(*beta);
      for(i = 0; i < *n; i++)
        if(d[i] < (*beta))
          res[i] = 1.0 - di*pow(de - pow((*beta) - d[i],ep),*alpha);
        else
          res[i] = 0.0;
    } else
      for(i = 0; i < *n; i++)
        if(d[i] == 0.0)
          res[i] = 1.0;
        else
          res[i] = 0.0;
  } else {
    for(i = 0; i < *n; i++)
      if(d[i] < beta[i]) {
        if(alpha[i] > 0.0) {
          ep = 1.0/alpha[i];
          de = pow(beta[i],ep);
          res[i] = 1.0 - pow(de - pow(beta[i] - d[i],ep),alpha[i])/beta[i];
        } else
          if(d[i] == 0.0)
            res[i] = 1.0;
          else
            res[i] = 0.0;
      } else
        res[i] = 0.0;
  }
  return;
}

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
