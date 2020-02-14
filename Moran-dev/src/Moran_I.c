/*
 * From Guillaume Guénard - Université de Montréal
 * August 2019
 * License: GPL version 3
 * C code
 */

#include"Moran_I.h"

void moran(double* x, double* w, int* n, double* I,
           int* nperm, int* pv, int* retIp, double* Ip_ret) {
  double meanx = 0.0, ssx = 0.0, W = 0.0, vp, cp, Ip, absI, mabsI, tmp;
  double EI = -1.0/(*n - 1);
  int i, j, k, p, *pidx, intmod;
  for(i = 0; i < *n; i++)
    meanx += x[i];
  meanx /= *n;
  for(i = 0; i < *n; i++) {
    x[i] -= meanx;
    ssx += x[i]*x[i];
  }
  i = 0;
  j = 0;
  vp = 0.0;
  for(k = 0; k < *n*(*n); k++) {
    if(i!=j) {
      W += w[k];
      vp += w[k]*x[i]*x[j];
    }
    j++;
    if(j == *n) {
      i++;
      j = 0;
    }
  }
  cp = (*n)/(ssx*W);
  *I = vp*cp - EI;
  if(*nperm) {
    intmod = INTNUM / (*n);
    intmod = intmod ? intmod : 1;
    if(*I > 0.0) {
      absI = *I;
      mabsI = -(*I);
    } else {
      absI = -(*I);
      mabsI = *I;
    }
    pv[2] = 1;
    pidx = (int *)Calloc(*n, int);
    if(pidx == NULL)
      error("Cannot allocate memory");
    for(i = 0; i < *n; i++)
      pidx[i] = i;
    GetRNGstate();
    for(p = 0; p < *nperm; p++) {
      permute_idx(pidx, *n);
      i = 0;
      j = 0;
      vp = 0.0;
      for(k = 0; k < *n*(*n); k++) {
        if(i!=j) {
          vp += w[k]*x[pidx[i]]*x[pidx[j]];
        }
        j++;
        if(j == *n) {
          i++;
          j = 0;
        }
      }
      Ip = vp*cp - EI;
      if(*retIp)
        Ip_ret[p] = Ip;
      if(Ip <= mabsI)
        pv[0]++;
      else if(Ip >= absI)
        pv[2]++;
      else
        pv[1]++;
      if(!(p%intmod))
        R_CheckUserInterrupt();
    }
    PutRNGstate();
    Free(pidx);
  }
  return;
}

// Caller must initialize the rng before calling this function.
void permute_idx(int* idx, int n) {
  int i, ii, tmp;
  double rn;
  for(i = 0; i < n; i++) {
    do {
      rn = unif_rand();
    } while (rn == 1.0);
    ii = (int)(rn * n);
    tmp = idx[i];
    idx[i] = idx[ii];
    idx[ii] = tmp;
  }
  return;
}
