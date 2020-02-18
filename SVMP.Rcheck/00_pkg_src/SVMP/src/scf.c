/*************************************************************************
 
 (c) 2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Spatial covariance functions**
 
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
 
 C function definitions
 
 *************************************************************************/

#include"scf.h"

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
