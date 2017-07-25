#ifndef FORWARD_H
#define FORWARD_H

void MinMaxVel (float *vmin, float *vmax, float *lam, float *mu, float *rho, int n); 

void StackSource (float *u, float *wv, float p, float xo, float dx, float dt, int nx, int nt);


void MakeDirect (float *ud, float *lam, float *mu, float *rho, Shot src, Fwdpar par, Geometry geo, int, int, int);

void forward (float *uo, float *lam, float *mu, float *rho, Shot src, Fwdpar par, Geometry geo, int, int, int);

#endif
