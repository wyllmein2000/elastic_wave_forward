#ifndef ACOUSTIC_H
#define ACOUSTIC_H

void MinMaxVel (float *vmin, float *vmax, float *v, int n); 

void StackSource (float *u, float *wv, float p, float xo, float dx, float dt, int nx, int nt);

void MakeDirect (float *ud, float *vp, float *rho, Shot src, Fwdpar par, Geometry geo, int, int, int);

void forward (float *uo, float *vp, float *rho, Shot src, Fwdpar par, Geometry geo, int, int, int);

#endif
