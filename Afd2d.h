#ifndef AFD2D_H
#define AFD2D_H

#define IX(i,j) ((i) * nx + (j))
#define IV(i,j) ((i) + (j) * nz)
#define IU(i,j) ((i + nor) + (j + nor) * nzp)

#define KTM 10

#include "para.h"

class fd2d {
  private:
    int nx, ny, nz, nt;     // model size including abs boundary
    int nxp, nyp, nzp;      // model size including abs + padding
    int nxb, nyb, nzb;      // pml boundary size
    int mwf, mxp, mzp;
    int nor;
    int ifs, izf;           // free surface
    int iztop, nbtop;
    float dx, dy, dz, dt;

    int snap_nx, snap_ny, snap_nz, snap_nt;    // save size
    float snap_dx, snap_dy, snap_dz, snap_dt;
    float snap_fx, snap_fy, snap_fz;

    float coef[20];

    float *p1, *p2, *p3;
    float *DPxx, *DPzz;
    float *dpdx, *dpdz, *dedx, *dedz;
    float *Expml, *Hxpml, *Ezpml, *Hzpml;
    float *axpml, *azpml;
  public:
    fd2d (int, int, int, float, float, float, int, int, int nor);
    void initvar();
    void freevar();
    void getfdcoef();
    void getpmlcoef(float *);
    void snapset (float, float, float, float, int, int, int, int);
    void initabc (int, float, float *);
    void forward(float *, float *, int, int, int *, int);
    void update (float *, int);
    void update_fs ();
    void update_pml (int, float *);
    //float vpmax(float *, float *, float *, int);
    //float vsmin(float *, float *, int);
};
#endif
