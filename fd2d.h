#ifndef FD2D_H
#define FD2D_H

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

    float *Vx, *Vz;
    float *Txx, *Txz, *Tzz;

    float *DVxx, *DVxz, *DVzx, *DVzz;
    float *DTxx, *DTxz, *DTzx, *DTzz;

    float *pml_ax, *pml_az, *pml_axh, *pml_azh;
    float *pml_bx, *pml_bz, *pml_bxh, *pml_bzh;
    float *pml_cx, *pml_cz, *pml_cxh, *pml_czh;

    float *pml_Qxx, *pml_Qxz, *pml_Qzx, *pml_Qzz;
    float *pml_Pxx, *pml_Pxz, *pml_Pzx, *pml_Pzz;
  public:
    fd2d (int, int, int, float, float, float, int, int, int nor);
    void initvar();
    void freevar();
    void getfdcoef();
    void getpmlcoef(float, float, float);
    void snapset (float, float, float, float, int, int, int, int);
    void initabc (int, float, float *, float *, float *);
    void forward(float *, float *, float *, float *, int, int *, int, int);
    void load_source (float *, float *, int, int, int);
    void update_stress (float *lam, float *mu, float dtdx, float dtdz);
    void update_stress_fs (float *lam, float *mu, float dtdx, float dtdz);
    void update_stress_pml (int, float *lam, float *mu);
    void update_velocity (float *rho, float dtdx, float dtdz);
    void update_velocity_fs (float *rho, float dtdx, float dtdz);
    void update_velocity_pml (int, float *rho);
    float vpmax(float *, float *, float *, int);
    float vsmin(float *, float *, int);
};
#endif
