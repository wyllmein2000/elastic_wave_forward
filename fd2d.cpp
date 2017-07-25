#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "math.h"
#include "numpy.h"
#include "inout.h"
#include "fd2d.h"
using namespace std;

/******************************************
 *     input size: nx = nxi + 2 * nxpml
 * wavefield size: nxp = nx + 2 * nor;
 *  pml w.f. size: nxpml
 * ***************************************/
fd2d::fd2d (int nx, int nz, int nt, float dx, float dz, float dt, int nxb, int nzb, int order) {
   this->nor = order / 2;              // computation grid
   this->nx = nx;
   this->nz = nz;
   this->nt = nt;
   this->dx = dx;
   this->dz = dz;
   this->dt = dt;
   this->nxb = nxb;
   this->nzb = nzb;
   this->nxp = nx + 2 * nor;
   this->nzp = nz + 2 * nor;

   this->snap_fx = nor * dx;           // snapshot
   this->snap_fy = nor * dy;
   this->snap_fz = nor * dz;
   this->snap_dx = dx;
   this->snap_dy = dy;
   this->snap_dz = dz;
   this->snap_dt = dt;
   this->snap_nx = nx;
   this->snap_ny = ny;
   this->snap_nz = nz;
   this->snap_nt = nt;

   this->mwf = nxp * nzp;             // total wavefield size
   this->mxp = 2 * nxb * nz;          // pml wavefield size
   this->mzp = 2 * nx * nzb;

   this->Vx = new float[mwf];
   this->Vz = new float[mwf];
   this->Txx = new float[mwf];
   this->Txz = new float[mwf];
   this->Tzz = new float[mwf];

   this->DVxx = new float[mxp];
   this->DVzx = new float[mxp];
   this->DVxz = new float[mzp];
   this->DVzz = new float[mzp];
   this->DTxx = new float[mxp];
   this->DTzx = new float[mxp];
   this->DTxz = new float[mzp];
   this->DTzz = new float[mzp];

   this->pml_ax = new float[nxb];
   this->pml_bx = new float[nxb];
   this->pml_cx = new float[nxb];
   this->pml_axh = new float[nxb];
   this->pml_bxh = new float[nxb];
   this->pml_cxh = new float[nxb];
   this->pml_az = new float[nzb];
   this->pml_bz = new float[nzb];
   this->pml_cz = new float[nzb];
   this->pml_azh = new float[nzb];
   this->pml_bzh = new float[nzb];
   this->pml_czh = new float[nzb];

   this->pml_Qxx = new float[mxp];
   this->pml_Qxz = new float[mzp];
   this->pml_Qzx = new float[mxp];
   this->pml_Qzz = new float[mzp];
   this->pml_Pxx = new float[mxp];
   this->pml_Pxz = new float[mzp];
   this->pml_Pzx = new float[mxp];
   this->pml_Pzz = new float[mzp];
}


void fd2d::initvar() {
   memset(Vx, 0, mwf * sizeof(float));
   memset(Vz, 0, mwf * sizeof(float));
   memset(Txx, 0, mwf * sizeof(float));
   memset(Txz, 0, mwf * sizeof(float));
   memset(Tzz, 0, mwf * sizeof(float));

   memset(DVxx, 0, mxp * sizeof(float));
   memset(DVzx, 0, mxp * sizeof(float));
   memset(DTxx, 0, mxp * sizeof(float));
   memset(DTzx, 0, mxp * sizeof(float));
   memset(DVxz, 0, mzp * sizeof(float));
   memset(DVzz, 0, mzp * sizeof(float));
   memset(DTxz, 0, mzp * sizeof(float));
   memset(DTzz, 0, mzp * sizeof(float));

   memset(pml_Qxx, 0, mxp * sizeof(float));
   memset(pml_Qzx, 0, mxp * sizeof(float));
   memset(pml_Pxx, 0, mxp * sizeof(float));
   memset(pml_Pzx, 0, mxp * sizeof(float));
   memset(pml_Qxz, 0, mzp * sizeof(float));
   memset(pml_Qzz, 0, mzp * sizeof(float));
   memset(pml_Pxz, 0, mzp * sizeof(float));
   memset(pml_Pzz, 0, mzp * sizeof(float));
}

void fd2d::initabc (int iz_free_surface, float fp, float *lam, float *mu, float *rho) {
  float v1, v2;
  v1 = vpmax(lam, mu, rho, nx * nz);
  v2 = vsmin(mu, rho, nx * nz);
  //cout << " vpmax=" << v1 << ", vsmin=" << v2 << endl;
  getpmlcoef (v1, v2, fp);     // pml

  this->ifs = 0;
  this->izf = iz_free_surface;
  this->iztop = 0;
  this->nbtop = nzb;

  if (this->izf >= 0) {              // free surface
      this->ifs == 1;
      this->iztop = izf + 1;
      this->nbtop = 0;
  }
}

void fd2d::snapset (float dxo, float dyo, float dzo, float dto, int nxo, int nyo, int nzo, int nto) {
  this->snap_fx = (nxb + nor) * dx;
  this->snap_fy = (nyb + nor) * dy;
  this->snap_fz = (nzb + nor) * dz;
  this->snap_dx = dxo;
  this->snap_dy = dyo;
  this->snap_dz = dzo;
  this->snap_dt = dto;
  this->snap_nx = nxo;
  this->snap_ny = nyo;
  this->snap_nz = nzo;
  this->snap_nt = nto;
}

void fd2d::forward(float *u, float *lam, float *mu, float *rho, int izs, int *izr, int isrc, int iw) {
  int jx, jz;
  int ktm = snap_dt / dt;
  float dtdx = dt / dx;
  float dtdz = dt / dz;
  float dtxz = dt / (dx * dz);
  float *snapshot;

  getfdcoef ();
  initvar ();

  ofstream fsnap1, fsnap2;
  if (iw == 1) {
      snapshot = new float[snap_nx * snap_nz];
      fsnap1.open("temp/txx.bin", ios::out|ios::binary);
      fsnap2.open("temp/tzz.bin", ios::out|ios::binary);
      cout << " model size: nx=" << nx << " nz=" << nz << endl;
      cout << "             dx=" << dx << " dz=" << dz << " dt=" << dt << endl;
      cout << " save u(t,z,x): nx=" << snap_nx << " nz=" << snap_nz << " nt=" << nt/ktm << endl;
      cout << " FD COEF: c1=" << coef[1] << " c2=" << coef[2] << " c3=" << coef[3] << endl;
  }


  for (int it = 0; it < nt; it ++) {

       if (iw == 1 && it % 500 == 0)
           cout << " *** it = " << it << "/" << nt << endl;

       update_stress(lam, mu, dtdx, dtdz);
       update_stress_pml (it, lam, mu);
       if (ifs == 1) 
           update_stress_fs (lam, mu, dtdx, dtdz);
       update_velocity(rho, dtdx, dtdz);
       update_velocity_pml (it, rho);
       if (ifs == 1) 
           update_velocity_fs (rho, dtdx, dtdz);
       load_source (u, rho, it, izs, isrc);

       // output trace
       for (int ix = 0; ix < nx; ix ++) {
            jz = izr[ix];            
            //u[ix*nt+it] = Vz[IU(jz,ix)];
            u[ix*nt+it] = 0.5 * (Vz[IU(jz,ix)] + Vz[IU(jz-1,ix)]);
       }

       // output snapshot
       if (iw == 1 && it % ktm == 0) {
           interp2d (snapshot, snap_fx, snap_fz, snap_dx, snap_dz,
                     snap_nx, snap_nz, Txx, 0.0, 0.0, dx, dz, nxp, nzp);
           fsnap1.write((char *)snapshot, snap_nx * snap_nz * sizeof(float));
           interp2d (snapshot, snap_fx, snap_fz, snap_dx, snap_dz,
                     snap_nx, snap_nz, Tzz, 0.0, 0.0, dx, dz, nxp, nzp);
           fsnap2.write((char *)snapshot, snap_nx * snap_nz * sizeof(float));
       }
  }

  if (iw == 1) {
      fsnap1.close();
      fsnap2.close();
  }

  freevar ();
}


void fd2d::load_source (float *u, float *rho, int it, int izs, int isrc) {
  if (isrc == 0) 
      for (int ix = 0; ix < nx; ix ++) {
           Txx[IU(izs,ix)] += u[ix*nt+it] * dt;
           Tzz[IU(izs,ix)] += u[ix*nt+it] * dt;
      }
  else if (isrc == 1)
      for (int ix = 0; ix < nx; ix ++) {
           Vz[IU(izs,  ix)] += u[ix*nt+it]/(rho[IV(izs,ix)] + rho[IV(izs+1,ix)]) * dt;
           Vz[IU(izs-1,ix)] += u[ix*nt+it]/(rho[IV(izs,ix)] + rho[IV(izs-1,ix)]) * dt;
      }
}

/****************************************************
 * (Txx,Tzz)
 *     0           1
 *   0 *-----Vx----*
 *     |     |     |
 *     |     |     |
 *     Vz---Txz----Vz
 *     |     |     |
 *     |     |     |
 *   1 *-----Vx----*
 *     |
 **************************************************** */
void fd2d::update_stress(float *lam, float *mu, float dtdx, float dtdz) {

  //cout << "Txx: dtdx=" << dtdx << ", dtdz=" << dtdz << endl;

  #pragma omp parallel for
  for (int j = 0; j < nx; j ++)
  for (int i = iztop; i < nz; i ++) {

       int iv1, iv2, iv3, iv4;
       float C11, C13, C55;

       float vx_x = 0;   // vx_x(i,j)
       float vz_z = 0;   // vz_z(i,j)
       float vx_z = 0;   // vx_z(i+0.5,j+0.5)
       float vz_x = 0;   // vz_x(i+0.5,j+0.5)

       for (int k = 1; k <= nor; k ++) {
            vx_x += coef[k] * (Vx[IU(i,j+k-1)] - Vx[IU(i,j-k)]);
            vx_z += coef[k] * (Vx[IU(i+k,j)] - Vx[IU(i-k+1,j)]);
            vz_x += coef[k] * (Vz[IU(i,j+k)] - Vz[IU(i,j-k+1)]);
            vz_z += coef[k] * (Vz[IU(i+k-1,j)] - Vz[IU(i-k,j)]);
       }


       iv1 = IV(i, j);
       iv2 = IV(i+1, j);
       iv3 = IV(i, j+1);
       iv4 = IV(i+1, j+1);

       C11 = lam[iv1] + 2.0 * mu[iv1];
       C13 = lam[iv1];

       Txx[IU(i,j)] += C11 * dtdx * vx_x + C13 * dtdz * vz_z;
       Tzz[IU(i,j)] += C13 * dtdx * vx_x + C11 * dtdz * vz_z;

       if (j < nx - 1 && i < nz - 1) {
           C55 = 0.25 * (mu[iv1] + mu[iv2] + mu[iv3] + mu[iv4]);
           Txz[IU(i,j)] += C55 * (dtdx * vz_x + dtdz * vx_z);
       }

       if (j < nxb) {
           DVxx[iv1] = vx_x / dx;
           DVzx[iv1] = vz_x / dx;
       }
       if (j >= nx - nxb)
          DVxx[IV(i, j - (nx - 2 * nxb))] = vx_x / dx;
       if (j >= nx - 1 - nxb && j < nx - 1)
          DVzx[IV(i, j - (nx - 1 - 2 * nxb))] = vz_x / dx;

       if (i < nbtop) {
          DVxz[IX(i, j)] = vx_z / dz;
          DVzz[IX(i, j)] = vz_z / dz;
       }
       if (i >= nz - nzb)
          DVzz[IX(i - (nz - 2 * nzb), j)] = vz_z / dz;
       if (i >= nz - 1 - nzb && i < nz - 1)
          DVxz[IX(i - (nz - 1 - 2 * nzb), j)] = vx_z / dz;
  }
}

void fd2d::update_velocity (float *rho, float dtdx, float dtdz) {

  #pragma omp parallel for
  for (int j = 0; j < nx; j ++)
  for (int i = iztop; i < nz; i ++) {

       float rhox, rhoz;
       int iv = IV(i, j);
       int iv1 = IV(i+1, j);
       int iv2 = IV(i, j+1);

       float txx_x = 0;   // txx_x(i,j+0.5);
       float txz_z = 0;   // txz_z(i,j+0.5);
       float txz_x = 0;   // txz_z(i+0.5,j);
       float tzz_z = 0;   // txz_z(i+0.5,j);

       /*
       for (int k = 1; k <= nor; k ++) {
            txx_x += coef[k] * (Txx[IU(i,j+k)] - Txx[IU(i,j-k+1)]);
            txz_z += coef[k] * (Txz[IU(i+k-1,j)] - Txz[IU(i-k,j)]);
            txz_x += coef[k] * (Txz[IU(i,j+k-1)] - Txz[IU(i,j-k)]);
            tzz_z += coef[k] * (Tzz[IU(i+k,j)] - Tzz[IU(i-k+1,j)]);
       }*/

       if (j < nx - 1) {
           rhox = 0.5 * (rho[iv] + rho[iv2]);
           for (int k = 1; k <= nor; k ++) {
                txx_x += coef[k] * (Txx[IU(i,j+k)] - Txx[IU(i,j-k+1)]);
                txz_z += coef[k] * (Txz[IU(i+k-1,j)] - Txz[IU(i-k,j)]);
           }
           Vx[IU(i,j)] += (dtdx * txx_x + dtdz * txz_z) / rhox;
       }

       if (i < nz - 1) {
           rhoz = 0.5 * (rho[iv] + rho[iv1]);
           for (int k = 1; k <= nor; k ++) {
                txz_x += coef[k] * (Txz[IU(i,j+k-1)] - Txz[IU(i,j-k)]);
                tzz_z += coef[k] * (Tzz[IU(i+k,j)] - Tzz[IU(i-k+1,j)]);
           }
           Vz[IU(i,j)] += (dtdx * txz_x + dtdz * tzz_z) / rhoz;
       }

       if (j < nxb) {
          DTxx[iv] = txx_x / dx;
          DTzx[iv] = txz_x / dx;
       }
       if (j >= nx - 1 - nxb && j < nx - 1)
          DTxx[IV(i, j - (nx - 1 - 2 * nxb))] = txx_x / dx;
       if (j >= nx - nxb)
          DTzx[IV(i, j - (nx - 2 * nxb))] = txz_x / dx;

       if (i < nbtop) {
          DTxz[IX(i, j)] = txz_z / dz;
          DTzz[IX(i, j)] = tzz_z / dz;
       }
       if (i >= nz - 1 - nzb && i < nz - 1)
          DTzz[IX(i - (nz - 1 - 2 * nzb), j)] = tzz_z / dz;
       if (i >= nz - nzb)
          DTxz[IX(i - (nz - 2 * nzb), j)] = txz_z / dz;
  }
}

void fd2d::update_stress_fs (float *lam, float *mu, float dtdx, float dtdz) {

  //cout << "Txx: dtdx=" << dtdx << ", dtdz=" << dtdz << endl;

  #pragma omp parallel for
  for (int j = 0; j < nx; j ++) {
       int i, k;
       int iv1, iv2, iv3, iv4;
       float C11, C13, C33, C55;

       float vx_x = 0;   // vx_x(i,j)
       float vz_z = 0;   // vz_z(i,j)
       float vx_z = 0;   // vx_z(i+0.5,j+0.5)
       float vz_x = 0;   // vz_x(i+0.5,j+0.5)

       // virtual region
       for (i = 0; i < izf; i ++) {
            k = 2 * izf - i;
            //Txx[IU(i,j)] is not used
            Tzz[IU(i,j)] = Tzz[IU(k,j)];
            Txz[IU(i,j)] = Txz[IU(k-1,j)];
       }          
  
       // free surface
       i = izf;

       for (int k = 1; k <= nor; k ++) {
            vx_x += coef[k] * (Vx[IU(i,j+k-1)] - Vx[IU(i,j-k)]);
            vx_z += coef[k] * (Vx[IU(i+k,j)] - Vx[IU(i-k+1,j)]);
            vz_x += coef[k] * (Vz[IU(i,j+k)] - Vz[IU(i,j-k+1)]);
            vz_z += coef[k] * (Vz[IU(i+k-1,j)] - Vz[IU(i-k,j)]);
       }

       iv1 = IV(i, j);
       iv2 = IV(i+1, j);
       iv3 = IV(i, j+1);
       iv4 = IV(i+1, j+1);

       C11 = lam[iv1] + 2.0 * mu[iv1];
       C13 = lam[iv1];
       C33 = C11;

       Txx[IU(i,j)] += (C11 - C13 * C13 / C33) * dtdx * vx_x;
       Tzz[IU(i,j)] = 0;

       if (j < nx - 1 && i < nz - 1) {
           C55 = 0.25 * (mu[iv1] + mu[iv2] + mu[iv3] + mu[iv4]);
           Txz[IU(i,j)] += C55 * (dtdx * vz_x + dtdz * vx_z);
       }

       if (j < nxb) {
           DVxx[iv1] = vx_x / dx;
           DVzx[iv1] = vz_x / dx;
       }
       if (j >= nx - nxb)
          DVxx[IV(i, j - (nx - 2 * nxb))] = vx_x / dx;
       if (j >= nx - 1 - nxb && j < nx - 1)
          DVzx[IV(i, j - (nx - 1 - 2 * nxb))] = vz_x / dx;

       if (i < nzb) {
          DVxz[IX(i, j)] = vx_z / dz;
          DVzz[IX(i, j)] = vz_z / dz;
       }
       if (i >= nz - nzb)
          DVzz[IX(i - (nz - 2 * nzb), j)] = vz_z / dz;
       if (i >= nz - 1 - nzb && i < nz - 1)
          DVxz[IX(i - (nz - 1 - 2 * nzb), j)] = vx_z / dz;
  }
}

void fd2d::update_velocity_fs (float *rho, float dtdx, float dtdz) {

  #pragma omp parallel for
  for (int j = 0; j < nx; j ++) {

       int i, k;
       float rhox, rhoz;

       float txx_x = 0;   // txx_x(i,j+0.5);
       float txz_z = 0;   // txz_z(i,j+0.5);
       float txz_x = 0;   // txz_z(i+0.5,j);
       float tzz_z = 0;   // txz_z(i+0.5,j);

       // above the free-surface
       for (i = 0; i < izf; i ++) {
            Vx[IU(i,j)] = 0;
            Vz[IU(i,j)] = 0;
       }
 
       // on the free surface
       i = izf;

       int iv = IV(i, j);
       int iv1 = IV(i+1, j);
       int iv2 = IV(i, j+1);

       /*
       for (int k = 1; k <= nor; k ++) {
            txx_x += coef[k] * (Txx[IU(i,j+k)] - Txx[IU(i,j-k+1)]);
            txz_z += coef[k] * (Txz[IU(i+k-1,j)] - Txz[IU(i-k,j)]);
            txz_x += coef[k] * (Txz[IU(i,j+k-1)] - Txz[IU(i,j-k)]);
            tzz_z += coef[k] * (Tzz[IU(i+k,j)] - Tzz[IU(i-k+1,j)]);
       }*/

       if (j < nx - 1) {
           rhox = 0.5 * (rho[iv] + rho[iv2]);
           for (int k = 1; k <= nor; k ++) {
                txx_x += coef[k] * (Txx[IU(i,j+k)] - Txx[IU(i,j-k+1)]);
                //txz_z += coef[k] * (Txz[IU(i+k-1,j)] - Txz[IU(i-k,j)]);
           }
           Vx[IU(i,j)] += dtdx * txx_x / rhox;
       }

       if (i < nz - 1) {
           rhoz = 0.5 * (rho[iv] + rho[iv1]);
           for (int k = 1; k <= nor; k ++) {
                txz_x += coef[k] * (Txz[IU(i,j+k-1)] - Txz[IU(i,j-k)]);
                tzz_z += coef[k] * (Tzz[IU(i+k,j)] - Tzz[IU(i-k+1,j)]);
           }
           Vz[IU(i,j)] += (dtdx * txz_x + dtdz * tzz_z) / rhoz;
       }

       if (j < nxb) {
          DTxx[iv] = txx_x / dx;
          DTzx[iv] = txz_x / dx;
       }
       if (j >= nx - 1 - nxb && j < nx - 1)
          DTxx[IV(i, j - (nx - 1 - 2 * nxb))] = txx_x / dx;
       if (j >= nx - nxb)
          DTzx[IV(i, j - (nx - 2 * nxb))] = txz_x / dx;

       if (i < nzb) {
          DTxz[IX(i, j)] = txz_z / dz;
          DTzz[IX(i, j)] = tzz_z / dz;
       }
       if (i >= nz - 1 - nzb && i < nz - 1)
          DTzz[IX(i - (nz - 1 - 2 * nzb), j)] = tzz_z / dz;
       if (i >= nz - nzb)
          DTxz[IX(i - (nz - 2 * nzb), j)] = txz_z / dz;
  }
}

/* Note 
 */
void fd2d::update_velocity_pml (int istep, float *rho) {

  float t, a, b;
  float th, ah, bh;

  // x - direction
  #pragma omp parallel for
  for (int iz = 0; iz < nz; iz ++)
  for (int ix = 0; ix < nxb; ix ++) {

       float c, ch, rhox, rhoz;

       // left side
       int iq = IV(iz, ix);
       int ip = IU(iz, ix);
       int ir = IV(iz, ix);
       int ir1 = IV(iz, ix + 1);
       int ir2 = IV(iz + 1, ix);

       c  = 1.0 - pml_cx [ix];
       ch = 1.0 - pml_cxh[ix];

       rhox = 2.0 * dt / (rho[ir] + rho[ir1]);
       Vx[ip] -= rhox * (ch * pml_Qxx[iq] + pml_cxh[ix] * DTxx[iq]);

       if (iz < nz - 1) {
           rhoz = 2.0 * dt / (rho[ir] + rho[ir2]);
           Vz[ip] -= rhoz * (c * pml_Qzx[iq] + pml_cx[ix] * DTzx[iq]);
       }
   
       pml_Qxx[iq] = pml_axh[ix] * pml_Qxx[iq] + pml_bxh[ix] * DTxx[iq];
       pml_Qzx[iq] = pml_ax [ix] * pml_Qzx[iq] + pml_bx [ix] * DTzx[iq];


       // right side
       int jx = nxb - 1 - ix;
       int jq = IV(iz, nxb + jx);

       int jx1 = nx - 2 - ix;
       int jx2 = nx - 1 - ix;
       int jp1 = IU(iz, jx1);
       int jp2 = IU(iz, jx2);

       int jr1 = IV(iz, jx1);
       int jr2 = IV(iz, jx2);
       int jr3 = IV(iz + 1, jx2);

       rhox = 2.0 * dt / (rho[jr1] + rho[jr2]);
       Vx[jp1] -= rhox * (ch * pml_Qxx[jq] + pml_cxh[ix] * DTxx[jq]);

       if (iz < nz - 1) {
           rhoz = 2.0 * dt / (rho[jr2] + rho[jr3]);
           Vz[jp2] -= rhoz * (c * pml_Qzx[jq] + pml_cx[ix] * DTzx[jq]);
       }

       pml_Qxx[jq] = pml_axh[ix] * pml_Qxx[jq] + pml_bxh[ix] * DTxx[jq];
       pml_Qzx[jq] = pml_ax [ix] * pml_Qzx[jq] + pml_bx [ix] * DTzx[jq];
  }//

  // z - direction
  #pragma omp parallel for
  for (int ix = 0; ix < nx; ix ++)
  for (int iz = 0; iz < nzb; iz ++) {

       // top side
       int iq = IX(iz, ix);
       int ip = IU(iz, ix);

       int ir = IV(iz, ix);
       int ir1 = IV(iz, ix + 1);
       int ir2 = IV(iz + 1, ix);
       
       //
       float c = 1.0 - pml_cz[iz];
       float ch = 1.0 - pml_cz[iz];

       if (ifs == 0) {

       if (ix < nx - 1)
       Vx[ip] -= (c  * pml_Qxz[iq] + pml_cz [iz] * DTxz[iq]) * dt * 2.0 / (rho[ir] + rho[ir1]);
       Vz[ip] -= (ch * pml_Qzz[iq] + pml_czh[iz] * DTzz[iq]) * dt * 2.0 / (rho[ir] + rho[ir2]);
   
       pml_Qxz[iq] = pml_az [iz] * pml_Qxz[iq] + pml_bz [iz] * DTxz[iq];
       pml_Qzz[iq] = pml_azh[iz] * pml_Qzz[iq] + pml_bzh[iz] * DTzz[iq];
       }
       //

       // bottom side
       int jz = nzb - 1 - iz;
       int jq = IX(nzb + jz, ix);

       int jz1 = nz - 1 - iz;
       int jz2 = nz - 2 - iz;
       int jp1 = IU(jz1, ix);
       int jp2 = IU(jz2, ix);

       int jr = IV(jz1, ix);
       int jr1 = IV(jz1, ix + 1);
       int jr2 = IV(jz2, ix);

       if (ix < nx - 1)
       Vx[jp1] -= (c  * pml_Qxz[jq] + pml_cz [iz] * DTxz[jq]) * dt * 2.0 / (rho[jr] + rho[jr1]);
       Vz[jp2] -= (ch * pml_Qzz[jq] + pml_czh[iz] * DTzz[jq]) * dt * 2.0 / (rho[jr] + rho[jr2]);

       pml_Qxz[jq] = pml_az [iz] * pml_Qxz[jq] + pml_bz [iz] * DTxz[jq];
       pml_Qzz[jq] = pml_azh[iz] * pml_Qzz[jq] + pml_bzh[iz] * DTzz[jq];
    
   }
}

void fd2d::update_stress_pml (int istep, float *lam, float *mu) {

  // x - direction
  #pragma omp parallel for
  for (int iz = 0; iz < nz; iz ++)
  for (int ix = 0; ix < nxb; ix ++) {

       float C11, C13, C55;
       float c, ch;

       // left side
       int iq = IV(iz, ix);
       int ip = IU(iz, ix);
       int ir = IV(iz, ix);
       int ir1 = IV(iz, ix + 1);
       int ir2 = IV(iz + 1, ix);
       int ir3 = IV(iz + 1, ix + 1);
       
       C11 = lam[ir] + 2 * mu[ir];
       C13 = lam[ir];

       c  = 1.0 - pml_cx [ix];
       ch = 1.0 - pml_cxh[ix];


       Txx[ip] -= C11 * dt * (c * pml_Pxx[iq] + pml_cx[ix] * DVxx[iq]);
       Tzz[ip] -= C13 * dt * (c * pml_Pxx[iq] + pml_cx[ix] * DVxx[iq]);

       if (iz < nz - 1) {
           C55 = 0.25 * (mu[ir] + mu[ir1] + mu[ir2] + mu[ir3]);
           Txz[ip] -= C55 * dt * (ch * pml_Pzx[iq] + pml_cxh[ix] * DVzx[iq]);
       }
   
       pml_Pxx[iq] = pml_ax [ix] * pml_Pxx[iq] + pml_bx [ix] * DVxx[iq];
       pml_Pzx[iq] = pml_axh[ix] * pml_Pzx[iq] + pml_bxh[ix] * DVzx[iq];


       // right side
       int jx = nxb - 1 - ix;
       int jq = IV(iz, nxb + jx);

       int jx1 = nx - 1 - ix;
       int jx2 = nx - 2 - ix;
       int jp1 = IU(iz, jx1);
       int jp2 = IU(iz, jx2);

       int jr1 = IV(iz, jx1);
       int jr2 = IV(iz, jx2);
       int jr3 = IV(iz + 1, jx2);
       int jr4 = IV(iz + 1, jx1);

       C11 = lam[jr1] + 2 * mu[jr1];
       C13 = lam[jr1];

       Txx[jp1] -= C11 * dt * (c * pml_Pxx[jq] + pml_cx[ix] * DVxx[jq]);
       Tzz[jp1] -= C13 * dt * (c * pml_Pxx[jq] + pml_cx[ix] * DVxx[jq]);

       if (iz < nz - 1) {
           C55 = 0.25 * (mu[jr1] + mu[jr2] + mu[jr3] + mu[jr4]);
           Txz[jp2] -= C55 * dt * (ch * pml_Pzx[jq] + pml_cxh[ix] * DVzx[jq]);
       }

       pml_Pxx[jq] = pml_ax [ix] * pml_Pxx[jq] + pml_bx [ix] * DVxx[jq];
       pml_Pzx[jq] = pml_axh[ix] * pml_Pzx[jq] + pml_bxh[ix] * DVzx[jq];
  }

  
  // z - direction
  #pragma omp parallel for
  for (int ix = 0; ix < nx; ix ++)
  for (int iz = 0; iz < nzb; iz ++) {

       float C11, C13, C55;
       float t, a, b;
       float th, ah, bh;

       float c  = 1.0 - pml_cz [iz];
       float ch = 1.0 - pml_czh[iz];

       // top side
       int iq = IX(iz, ix);
       int ip = IU(iz, ix);
       int ir = IV(iz, ix);
       int ir1 = IV(iz, ix + 1);
       int ir2 = IV(iz + 1, ix);
       int ir3 = IV(iz + 1, ix + 1);

       C11 = lam[ir] + 2 * mu[ir];
       C13 = lam[ir];

       
       if (ifs == 0) {

       Txx[ip] -= C13 * dt * (c * pml_Pzz[iq] + pml_cz[iz] * DVzz[iq]);
       Tzz[ip] -= C11 * dt * (c * pml_Pzz[iq] + pml_cz[iz] * DVzz[iq]);

       if (ix < nx - 1) {
           C55 = 0.25 * (mu[ir] + mu[ir1] + mu[ir2] + mu[ir3]);
           Txz[ip] -= C55 * dt * (ch * pml_Pxz[iq] + pml_czh[iz] * DVxz[iq]);
       }
   
       pml_Pzz[iq] = pml_az [iz] * pml_Pzz[iq] + pml_bz [iz] * DVzz[iq];
       pml_Pxz[iq] = pml_azh[iz] * pml_Pxz[iq] + pml_bzh[iz] * DVxz[iq];
       }
       //

       // bottom side
       int jz = nzb - 1 - iz;
       int jq = IX(nzb + jz, ix);

       int jz1 = nz - 1 - iz;
       int jz2 = nz - 2 - iz;
       int jp1 = IU(jz1, ix);
       int jp2 = IU(jz2, ix);

       int jr1 = IV(jz1, ix);
       int jr2 = IV(jz2, ix);
       int jr3 = IV(jz2, ix + 1);
       int jr4 = IV(jz1, ix + 1);

       C11 = lam[jr1] + 2 * mu[jr1];
       C13 = lam[jr1];

       Txx[jp1] -= C13 * dt * (c * pml_Pzz[jq] + pml_cz[iz] * DVzz[jq]);
       Tzz[jp1] -= C11 * dt * (c * pml_Pzz[jq] + pml_cz[iz] * DVzz[jq]);

       if (ix < nx - 1) {
           C55 = 0.25 * (mu[jr1] + mu[jr2] + mu[jr3] + mu[jr4]);
           Txz[jp2] -= C55 * dt * (ch * pml_Pxz[jq] + pml_czh[iz] * DVxz[jq]);
       }

       pml_Pzz[jq] = pml_az [iz] * pml_Pzz[jq] +  pml_bz [iz] * DVzz[jq];
       pml_Pxz[jq] = pml_azh[iz] * pml_Pxz[jq] +  pml_bzh[iz] * DVxz[jq];
  }
  
}

float fd2d::vpmax (float *lam, float *mu, float *rho, int n) {
  float v, vm = sqrt((lam[0] + 2.0 * mu[0]) / rho[0]);
  for (int i = 1; i < n; i ++) {
       v = sqrt((lam[i] + 2.0 * mu[i]) / rho[i]);
       if (vm < v) vm = v;
  }
  return vm;
}

float fd2d::vsmin (float *mu, float *rho, int n) {
  float v, vm = sqrt(mu[0] / rho[0]);
  for (int i = 1; i < n; i ++) {
       v = sqrt(mu[i] / rho[i]);
       if (vm > v) vm = v;
  }
  return vm;
}

void fd2d::getpmlcoef (float vpmax, float vsmin, float fp) {
  int flag = 1;
  float pd = 2.0, dfac = 3.0;
  float pa = 2.0, afac = 1.0;
  float pb = 2.0, PPW0 = 6.0;    // the 4-th order staggered FD
  float pmlthickx = nxb * dx;
  float pmlthickz = nzb * dz;
  float lnRx = log(10.0) * (-3.0 - (log10(nxb*1.0) - 1.0)/log10(2.0));
  float lnRz = log(10.0) * (-3.0 - (log10(nzb*1.0) - 1.0)/log10(2.0));
  float d0x = -0.5 * (pd + 1.0) * (vpmax / pmlthickx) * lnRx * dfac;
  float d0z = -0.5 * (pd + 1.0) * (vpmax / pmlthickz) * lnRz * dfac;
  float alpha0 = PI * fp * afac;
  float beta0 = vsmin / (0.5 * PPW0 * fp * dx);
  if (beta0 < 1.0) beta0 = 1.0;
  float p, q, d, alpha, beta;
  // alpha = 0.0, beta = 1.0

  for (int i = 0; i < nxb; i++) {
       p = (nxb * 1.0 - i) / nxb;
       d = d0x * pow(p, pd);
       alpha = alpha0 * (1.0 - pow(p, pa));
       beta = 1.0 + (beta0 - 1.0) * pow(p, pb);
       q = alpha + d / beta;
       pml_ax[i]  = exp(-q * dt);
       pml_bx[i]  = (1.0 - pml_ax[i]) * (d / beta) / q;
       pml_cx[i]  = 1.0 - 1.0 / beta;

       p = (nxb - 0.5 - i) / nxb;
       d = d0x * pow(p, pd);
       alpha = alpha0 * (1.0 - pow(p, pa));
       beta = 1.0 + (beta0 - 1.0) * pow(p, pb);
       q = alpha + d / beta;
       pml_axh[i]  = exp(-q * dt);
       pml_bxh[i]  = (1.0 - pml_axh[i]) * (d / beta) / q;
       pml_cxh[i]  = 1.0 - 1.0 / beta;

       //cout << " i=" << i << " pml_c=" << pml_ax[i] << " " << pml_axh[i] << endl;
  }

  for (int i = 0; i < nzb; i++) {
       p = (nzb * 1.0 - i) / nzb;
       d = d0z * pow(p, pd);
       alpha = alpha0 * (1.0 - pow(p, pa));
       beta = 1.0 + (beta0 - 1.0) * pow(p, pb);
       q = alpha + d / beta;
       if (flag == 1) {
           pml_az[i]  = exp(-q * dt);
           pml_bz[i]  = (1.0 - pml_az[i]) * (d / beta) / q;
           pml_cz[i]  = 1.0 - 1.0 / beta;
       }
       else {
           pml_az[i] = q * dt;
           pml_bz[i] = (d / beta) * dt / dz;
           pml_cz[i] = 1.0 / beta;
       }
       //cout << " pml i=" << i << " alpha=" << alpha << " beta=" << beta << " d=" << d << endl;

       p = (nzb - 0.5 - i) / nzb;
       d = d0z * pow(p, pd);
       alpha = alpha0 * (1.0 - pow(p, pa));
       beta = 1.0 + (beta0 - 1.0) * pow(p, pb);
       q = alpha + d / beta;
       if (flag == 1) {
           pml_azh[i]  = exp(-q * dt);
           pml_bzh[i]  = (1.0 - pml_azh[i]) * (d / beta) / q;
           pml_czh[i]  = 1.0 - 1.0 / beta;
       }
       else {
           pml_azh[i] = q * dt;
           pml_bzh[i] = (d / beta) * dt / dz;
           pml_czh[i] = 1.0 / beta;
       }

       //cout << " pml i=" << i << " alpha=" << alpha << " beta=" << beta << " d=" << d << endl;
  }

  /*
  for (int i = 0; i < nzb; i ++)
       cout << " pml i=" << i << " az=" << pml_az[i] << " " << pml_azh[i] << " bz=" << pml_bz[i] << " " << pml_bzh[i] << " cz=" << pml_cz[i] << " " << pml_czh[i] << endl;
  */
        
}

void fd2d::getfdcoef() {
  switch (nor) {
    case 1:
       coef[1] = 1;
       break;
    case 2:
       coef[1] =  9.0/8.0;
       coef[2] = -1.0/24.0;
       break;
    case 3:
       coef[1] =  75.0/ 64.0;
       coef[2] = -25.0/384.0;
       coef[3] =   3.0/640.0;
       break;
    case 4:
       coef[1] =  1.196289;
       coef[2] = -0.0797526;
       coef[3] =  0.009570313;
       coef[4] = -0.0006975447;
       break;
    case 5:
       coef[1] =  1.211243;
       coef[2] = -0.08972168;
       coef[3] =  0.001384277;
       coef[4] = -0.00176566;
       coef[5] = -0.0001186795;
       break;
    default:
       cout << " Warning: no fd coefficients are given ...\n";
       exit(0);
  }
}

void fd2d::freevar() {
   delete [] Vx;
   delete [] Vz;
   delete [] Txx;
   delete [] Txz;
   delete [] Tzz;

   delete [] DVxx;
   delete [] DVxz;
   delete [] DVzx;
   delete [] DVzz;
   delete [] DTxx;
   delete [] DTxz;
   delete [] DTzx;
   delete [] DTzz;

   delete [] pml_ax;
   delete [] pml_bx;
   delete [] pml_cx;
   delete [] pml_az;
   delete [] pml_bz;
   delete [] pml_cz;
   delete [] pml_axh;
   delete [] pml_bxh;
   delete [] pml_cxh;
   delete [] pml_azh;
   delete [] pml_bzh;
   delete [] pml_czh;

   delete [] pml_Qxx;
   delete [] pml_Qxz;
   delete [] pml_Qzx;
   delete [] pml_Qzz;
   delete [] pml_Pxx;
   delete [] pml_Pxz;
   delete [] pml_Pzx;
   delete [] pml_Pzz;
}
