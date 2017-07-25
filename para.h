#ifndef PARA_H
#define PARA_H
#define PI 3.1415926

class Geometry 
{
  public:
    int nx,ny,nz,nt;
    int nb,noff;           // padding length
    float bx,by,bz;
    float dx,dy,dz,dt;
    int kzmin,kzmax;
    long int ns,nr;        // total number of shots and traces
};

class Shot
{
  public:
    float xo,p;           // plane gather parameter
    float xs,ys,zs;       // source location
    float xg,yg,zg;       // first receiver location
    float dx,dy,dz,dt;
    int ishot;            // shot number
    int nx,ny,nz,nt;      // number of receivers
    int nxb,nyb,nzb;
    int ixs,iys,izs;
    int ixg,iyg,izg;
    int topo;
    int *ixr, *iyr, *izr;
    float *xr, *yr, *zr;
};

class Fwdpar
{
  public:
    int iwvl;                     // type of wavelet
    int isrc;                     // type of source
    float fp;                     // data dominant frequency
    float td;                     // data delay
    float f1,f2,f3,f4;
    float vmin,vmax;

    int iplw;                     // =0, shot; =1, plane-wave
    int ishot1,ishot2;            // shot/pw index
    int dim;                      // =1, 1d; =2, 2d
    int iway;                     // =0; =1, subtract; =2;


    int nproc,nthread;       // parallel computing
    int nxpml,nzpml;         // FD parameters
    int nxor, ntor;          // fd order
    int ifs;                 // =1, free surface; =0, else
};

#endif
