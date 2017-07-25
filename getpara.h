#ifndef GETPARA_H
#define GETPARA_H
#include "para.h"

class GetPara
{
  private:
    int nx, ny, nz, nt, ns;
    float dx, dy, dz, dt;
    float bx, by, bz;

    char vf1[100];
    char vf2[100];
    char vf3[100];
    char headf[100];
    char dataf[100];

    void readpar (char *, char *, char *, char *, char *, char *, Geometry *geomd, Fwdpar *par);    
    int next_atom( char **startp, char* atom);

  public:
    GetPara (char *filename, Geometry *geo,Fwdpar *par);
    void GetHeader (Shot *src, int myid);
    void GetAcouVel (float *v1, float *v2, float *v3);
    void GetElasVel (float *v1, float *v2, float *v3);
    void MapVpToLam (float *, float *, float *, float *, float *);
    void PrintPara (Fwdpar par);
};

#endif
