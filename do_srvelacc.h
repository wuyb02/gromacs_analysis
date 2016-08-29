#ifndef __DO_SRVELACC_H__
#define __DO_SRVELACC_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "do_xxx.h"

enum {
  emoSOL, emoNR
};

void do_srvelacc(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, bool bMW, float xc, float yc, float zc, float r1, float r2, int nskipt0);

int SRVELACCComputing(MMatrix_t Wt, MMatrix_t Wtv, int nWt, float *nhits, MMatrix_t msdWt, MMatrix_t msdWt0, MMatrix_t msdWtTemp,
                 int origins,int nmax, int nskipt0,
                 float xc, float yc, float zc, float r1, float r2,
		 bool bRDF, real L[3]);
int Averaging(int nmax, MMatrix_t msdWt, MMatrix_t msdWt0, float nhits);
int Output(MMatrix_t msdWt, MMatrix_t msdWt0, float dt, int origins, char *fn, float nhits);


#endif
