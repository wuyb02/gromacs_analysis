/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_CA_H__
#define __DO_CA_H__

#include "do_xxx.h"
#include "Mgromacs.h"
#include "Mmatrix.h"

enum {
  emoPEG, emoOW, emoHW1, emoHW2, emoNR
};

void do_ca(char *fnTPS, char *fnTRX, char *fnRDF, char *fnGridR, char *fnGridZ,
           int ng, int *isize, char **grpname, atom_id **index,
	   int snp, real dT, real zg, real Lz, real R, real zbin, real dA);

int analyze_CA(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posHW1, MMatrix_t posHW2,
	       MMatrix_t wtrho,
	       int nWt, int section, float L[3], real zg, real Lz, real R, real zbin, real dA);
int Averaging(MMatrix_t wtrho, float dV, int totalStep);
int Output(MMatrix_t wtrho, float dA, float zbin, float zg, char *fn, char *fnGridR, char *fnGridZ);

#endif
