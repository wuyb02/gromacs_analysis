/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_WTCOMXY_H__
#define __DO_WTCOMXY_H__

#include "do_xxx.h"
#include "Mgromacs.h"
#include "Mmatrix.h"

enum {
  emoPEG, emoOW, emoHW1, emoHW2, emoNR
};

void do_wtcomxy(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real z1, real z2);

int analyze_WtCOMXY(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posHW1, MMatrix_t posHW2,
		    MMatrix_t Xcom, MMatrix_t Ycom,
		    int nOS, int nWt, int section, float L[3], real z1, real z2);
int Output(MMatrix_t Xcom, MMatrix_t Ycom, char *fn, real dT);

#endif
