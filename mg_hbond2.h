#ifndef __MG_HBOND_H__
#define __MG_HBOND_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoOW, emoH1, emoH2, emoNR
};

void do_hbond2(char *fnTPS, char *fnTRX, char *fnRDF, 
               int ng, int *isize, char **grpname, atom_id **index,
	       int snp, bool bMW, int nAng, real dT, real tUpdate, real rcut, real acut, real rlist);

#endif
