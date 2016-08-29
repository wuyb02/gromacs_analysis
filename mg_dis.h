#ifndef __MG_DIS_H__
#define __MG_DIS_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoSOL, emoNR
};

void do_dis(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real gridsize);

#endif
