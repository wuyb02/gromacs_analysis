#ifndef __MG_OHOH_H__
#define __MG_OHOH_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoOW, emoHW1, emoHW2, emoNR
};

void do_ohoh(char *fnTPS, char *fnTRX, char *fnRDF, 
             int ng, int *isize, char **grpname, atom_id **index,
	     int snp, real dT, real gridsize, real dtheta, real zdis);

#endif
