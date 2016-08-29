#ifndef __MG_MSD_H__
#define __MG_MSD_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoOW, emoHW1, emoHW2, emoNR
};

void do_hbtau2(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real chunk, int interval, real rcut, real acut, bool bMW);

#endif
