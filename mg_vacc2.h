#ifndef __MG_VACC_H__
#define __MG_VACC_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoNR
};

void do_vacc2(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real gridsize, real solutesize);

#endif
