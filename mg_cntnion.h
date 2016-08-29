#ifndef __MG_CNTNION_H__
#define __MG_CNTNION_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoOW, emoNR
};

void do_cntnion(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real zlcnt, real zucnt);

#endif
