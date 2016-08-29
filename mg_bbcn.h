#ifndef __MG_HBOND_H__
#define __MG_HBOND_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoNR
};

void do_bbcn(char *fnTPS, char *fnTRX, char *fnRDF, 
             int ng, int *isize, char **grpname, atom_id **index,
	     int snp);

#endif
