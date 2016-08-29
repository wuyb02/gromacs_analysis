#ifndef __MG_MSD_H__
#define __MG_MSD_H__


#include "Mgromacs.h"
#include "Mmatrix.h"
#include "read_index.h"

enum {
  emoPEG, emoSOL, emoNR
};
//enum {
//  emoPEG, emoSOL, emoNR
//};

void do_msd2(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, bool bMW, int interval);

#endif
