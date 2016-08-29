#ifndef __MAIN_H__
#define __MAIN_H__


#include "Mgromacs.h"
#include "Mmatrix.h"

void read_index(char *fnNDX, int *ng, int *isize, atom_id **index, char **grpname);
void do_rdf(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    real bin, real Lx, real Ly, real Lz, int natom, int snp);

#endif
