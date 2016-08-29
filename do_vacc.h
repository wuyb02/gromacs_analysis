/**
 ** do_hbond.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of hbondation
 **   StorePosition
 **   RemovePBCSep, pbc=full, Rhodamine is seperated by PBC
 **
 ** To Speed up
 **  Use bin for RDF NS
 **
 ** totalStep = snp;
 ** timeInterval = dT;
 ** sectionSize  = (int)(chunk/timeInterval);
 ** effSectSize  = sectionSize;       
 ** totalSection = totalStep/sectionSize;
 ** origins = effSectSize/2;        
 ** nmax = origins;              
 **
 ** nrbin	: # of bins
 ** rbin	: matrix, bin radius, nrbin*1
 ** bin1	: matrix, hits in bin, nrbin*1
 **
 ** nOS		: # of PEG atoms
 ** posOS	: matrix, PEG atoms coordinates, (nOS*effSectSize)*DIM
 ** nWt		: bMW=true, # of water molecules
 **		: bMW=false, # of water atoms
 ** posOW	: matrix, bMW=true, center of mass of water, (nWt*effSectSize)*DIM
 **		: matrix, bMW=false, water atoms coordinates, (nWt*effSectSize)*DIM
 **
 ** hbondWtTemp	: matrix, Temp MSD, DIM*origins
 ** hbondWt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_vacc.h"

int analyze_VACC(MMatrix_t posOS, MMatrix_t posOW, 
                 MMatrix_t grids, MMatrixInt_t vacc, MMatrix_t fvacc,
		 t_topology top, int size, atom_id *index, char *grpname,
                 int nOS, real solutesize, int step, float L[3], int Ln[3]);
//int Averaging(MMatrix_t densWt, MMatrix_t hBond, MMatrix_t h2Bond, MMatrix_t OH_Oangle, int nbin, int nAng);
int Output(MMatrix_t fvacc, int nOS, char *fn, int totalStep);