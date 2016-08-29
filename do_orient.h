/**
 ** do_orient.c
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007
 **
 ** Calculate Autocorrelation of orientation
 **   StorePosition
 **   RemovePBCSep, pbc=full, Rhodamine is seperated by PBC
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
 ** binRC	: matrix, hits in bin, nrbin*1
 **
 ** nOS		: # of PEG atoms
 ** posOS	: matrix, PEG atoms coordinates, nOS*DIM
 ** nWt		: bMW=true, # of water molecules
 **		: bMW=false, # of water atoms
 ** posOW	: matrix, bMW=true, center of mass of water, nWt*DIM
 **		: matrix, bMW=false, water atoms coordinates, nWt*DIM
 **
 ** msdWtTemp	: matrix, Temp MSD, DIM*origins
 ** msdWt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_orient.h"

#define PBCSEP 0.6

void StoreVec(MMatrix_t vecRC, MMatrix_t posOW, int nMO, int nAT, int i);
void autoCorrelationRC(MMatrix2_t corRC, MMatrix_t corRCTemp, MMatrix_t vecRC, MMatrix_t binRC, 
                       int nMO, int origins, int nmax, int interval);
int Averaging(MMatrix_t bin1, int nmax, MMatrix2_t msdWt, int nbin);
int Output(MMatrix2_t msdWt, MMatrix_t rbin, float dt, int origins, int nbin, int nWt, char *fn);
