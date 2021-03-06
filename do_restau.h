/**
 ** do_hbtau.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of hbtauation
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
 ** hbtauWtTemp	: matrix, Temp MSD, DIM*origins
 ** hbtauWt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_restau.h"

int accumulate_RESTAU(MMatrix_t posOS, MMatrix_t posO, 
		      MMatrix_t resLife_bin,
		      int nbin, MMatrix_t rbin, int nWt, int nOS,
		      int nmax, int origins, int interval, float L[3]);
int Averaging(MMatrix_t resLife_bin, int nWt, int nbin, int nmax);
int Output(MMatrix_t resLife, MMatrix_t rbin, real dt, int nbin, int nmax, int nWt, char *fn);
