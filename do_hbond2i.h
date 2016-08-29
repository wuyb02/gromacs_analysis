/**
 ** do_hbond2i.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of hbond2ation
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
 ** hbond2WtTemp	: matrix, Temp MSD, DIM*origins
 ** hbond2Wt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_hbond2i.h"

int analyze_Hbonding(MMatrix_t posOS, MMatrix_t posO,MMatrix_t posH1,MMatrix_t posH2,
                     MMatrix_t densWt, MMatrix_t hBond, MMatrix2_t OH_Oangle, MMatrix_t hBond2, MMatrixInt_t list,
                     int nWt, int nOS, int nbin, MMatrix_t rbin, float L[3], float rcut,float acut, int nAng);
int Averaging(MMatrix_t densWt, MMatrix_t hBond, MMatrix_t h2Bond, MMatrix2_t OH_Oangle, int nbin, int nAng, int nOS);
int Output(MMatrix_t rbin, MMatrix_t densW, MMatrix_t hBond, MMatrix_t hBond2, MMatrix2_t OH_Oangle,
           int nbin, int nWt, char *fn, int nAng, int nOS);
