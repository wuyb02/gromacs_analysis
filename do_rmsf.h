/**
 ** do_rmsf.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of rmsfation
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
 ** rmsfWtTemp	: matrix, Temp MSD, DIM*origins
 ** rmsfWt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_rmsf.h"

int analyze_RMSF(MMatrix_t posOS, MMatrix_t posOS_old, MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, int nOS, float L[3]);
int Averaging(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, int nOS);
int Output(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, MMatrix2_t rmsf, int nOS, int totalSection, char *fnRDF);
int CalcRMSF(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, double **rmsf, int nOS);
