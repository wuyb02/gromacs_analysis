/**
 ** do_msd.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of msdation
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
 ** msdWtTemp	: matrix, Temp MSD, DIM*origins
 ** msdWt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_msd0.h"

#ifdef __MG_DEBUG_VMD_REGION__
int MSDComputing2(MMatrix_t Wt, int nWt, MMatrix_t bin1, MMatrix2_t msdWt, MMatrix_t msdWtTemp,
		 MMatrix_t posOS, int nOS,
                 int origins,int nmax,
		 float rSpacing,int nbin,float R, 
		 bool bRDF, real Lx, real Ly, real Lz,
		 atom_id **index);
#endif
