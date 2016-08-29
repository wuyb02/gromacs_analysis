/**
 ** do_hbtau2.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/28/2007, increasing statistics by 2
 **
 ** Calculate Autocorrelation of hbtau2ation
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
 ** hbtau2WtTemp	: matrix, Temp MSD, DIM*origins
 ** hbtau2Wt	: matrix2, MSD, nrbin*DIM*origins
 **/

#include "do_xxx.h"
#include "mg_hbtau2.h"

int accumulate_Hbonding(MMatrix_t posOS, MMatrix_t posO, MMatrix_t posH1, MMatrix_t posH2,
                        MMatrix_t donorLife_bin, MMatrix_t acceptorLife_bin,
			float rcut, float acut,
			int nbin,MMatrix_t rbin, int nWt, int nOS, 
			int nmax, int origin, int interval, float L[3]);
int Averaging(MMatrix_t acptLife_bin, MMatrix_t donorLife_bin, int nWt, int nbin, int nmax);
int Output(MMatrix_t acptrLife, MMatrix_t donorLife, MMatrix_t rbin, real dt, int nbin, int nmax, int nWt, char *fn);
