#include "do_msd2.h"

void do_msd2(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, bool bMW, int interval)
{
  int i, j, k;

  int g,natoms,nbin,j0,j1,n,nframes;
  int **count;
  int isize_cm=0,nrdf=0,max_i;
  atom_id *index_cm=NULL;
  unsigned long int *sum;

/**
 ** Topology
 **   top.atoms.nr,
 **   top.atoms.atom[i].m,
 **   top.atoms.atom[i].q,
 **   top.atoms.atom[i].type,
 **   top.atoms.atom[i].resnr
 **   top.atoms.nres
 **/
  t_topology top;
  char title[256];
  rvec *xdum;
  matrix box;
  bool bTop;

  bTop=read_tps_conf(fnTPS,title,&top,&xdum,NULL,box,bMW);

  int nWt, nOS;
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoSOL], index[emoSOL], grpname[emoSOL]);
  }else
  {
    nWt = isize[emoSOL];
  }
  nOS = isize[emoPEG];
  
/**
 ** Matrix
 **/
  MMatrix_t posOS;
  MMatrix_t posOW;

  MMatrix_t rbin;
  int nrbin = ReadRBin(&rbin, "MSD2BIN.ppa");
  MMatrix_t bin1 = CreateMatrix(nrbin, 1);

  int totalStep = snp;
  real timeInterval = dT;
  int sectionSize  = (int)(chunk/timeInterval);
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int origins = effSectSize/2;        
  int nmax = origins;              
  int hsection;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);
  printf("# of origins : %d\n",origins);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);
  MMatrix_t msdWtTemp = CreateMatrix(3,origins);
  MMatrix2_t msdWt = CreateMatrix2(nrbin,3,origins);

/**
 ** T_BLOCK
 **/
  t_block    *excl;
  excl=NULL;

/**
 ** Trajectory
 **/
  t_trxframe fr;
  int status;
  int flags = TRX_READ_X;
  read_first_frame(&status, fnTRX, &fr, flags);

  real L[3]={fr.box[0][0], fr.box[1][1], fr.box[2][2]};

/**
 ** ANALYSIS
 **/
  if(sectionSize%2)
  {
    printf("sectionSize%%2!=0, %d\n", sectionSize);
    exit(0);
  }
  for(hsection=0;hsection<totalSection*2;hsection++)
  {
    for(i=0;i<sectionSize/2;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize/2)
      {
	StorePosition(fr, top, posOW, index[emoSOL], isize[emoSOL], nWt, (i+sectionSize/2)*nWt, bMW);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, (i+sectionSize/2)*nOS, FALSE);
      }
    }
    printf("Section %3.1f of total section %d loaded\n", 0.5*hsection, totalSection);
    if(hsection)
    {
      MSDComputing(posOW, nWt, bin1, msdWt, msdWtTemp,
		   posOS, nOS,
		   origins, nmax, rbin, nrbin, FALSE, L, interval);
    }
    ShiftPosition(posOS, sectionSize/2*nOS);
    ShiftPosition(posOW, sectionSize/2*nWt);
  }
  Averaging(bin1, nmax, msdWt, nrbin);

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(msdWt, rbin, bin1, timeInterval, origins, nrbin, nWt, fnRDF);

  FreeMatrix(rbin);
  FreeMatrix(bin1);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(msdWtTemp);
  FreeMatrix2(msdWt);
}
