//#define __MG_DEBUG__
//# define __MG_DEBUG_VMD_REGION__

#include "do_msd.h"

void do_msd(char *fnTPS, char *fnTRX, char *fnRDF, 
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
  printf("nWt=%d, nOS=%d\n", nWt, nOS);
  
/**
 ** Matrix
 **/
  MMatrix_t posOS;
  MMatrix_t posOW;

  MMatrix_t rbin;
  int nrbin = ReadRBin(&rbin, "MSDBIN.ppa");
  MMatrix_t bin1 = CreateMatrix(nrbin, 1);

  int totalStep = snp;
  real timeInterval = dT;
  int sectionSize  = (int)(chunk/timeInterval);
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int origins = effSectSize/2;        
  int nmax = origins;              
  int section;

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
  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize)
      {
	StorePosition(fr, top, posOW, index[emoSOL], isize[emoSOL], nWt, i*nWt, bMW);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    MSDComputing(posOW, nWt, bin1, msdWt, msdWtTemp,
		 posOS, nOS,
                 origins, nmax, rbin, nrbin, FALSE, L, interval);
#ifdef __MG_DEBUG_VMD_REGION__
    if(section==0 || section==3)
    {
      MSDComputing2(posOW, nWt, bin1, msdWt, msdWtTemp,
		   posOS, nOS,
		   origins, nmax, rSpacing, nrbin, R, FALSE, fr.box[0][0], fr.box[1][1], fr.box[2][2], index);
    }
#endif
//    printf("MSD Computing, Section %d, done\n", section);
  }
  Averaging(bin1, nmax, msdWt, nrbin);
//  printf("Averaging done\n");

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

int MSDComputing2(MMatrix_t Wt, int nWt, MMatrix_t bin1, MMatrix2_t msdWt, MMatrix_t msdWtTemp,
		 MMatrix_t posOS, int nOS,
                 int origins,int nmax,
		 float rSpacing,int nbin,float R, 
		 bool bRDF, real Lx, real Ly, real Lz,
		 atom_id **index)
{
  int i,j,k,me,me2,l,m,n,indicator;
  int jstart,kend;

  real L[3] = {Lx, Ly, Lz};

  real dis1, dis2;

  for(i=0;i<nWt;i++)
  {
    //Position at t=0
    dis1 = GetDisOSOW(Wt, i, posOS, nOS, 0, L);

    me = (int)(dis1/rSpacing);
    if(me >= nbin) me = nbin-1;

    printf("\nme=%10d, index=%10d, dis1=%10.3f\n", me, index[emoSOL][i]+1, dis1);
  }

  return 0;
}
