//#define __MG_DEBUG__
//# define __MG_DEBUG_VMD_REGION__

#include "do_msd.h"

void do_msd(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, bool bMW)
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
                 origins, nmax, rbin, nrbin, FALSE, fr.box[0][0], fr.box[1][1], fr.box[2][2]);
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
}

int MSDComputing(MMatrix_t Wt, int nWt, MMatrix_t bin1, MMatrix2_t msdWt, MMatrix_t msdWtTemp,
		 MMatrix_t posOS, int nOS,
                 int origins,int nmax,
		 MMatrix_t rbin,int nbin,
		 bool bRDF, real Lx, real Ly, real Lz)
{
  int i,j,k,me,me2,l,m,n,indicator;
  int jstart,kend;

  real L[3] = {Lx, Ly, Lz};

  real dis1, dis2;

  for(i=0;i<nWt;i++)
  {
    //Position at t=0
    dis1 = GetDisOSOW(Wt, i, posOS, nOS, 0, L);
    me = WhichBin(rbin, nbin, dis1);

    indicator = 0;
//    if(dis1<R)
//    {
      indicator = 1;

      ResetMatrix(msdWtTemp);

      for(j=0;j<origins;j++)
      {
	jstart = j*nWt+i;
	for(k=0;k<nmax;k++)
	{
	  kend = jstart+k*nWt;
	  for(l=0; l<DIM; l++)
	  {
	    msdWtTemp->m[l][k] += pow(periodicity(Wt->m[kend][l]-Wt->m[jstart][l],L[l]),2);
	  }
	}
      }

//      dis2 = GetDisOSOW(Wt, (origins+nmax-1)*nWt+i, posOS, nOS, (origins+nmax-1)*nOS, L);
//      me2 = WhichBin(rbin, nbin, dis2);
//      if(me2>=nbin) me2 = nbin-1;
//      if(me2!=me)
//      {
//	printf("dis2/rSpacing!=me, %d, %d, %10.3f, %10.3f\n", (int)(dis2/rSpacing), me, dis1, dis2);
//	indicator=0;
//      }

//      printf("dis1=%15.3f, me=%d, dis2=%15.3f, me2=%d\n", dis1, me, dis2, me2);

      if(indicator == 1)
      {
	bin1->m[me][0] += 1.0;
	for(k=0;k<nmax;k++)
	{
	  for(l=0; l<DIM; l++)
	  {
	    msdWt->m[me][l][k] += msdWtTemp->m[l][k];
	  }
	}
      }

//      printf("indicator=%d\n", indicator);
//    }
  }

  return 0;
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

int Averaging(MMatrix_t bin1, int nmax, MMatrix2_t msdWt, int nbin)
{
  int i,j,k;

  for(i=0;i<nbin;i++)
  {
    for(j=0;j<3;j++)
    {
      for(k=0;k<nmax;k++)
      {
	if(bin1->m[i][0]!=0)
	{
	  msdWt->m[i][j][k] /= (bin1->m[i][0]*nmax);
	}
      }
    }
  }
  
  return 0;
}

int Output(MMatrix2_t msdWt, MMatrix_t rbin, MMatrix_t bin1, float dt, int origins, int nbin, int nWt, char *fn)
{
  FILE *fp;
  int i,j,k;

  if(nWt>0)
  {
    fp = fopen(fn,"w");
    for(i=0;i<nbin;i++)
    {
      for(j=0;j<origins;j++)
      {
	fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
		rbin->m[i][0],
		j*dt,
		msdWt->m[i][0][j]+msdWt->m[i][1][j]+msdWt->m[i][2][j],
		msdWt->m[i][0][j], 
		msdWt->m[i][1][j],
		msdWt->m[i][2][j],
		bin1->m[i][0]);
      }
    }
    fclose(fp);
  }
  
  return 0;
}
