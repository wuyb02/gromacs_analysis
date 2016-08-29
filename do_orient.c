#include "do_orient.h"

void do_orient(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, int interval)
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
  bool bMW=FALSE, bTop;

  bTop=read_tps_conf(fnTPS,title,&top,&xdum,NULL,box,bMW);

  int nWt, nOS, nAT, nMO;
  nWt = isize[emoSOL];
  nMO = MolInGRP(top, isize[emoSOL], index[emoSOL], grpname[emoSOL]);
  nAT = nWt/nMO;
  nOS = isize[emoPEG];
  printf("\nPEG, %d, nWt, %d, nMO, %d, nAT, %d\n", nOS, nWt, nMO, nAT);
  
/**
 ** Matrix
 **/
  MMatrix_t posOS;
  MMatrix_t posOW;

  MMatrix_t rbin;
  int nrbin = ReadRBin(&rbin, "ORIENTBIN.ppa");
  MMatrix_t bin1 = CreateMatrix(nrbin, 1);
  MMatrix_t binRC = CreateMatrix(nrbin, 1);

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

  posOS = CreateMatrix(nOS,3);
  posOW = CreateMatrix(nWt,3);
  MMatrix_t vecRC = CreateMatrix(effSectSize*nWt, 3); //Rotational Correction
  MMatrix_t corRCTemp = CreateMatrix(3,origins);
  MMatrix2_t corRC = CreateMatrix2(nrbin, 3, nmax);

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

  for(section=0;section<totalSection;section++) 
  {
    for(i=0;i<sectionSize;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize)
      { 
//	int num = RemovePBCSep2(fr, top, isize[emoSOL], index[emoSOL], grpname[emoSOL]);
//	printf("step=%d, PBC Shift, num=%d\n", section*totalSection+i, num);
	StorePosition(fr, top, posOW, index[emoSOL], isize[emoSOL], nWt, 0, FALSE);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, 0, FALSE);
	StoreVec(vecRC, posOW, nMO, nAT, i);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    autoCorrelationRC(corRC,corRCTemp,vecRC,binRC,nMO,origins,nmax,interval);
  }

//  printf("binRC, %5.3f\n", binRC->m[0][0]);
  Averaging(binRC, nmax, corRC, nrbin);

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(corRC, rbin, dT, origins, nrbin, nWt, fnRDF);

  FreeMatrix(bin1);
  FreeMatrix(binRC);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(vecRC);
  FreeMatrix(corRCTemp);
  FreeMatrix2(corRC);
}

/*
int StorePosition(t_trxframe fr, MMatrix_t posOW, int *OWIndex, int nWt)
{
  int i, j, k;

  if(nWt==0)
  {
    return 0;
  }
  if(posOW->rows<nWt)
  {
    printf("do_orient.c, posOW->rows<nWt\n");
    exit(0);
  }
  if(posOW->cols<3)
  {
    printf("do_orient.c, posOW->cols<3\n");
    exit(0);
  }

  for(i=0;i<nWt;i++)
  {
    posOW->m[i][0] = fr.x[OWIndex[i]][0];
    posOW->m[i][1] = fr.x[OWIndex[i]][1];
    posOW->m[i][2] = fr.x[OWIndex[i]][2];
  }

  return 0;
}
*/

void StoreVec(MMatrix_t vecRC, MMatrix_t posOW, int nMO, int nAT, int step)
{
  int i, j, k, ll, mm, nn;

  MMatrix_t aa=CreateMatrix(3,1); //Rotational Vector, Characteristic direction
  MMatrix_t A=CreateMatrix(nAT,3);

  if(vecRC->cols!=aa->rows)
  {
    printf("do_orient.c, StoreVec(), vecRC->cols!=aa->rows, %d, %d\n", vecRC->cols, aa->rows);
    exit(0);
  }

  for(i=0; i<nMO; i++)
  {
    for(j=0; j<A->rows; j++)
    {
      for(k=0; k<A->cols; k++)
      {
	A->m[j][k] = posOW->m[i*nAT+j][k];
      }
    }
    GetOrientation(A, aa);
//    GetEEOrientation(A, aa);

    for(j=0; j<aa->rows; j++)
    {
      ll = step*nMO+i;
      vecRC->m[ll][j] = aa->m[j][0];
    }
  }

  FreeMatrix(aa);
  FreeMatrix(A);
}

void autoCorrelationRC(MMatrix2_t corRC, MMatrix_t corRCTemp, MMatrix_t vecRC, MMatrix_t binRC, 
                       int nMO, int origins, int nmax, int interval)
{
  int i,j,k,l,m,mz;
  int jstart,kend;
  int indicator;

  if(nMO == 0)
  {
    return;
  }

  for(j=0;j<origins;j=j+interval)
  {
    for(i=0; i<nMO; i++)
    {
	jstart = j*nMO+i;

	indicator = 1;
	ResetMatrix(corRCTemp);
  //      if(DA[i][2]>z1 && DA[i][2]<z2) 
  //      {
  //        mz = (DA[i][2]-z1)/zsize;
	mz=0;
      
	for(k=0; k<nmax; k++)
	{
	  kend = jstart+k*nMO;
	  for(l=0; l<3; l++)
	  {
	    corRCTemp->m[l][k] += vecRC->m[jstart][l]*vecRC->m[kend][l];
	  }
	}

	if(indicator==1)
	{
	  binRC->m[mz][0]+=1.0;
	  for(k=0; k<nmax; k++)
	  {
	    for(l=0; l<3; l++)
	    {
	      corRC->m[mz][l][k] += corRCTemp->m[l][k];
	    }
	  }
	}

  //      }
    }
    printf("Orientation Autocorrelation for frame %d\n", j);
  }
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
	  msdWt->m[i][j][k] /= (bin1->m[i][0]);
	}
      }
    }
  }
  
  return 0;
}

int Output(MMatrix2_t msdWt, MMatrix_t rbin, float dt, int origins, int nbin, int nWt, char *fn)
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
	fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
		rbin->m[i][0],
		j*dt,
		msdWt->m[i][0][j]+msdWt->m[i][1][j]+msdWt->m[i][2][j],
		msdWt->m[i][0][j], 
		msdWt->m[i][1][j],
		msdWt->m[i][2][j]);
      }
    }
    fclose(fp);
  }
  
  return 0;
}
