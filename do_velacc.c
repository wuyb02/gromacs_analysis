#include "do_velacc.h"

void do_velacc(char *fnTPS, char *fnTRX, char *fnRDF, 
            int ng, int *isize, char **grpname, atom_id **index,
	    int snp, real dT, real chunk, bool bMW, float z1, float z2, int nskipt0)
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

  int nWt;
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoSOL], index[emoSOL], grpname[emoSOL]);
  }else
  {
    nWt = isize[emoSOL];
  }
  
/**
 ** Matrix
 **/
  MMatrix_t posOW, velOW;

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

  posOW = CreateMatrix(effSectSize*nWt,3);
  velOW = CreateMatrix(effSectSize*nWt,3);
  MMatrix_t msdWtTemp = CreateMatrix(3,origins);
  MMatrix_t msdWt = CreateMatrix(3,origins);
  MMatrix_t msdWt0 = CreateMatrix(3,1);
  float nhits = 0.0;

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
  int flags = TRX_READ_X | TRX_READ_V;
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
	StoreVelocity(fr, top, velOW, index[emoSOL], isize[emoSOL], nWt, (i+sectionSize/2)*nWt, bMW);
      }
    }
    printf("Section %3.1f of total section %d loaded, nhits=%10.1f\n", 0.5*hsection, totalSection, nhits);
    if(hsection)
    {
      VELACCComputing(posOW, velOW, nWt, &nhits, msdWt, msdWt0, msdWtTemp,
		   origins, nmax, nskipt0, z1, z2, FALSE, L);
    }
    ShiftPosition(posOW, sectionSize/2*nWt);
    ShiftPosition(velOW, sectionSize/2*nWt);
  }
  Averaging(nmax, msdWt, msdWt0, nhits);

  printf("Analysis done\n");
  Output(msdWt, msdWt0, dT, origins, fnRDF, nhits);

  FreeMatrix(posOW);
  FreeMatrix(velOW);
  FreeMatrix(msdWtTemp);
  FreeMatrix(msdWt);
  FreeMatrix(msdWt0);
}


int VELACCComputing(MMatrix_t Wt, MMatrix_t Wtv, int nWt, float *nhits, MMatrix_t msdWt, MMatrix_t msdWt0, MMatrix_t msdWtTemp,
                 int origins,int nmax, int nskipt0,
                 float z1, float z2,
		 bool bRDF, real L[3])
{
  int i,j,k,l,m,n,indicator;
  int jstart,kend;

  for(j=0;j<origins;j=j+nskipt0)
  {
    for(i=0;i<nWt;i++)
    {
      jstart = j*nWt+i;

      ResetMatrix(msdWtTemp);
      indicator = 1;
  
      for(k=0;k<nmax;k++)
      {
	kend = jstart+k*nWt;

	if(Wt->m[kend][2]<z1 || Wt->m[kend][2]>=z2)
	{
	  //printf("Wt->m[%d][2]=%10.3f, z1=%10.3f, z2=%10.3f\n", kend, Wt->m[kend][2], z1, z2);
	  indicator = 0;
	  break;
	}

	for(l=0; l<DIM; l++)
	{
	  msdWtTemp->m[l][k] = Wtv->m[kend][l]*Wtv->m[jstart][l];
	}
      }

      if(indicator==1)
      {
	*nhits += 1.0;
	for(k=0;k<nmax;k++)
	{
	  for(l=0; l<DIM; l++)
	  {
	    msdWt->m[l][k] += msdWtTemp->m[l][k];
	  }
	}
	for(l=0; l<DIM; l++)
	{
	  msdWt0->m[l][0] += msdWtTemp->m[l][0];
	}
      }
    }
    if(j%100==0) printf("j=%d out of %d, nhits=%10.1f\n", j, origins, *nhits);
  }

  return 0;
}

int Averaging(int nmax, MMatrix_t msdWt, MMatrix_t msdWt0, float nhits)
{
  int i,j,k;

  for(i=0;i<DIM;i++)
  {
    for(j=0;j<nmax;j++)
    {
      if((int)nhits!=0)
      {
	msdWt->m[i][j] /= nhits;
      }
    }
    msdWt0->m[i][0] /= nhits;
  }
  
  return 0;
}

int Output(MMatrix_t msdWt, MMatrix_t msdWt0, float dt, int origins, char *fn, float nhits)
{
  FILE *fp;
  int i,j,k;

  fp = fopen(fn,"w");
  for(j=0;j<origins;j++)
  {
    if((int)nhits>0)
    {
      fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",
	      j*dt,
	      (msdWt->m[0][j]+msdWt->m[1][j]+msdWt->m[2][j])/(msdWt->m[0][0]+msdWt->m[1][0]+msdWt->m[2][0]),
	      msdWt->m[0][j]/msdWt->m[0][0],
	      msdWt->m[1][j]/msdWt->m[1][0],
	      msdWt->m[2][j]/msdWt->m[2][0]);
    }else
    {
      fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e\n", j*dt, 0.0, 0.0, 0.0, 0.0);
    }
  }
  fclose(fp);
  
  return 0;
}
