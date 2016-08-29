//#define __MG_DEBUG__
//# define __MG_DEBUG_VMD_REGION__

#include "do_msd0.h"

int MSDComputing(MMatrix_t Wt, int nWt, MMatrix_t bin1, MMatrix2_t msdWt, MMatrix_t msdWtTemp,
		 MMatrix_t posOS, int nOS,
                 int origins,int nmax,
		 MMatrix_t rbin,int nbin,
		 bool bRDF, real L[3], int interval)
{
  int i,j,k,me,me2,l,m,n,indicator;
  int jstart,kend;
  real dis1, dis2;

  for(j=0;j<origins;j=j+interval)
  {
    for(i=0;i<nWt;i++)
    {
      jstart = j*nWt+i;

      //Position at t=0
      if(nbin>1)
      {
	dis1 = GetDisOSOW(Wt, jstart, posOS, nOS, j*nOS, L);
	me = WhichBin(rbin, nbin, dis1);
      }else
      {
	me=0;
      }

      indicator = 1;
  
      ResetMatrix(msdWtTemp);
  
      for(k=0;k<nmax;k++)
      {
	kend = jstart+k*nWt;
	for(l=0; l<DIM; l++)
	{
	  msdWtTemp->m[l][k] += pow(periodicity(Wt->m[kend][l]-Wt->m[jstart][l],L[l]),2);
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

      if(indicator==1)
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
//      printf("MSDComputing for frame %d, water %d\n", j, i);
    }
    printf("MSDComputing for frame %d\n", j);
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
	  msdWt->m[i][j][k] /= bin1->m[i][0];
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
