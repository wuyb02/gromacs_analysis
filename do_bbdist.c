#include "do_bbdist.h"

//#define __MG_DEBUG__

void do_bbdist(char *fnTPS, char *fnTRX, char *fnRDF, 
               int ng, int *isize, char **grpname, atom_id **index,
	       int snp)
{
  int i, j, k;
  int ll, mm, nn;

  int g,natoms,j0,j1,n,nframes;
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
  bool bMW = FALSE;

  bTop=read_tps_conf(fnTPS,title,&top,&xdum,NULL,box,bMW);

  int nWt, nOS;
  nOS = isize[emoPEG];

/**
 ** Matrix
 **/
  MMatrix_t posOS;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  MMatrix_t densOS = CreateMatrix(nOS-1,1);
  MMatrix_t mr = CreateMatrix(nOS-1,1);
  MMatrix_t mrt = CreateMatrix(totalSection, nOS-1);

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
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
//    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_BBDIS(posOS, densOS, mr, nOS, L);
    analyze_BBDIST(densOS, mr, mrt->m[section], nOS);
//    printf("MSD Computing, Section %d, done\n", section);
  }
//  Averaging(densOS, mr, nOS);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(densOS, mrt, nOS, totalSection, fnRDF);

  FreeMatrix(posOS);
  FreeMatrix(densOS);
  FreeMatrix(mr);
  FreeMatrix(mrt);
}

int analyze_BBDIS(MMatrix_t posOS, MMatrix_t densOS, MMatrix_t mr,
                  int nOS, float L[3])
{
  int i,j,k;
  real pos1[3], pos2[3], dis;
                                                                                                             
  for(i=0; i<nOS-1; i++)
  {
    for(k=0; k<3; k++)
    {
      pos1[k]=posOS->m[i][k];
      pos2[k]=posOS->m[i+1][k];
    }

    dis = dist1(pos1, pos2, L);
//    printf("i=%d, %10.5f\n", i, dis);

    densOS->m[i][0] += 1.0;
    mr->m[i][0] += dis;
  }

  return 0;
}

int Averaging(MMatrix_t densOS, MMatrix_t mr, int nOS)
{
  int i,j,k;

  for(i=0;i<nOS-1;i++)
  {
    if(densOS->m[i][0]>1e-4)
    {
      mr->m[i][0] /= densOS->m[i][0];
    }
  }

  return 0;
}

int Output(MMatrix_t densOS, MMatrix_t mrt, int nOS, int totalSection, char *fn)
{
  FILE *fp;
  int i,j,k;
  float bbdis;
  
  if(nOS-1>0)
  {
    fp = fopen(fn, "w");
    for(i=0;i<totalSection;i++)
    {
      for(j=0;j<nOS-1;j++)
      {
	fprintf(fp,"%10.6f", mrt->m[i][j]);
      }
      fprintf(fp,"\n");
    }

    fclose(fp);
  }

  return 0;
}

int analyze_BBDIST(MMatrix_t densOS, MMatrix_t mr, double *mrt, int nOS)
{
  int i, j, k;
  real fmr;

  if(nOS>0)
  {
    for(i=0;i<nOS-1;i++)
    {
      fmr=0.0;
      if(densOS->m[i][0]>1e-5)
      {
	fmr=mr->m[i][0]/densOS->m[i][0];
      }
      mrt[i] = fmr;
    }
  }

  return 0;
}
