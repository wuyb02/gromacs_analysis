#include "do_rmsf.h"

//#define __MG_RMSF_DEBUG__

void do_rmsf(char *fnTPS, char *fnTRX, char *fnRDF, 
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
  MMatrix_t posOS, posOS_old;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOS_old = CreateMatrix(effSectSize*nOS,3);
  MMatrix_t densOS = CreateMatrix(nOS,1);
  MMatrix_t mr = CreateMatrix(nOS,3);
  MMatrix_t mr2 = CreateMatrix(nOS,3);
  MMatrix2_t rmsf = CreateMatrix2(totalSection,nOS,7);

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
	if(section*sectionSize+i>0)
	{
	  MatrixCopy(posOS_old, posOS);
	}
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
//    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_RMSF(posOS, posOS_old, densOS, mr, mr2, nOS, L);
    CalcRMSF(densOS, mr, mr2, rmsf->m[section], nOS);
//    printf("MSD Computing, Section %d, done\n", section);
  }
  Averaging(densOS, mr, mr2, nOS);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(densOS, mr, mr2, rmsf, nOS, totalSection, fnRDF);
//  printf("Output done\n");

  FreeMatrix(posOS);
  FreeMatrix(posOS_old);
  FreeMatrix(densOS);
  FreeMatrix(mr);
  FreeMatrix(mr2);
  FreeMatrix2(rmsf);
}

int analyze_RMSF(MMatrix_t posOS, MMatrix_t posOS_old, MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, 
                 int nOS, float L[3])
{
  int i,j,k,me,meAng,meAccept,m;
  real pos[3], ppos[3], diff[3];

  for(k=0; k<3; k++)
    ppos[k] = 0.0;
                                                                                                             
  for(i=0; i<nOS; i++)
  {
    for(k=0; k<3; k++)
    {
      pos[k]=posOS->m[i][k];
    }

    if(densOS->m[i][0]>1e-4)
    {
      for(k=0; k<3; k++)
      {
	ppos[k]=posOS_old->m[i][k];
	diff[k]=periodicity(pos[k]-ppos[k], L[k]);
	pos[k]=ppos[k]+diff[k];
      }
    }

//    printf("%10d %15.4f %15.4f %15.4f, %15.4f %15.4f %15.4f\n", 
//	   i,
//	   pos[0], pos[1], pos[2],
//	   mpos[0], mpos[1], mpos[2]);
    densOS->m[i][0] += 1.0;
    for(k=0; k<3; k++)
    {
      mr->m[i][k] += pos[k];
      mr2->m[i][k] += pos[k]*pos[k];
      posOS->m[i][k]=pos[k];
    }
  }

  return 0;
}

int CalcRMSF(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, double **rmsf, int nOS)
{
  int i, j, k;
  real fmr, fmr2;

  if(nOS>0)
  {
    for(i=0;i<nOS;i++)
    {
      for(k=0; k<3; k++)
      {
	fmr=0.0; fmr2=0.0;
	if(densOS->m[i][0]>1e-5)
	{
	  fmr=mr->m[i][k]/densOS->m[i][0];
	  fmr2=mr2->m[i][k]/densOS->m[i][0];
	}
	rmsf[i][6] += (fmr2-fmr*fmr);
	rmsf[i][k] = fmr;
	rmsf[i][k+3] = fmr2;
      }
//      rmsf[i] = sqrt(rmsf[i]);
    }
  }

  return 0;
}

int Averaging(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, int nOS)
{
  int i,j,k;

  for(i=0;i<nOS;i++)
  {
    if(densOS->m[i][0]>1e-4)
    {
      for(k=0;k<3;k++)
      {
	mr->m[i][k] /= densOS->m[i][0];
	mr2->m[i][k] /= densOS->m[i][0];
      }
    }
  }

  return 0;
}

int Output(MMatrix_t densOS, MMatrix_t mr, MMatrix_t mr2, MMatrix2_t rmsf, int nOS, int totalSection, char *fn)
{
  FILE *fp;
  int i,j,k;
  float frmsf;
  
  if(nOS>0)
  {
    fp = fopen(fn, "w");
#ifndef __MG_RMSF_DEBUG__
    for(i=0;i<nOS;i++)
    {
      for(k=0, frmsf=0.0; k<3; k++)
      {
	frmsf += (mr2->m[i][k]-mr->m[i][k]*mr->m[i][k]);
      }
      frmsf = sqrt(frmsf);
      fprintf(fp,"%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
	      densOS->m[i][0],
	      frmsf,
	      pow(mr->m[i][0],2), 
	      pow(mr->m[i][1],2), 
	      pow(mr->m[i][2],2), 
	      mr2->m[i][0], 
	      mr2->m[i][1], 
	      mr2->m[i][2]);
    }
#else
    for(i=0; i<totalSection; i++)
    {
      for(j=0; j<nOS; j++)
      {
	for(k=0; k<7; k++)
	{
	  fprintf(fp, "%10.6f", rmsf->m[i][j][k]);
	}
      }
      fprintf(fp, "\n");
    }
#endif

    fclose(fp);
  }

  return 0;
}
