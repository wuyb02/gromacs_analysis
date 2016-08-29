#include "do_bbcn.h"

//#define __MG_DEBUG__

void do_bbcn(char *fnTPS, char *fnTRX, char *fnRDF, 
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
  bool bTop, bMW=FALSE;

  bTop=read_tps_conf(fnTPS,title,&top,&xdum,NULL,box,bMW);

  int nOS;
  nOS = isize[emoPEG];

/**
 ** Matrix
 **/
  MMatrix_t posOS;

  MMatrix_t rbin;
  int nrbin = 1;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS, 3);
  MMatrix_t dens = CreateMatrix(nOS-1,nOS-1);
  MMatrix_t vecOS = CreateMatrix(nOS-1, 3);
  MMatrix_t bbCn = CreateMatrix(nOS-1, nOS-1);

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
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, 0, FALSE);
	StoreVec(posOS, vecOS, nOS, L);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
//    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_BBCN(vecOS, dens, bbCn, nOS, L);
//    printf("MSD Computing, Section %d, done\n", section);
  }
  Averaging(dens, bbCn, nOS);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(dens, bbCn, nOS, fnRDF);

  FreeMatrix(posOS);
  FreeMatrix(dens);
  FreeMatrix(bbCn);
  FreeMatrix(vecOS);
}

int analyze_BBCN(MMatrix_t vecOS, MMatrix_t dens, MMatrix_t bbCn, int nOS, float L[3])
{
  int i,j,k;
  float vec1[3], vec2[3], result;
                                                                                                             
  for(i=0; i<nOS-1; i++)
  {
    for(k=0; k<3; k++)
    {
      vec1[k] = vecOS->m[i][k];
    }
    for(j=0; j<nOS-1; j++)
    {
      for(k=0; k<3; k++)
      {
	vec2[k] = vecOS->m[j][k];
      }

      result = i_prod(vec1, vec2);

      bbCn->m[i][j] += result;
      dens->m[i][j] += 1.0;
    }
  }

  return 0;
}

int Averaging(MMatrix_t dens, MMatrix_t bbCn, int nOS)
{
  int i,j;

  for(i=0; i<nOS-1; i++)
  {
    for(j=0; j<nOS-1; j++)
    {
      if(dens->m[i][j]>1e-5)
      {
	bbCn->m[i][j] /= dens->m[i][j];
      }
    }
  }

  return 0;
}

int Output(MMatrix_t dens, MMatrix_t bbCn, int nOS, char *fn)
{
  FILE *fp;
  int i,j,k;
  
  if(nOS>0)
  {
    fp = fopen(fn, "w");

    for(i=0;i<nOS-1;i++)
    {
      for(j=0; j<nOS-1; j++)
      {
	fprintf(fp,"%15.6f", bbCn->m[i][j]);
	printf("%15.6f", dens->m[i][j]);
      }
      fprintf(fp,"\n");
      printf("\n");
    }

    fclose(fp);
  }

  return 0;
}

void StoreVec(MMatrix_t posOS, MMatrix_t vecOS, int nOS, float L[3])
{
  int i, j, k, ll, mm, nn;
  real vec[3];

  for(i=0; i<nOS-1; i++)
  {
    ll = i;
    mm = i+1;

    for(k=0; k<3; k++)
    {
      vec[k]=periodicity(posOS->m[mm][k]-posOS->m[ll][k], L[k]);
    }
    NormVec(vec);

    for(k=0; k<3; k++)
    {
      vecOS->m[i][k] =vec[k];
    }
  }
}
