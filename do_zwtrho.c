#include "do_zwtrho.h"

//#define __MG_DEBUG__

void do_zwtrho(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real gridsize)
{
  int i, j, k;
  int ll, mm, nn;

  int g,natoms,j0,j1,n,nframes;
  int **count;
  int isize_cm=0,nrdf=0,max_i;
  atom_id *index_cm=NULL;

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
//  LJSigmaByAtomName(char *atomname);

  int nOS, nWt;
  nOS = isize[emoPEG];
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoOW], index[emoOW], grpname[emoOW]);
  }else
  {
    nWt = isize[emoOW];
  }
//  PrintAtomType(top, isize[emoPEG], index[emoPEG], grpname[emoPEG]);

/**
 ** Matrix
 **/
  MMatrix_t posOS, posOW;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);

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

  int nrbin=(int)(L[2]/gridsize)+1;
  MMatrix_t zwtrho = CreateMatrix(nrbin, 1);

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
	StorePosition(fr, top, posOW, index[emoOW], isize[emoOW], nWt, i*nWt, FALSE);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_ZWTRHO(posOS, posOW, zwtrho, nOS, nWt, gridsize, L);
  }
  Averaging(zwtrho, gridsize, L, totalStep);
            
  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(zwtrho, fnRDF, gridsize);

  FreeMatrix(zwtrho);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
}

int analyze_ZWTRHO(MMatrix_t posOS, MMatrix_t posOW,
		   MMatrix_t zwtrho,
		   int nOS, int nWt, real gridsize, float L[3])
{
  int i,j,k;
  float zSheet, zdis;
  int irbin;

//  printf("nOS=%d, nWt=%d, gridsize=%f\n", nOS, nWt, gridsize);

  if(posOS->m[0][2]-posOS->m[1][2]>1e-5)
  {
    printf("Plate z coordinate not the same\n");
    exit(0);
  }
  zSheet = posOS->m[0][2];

  for(i=0; i<nWt; i++)
  {
    zdis = posOW->m[i][2]-zSheet;
    if(zdis<0)
    {
      zdis = zdis + L[2];
    }
    irbin = (int)(zdis/gridsize);

    if(irbin>(int)(L[2]/gridsize))
    {
      printf("irbin=%d, %d\n", irbin, (int)(L[2]/gridsize));
      exit(0);
    }

//    printf("i=%d, zSheet=%f, posOW=%f, irbin=%d\n", i, zSheet, posOW->m[i][2], irbin);

    zwtrho->m[irbin][0] += 1.0;
  }

  return 0;
}

int Averaging(MMatrix_t zwtrho, real gridsize, float L[3], int totalStep)
{
  float dV = L[0]*L[1]*gridsize;
  int i, j;

  for(i=0; i<zwtrho->rows; i++)
  {
    zwtrho->m[i][0] /= (dV*totalStep);
  }

  return 0;
}

int Output(MMatrix_t zwtrho, char *fn, real gridsize)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  for(i=0;i<zwtrho->rows;i++)
  {
    fprintf(fp,"%15.6f",(i+0.5)*gridsize);
    for(j=0;j<zwtrho->cols;j++)
    {
      fprintf(fp,"%15.6f",zwtrho->m[i][j]);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  return 0;
}
