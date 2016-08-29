#include "do_zwtrhocom.h"

//#define __MG_DEBUG__

void do_zwtrhocom(char *fnTPS, char *fnTRX, char *fnRDF, 
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
  MMatrix_t posOS, posOW, posH1, posH2;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);
  posH1 = CreateMatrix(effSectSize*nWt,3);
  posH2 = CreateMatrix(effSectSize*nWt,3);

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
  MMatrix_t zwtrhocom = CreateMatrix(nrbin, 1);

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
	StorePosition(fr, top, posH1, index[emoHW1], isize[emoHW1], nWt, i*nWt, FALSE);
	StorePosition(fr, top, posH2, index[emoHW2], isize[emoHW2], nWt, i*nWt, FALSE);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_ZWTRHOCOM(posOS, posOW, posH1, posH2, zwtrhocom, nOS, nWt, gridsize, L);
  }
  Averaging(zwtrhocom, gridsize, L, totalStep);
            
  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(zwtrhocom, fnRDF, gridsize);

  FreeMatrix(zwtrhocom);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posH1);
  FreeMatrix(posH2);
}

int analyze_ZWTRHOCOM(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posH1, MMatrix_t posH2,
		      MMatrix_t zwtrhocom,
		      int nOS, int nWt, real gridsize, float L[3])
{
  int i,j,k;
  float zSheet, zdis, zpos;
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
    zpos = (posOW->m[i][2]*MASSO + posH1->m[i][2]*MASSH + posH2->m[i][2]*MASSH)/(MASSO+2*MASSH);
    zdis = zpos-zSheet;
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

    zwtrhocom->m[irbin][0] += 1.0;
  }

  return 0;
}

int Averaging(MMatrix_t zwtrhocom, real gridsize, float L[3], int totalStep)
{
  float dV = L[0]*L[1]*gridsize;
  int i, j;

  for(i=0; i<zwtrhocom->rows; i++)
  {
    zwtrhocom->m[i][0] /= (dV*totalStep);
  }

  return 0;
}

int Output(MMatrix_t zwtrhocom, char *fn, real gridsize)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  for(i=0;i<zwtrhocom->rows;i++)
  {
    fprintf(fp,"%15.6f",(i+0.5)*gridsize);
    for(j=0;j<zwtrhocom->cols;j++)
    {
      fprintf(fp,"%15.6f",zwtrhocom->m[i][j]);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  return 0;
}
