#include "do_wtcomxy.h"

//#define __MG_DEBUG__

void do_wtcomxy(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real z1, real z2)
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
  MMatrix_t posOS, posOW, posHW1, posHW2;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);
  posHW1 = CreateMatrix(effSectSize*nWt,3);
  posHW2 = CreateMatrix(effSectSize*nWt,3);

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

  MMatrix_t Xcom = CreateMatrix(1, totalSection);
  MMatrix_t Ycom = CreateMatrix(1, totalSection);

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
	StorePosition(fr, top, posHW1, index[emoHW1], isize[emoHW1], nWt, i*nWt, FALSE);
	StorePosition(fr, top, posHW2, index[emoHW2], isize[emoOW], nWt, i*nWt, FALSE);
      }
    }
//    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_WtCOMXY(posOS, posOW, posHW1, posHW2, Xcom, Ycom, nOS, nWt, section, L, z1, z2);
  }
            
  printf("Analysis done\n");
  Output(Xcom, Ycom, fnRDF, dT);

  FreeMatrix(Xcom);
  FreeMatrix(Ycom);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posHW1);
  FreeMatrix(posHW2);
}

int analyze_WtCOMXY(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posHW1, MMatrix_t posHW2,
		    MMatrix_t Xcom, MMatrix_t Ycom,
		    int nOS, int nWt, int section, float L[3], real z1, real z2)
{
  int i,j,k;
  real zcom[3], totalmass;

  totalmass=0.0;
  for(i=0; i<nWt; i++)
  {
    for(k=0; k<3; k++)
    {
      zcom[k] = (posOW->m[i][k]*MASSO + posHW1->m[i][k]*MASSH + posHW2->m[i][k]*MASSH) / (MASSO+2*MASSH);
      zcom[k] = ZRemovePBC(zcom[k], L[k]);
    }

    if(zcom[2]>z1 && zcom[2]<z2)
    {
      totalmass += 1.0;
      Xcom->m[0][section] += zcom[0];
      Ycom->m[0][section] += zcom[1];
    }
//    printf("step %d\n", section);
//    printf("posOW , %10.3f%10.3f%10.3f\n", posOW->m[i][0], posOW->m[i][1], posOW->m[i][2]);
//    printf("posHW1, %10.3f%10.3f%10.3f\n", posHW1->m[i][0], posHW1->m[i][1], posHW1->m[i][2]);
//    printf("posHW2, %10.3f%10.3f%10.3f\n", posHW2->m[i][0], posHW2->m[i][1], posHW2->m[i][2]);
//    printf("com   , %10.3f%10.3f%10.3f\n", zcom[0], zcom[1], zcom[2]);
//    printf("Xcom  , %10.3f%10.3f\n", Xcom->m[0][section], Ycom->m[0][section]);
//    printf("\n");
  }
  
  if(totalmass<1e-5)
  {
    printf("section=%d, totalmass=%10.5f\n", section, totalmass);
    exit(0);
  }
  Xcom->m[0][section] /= totalmass;
  Ycom->m[0][section] /= totalmass;

  return 0;
}

int Output(MMatrix_t Xcom, MMatrix_t Ycom, char *fn, real dT)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  for(i=0;i<Xcom->cols;i++)
  {
    fprintf(fp,"%15.6f%15.6f%15.6f\n",i*dT, Xcom->m[0][i], Ycom->m[0][i]);
  }

  fclose(fp);

  return 0;
}