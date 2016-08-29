#include "do_ca.h"

void do_ca(char *fnTPS, char *fnTRX, char *fnRDF, char *fnGridR, char *fnGridZ,
           int ng, int *isize, char **grpname, atom_id **index,
	   int snp, real dT, real zg, real Lz, real R, real zbin, real dA)
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

  int nzbin, nrbin;
  nzbin = (int)((Lz-zg)/zbin);
  Lz = zg + zbin*nzbin;
  nrbin = (int)(R*R/dA*PI);
  R = sqrt(nrbin*dA/PI);

/**
 ** Trajectory
 **/
  t_trxframe fr;
  int status;
  int flags = TRX_READ_X;
  read_first_frame(&status, fnTRX, &fr, flags);

  real L[3]={fr.box[0][0], fr.box[1][1], fr.box[2][2]};

  printf("%10.3f%10.3f%10.3f\n", L[0], L[1], L[2]);

  MMatrix_t wtrho = CreateMatrix(nzbin, nrbin);

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
    //printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_CA(posOS, posOW, posHW1, posHW2, wtrho, nWt, section, L, zg, Lz, R, zbin, dA);
  }
            
  printf("Analysis done\n");

  Averaging(wtrho, zbin*dA, totalStep);
  printf("Averaging done\n");

  Output(wtrho, dA, zbin, zg, fnRDF, fnGridR, fnGridZ);
  printf("Output done\n");

  FreeMatrix(wtrho);
  printf("Free wtrho\n");
  FreeMatrix(posOS);
  printf("Free posOS\n");
  FreeMatrix(posOW);
  printf("Free posOW\n");
  FreeMatrix(posHW1);
  printf("Free posHW1\n");
  FreeMatrix(posHW2);
  printf("Free posHW2\n");
}

int analyze_CA(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posHW1, MMatrix_t posHW2,
	       MMatrix_t wtrho,
	       int nWt, int section, float L[3], real zg, real Lz, real R, real zbin, real dA)
{
  int i,j,k;
  real zcom[3], r2, totalmass, Xcom, Ycom;
  int irbin, izbin;

  totalmass=0.0;
  Xcom = 0.0;
  Ycom = 0.0;
  for(i=0; i<nWt; i++)
  {
    for(k=0; k<3; k++)
    {
      zcom[k] = (posOW->m[i][k]*MASSO + posHW1->m[i][k]*MASSH + posHW2->m[i][k]*MASSH) / (MASSO+2*MASSH);
      zcom[k] = ZRemovePBC(zcom[k], L[k]);
    }

    r2 = (zcom[0]-L[0]/2)*(zcom[0]-L[0]/2) + (zcom[1]-L[1]/2)*(zcom[1]-L[1]/2);
    if(r2<R*R && zcom[2]<Lz && zcom[2]>zg)
    {
      Xcom += zcom[0];
      Ycom += zcom[1];
      totalmass += 1.0;
    }
  }

  if(totalmass<1e-5)
  {
    printf("section=%d, totalmass=%10.5f\n", section, totalmass);
    exit(0);
  }
  Xcom /= totalmass;
  Ycom /= totalmass;

//  printf("(Xcom, Ycom), (%10.3f, %10.3f)\n", Xcom, Ycom);

  for(i=0; i<nWt; i++)
  {
    for(k=0; k<3; k++)
    {
      zcom[k] = (posOW->m[i][k]*MASSO + posHW1->m[i][k]*MASSH + posHW2->m[i][k]*MASSH) / (MASSO+2*MASSH);
      zcom[k] = ZRemovePBC(zcom[k], L[k]);
    }

    r2 = (zcom[0]-Xcom)*(zcom[0]-Xcom) + (zcom[1]-Ycom)*(zcom[1]-Ycom);
    if(r2<R*R && zcom[2]<Lz && zcom[2]>zg)
    {
      izbin = (int)((zcom[2]-zg)/zbin);
      irbin = (int)(r2*PI/dA);

      wtrho->m[izbin][irbin] += 1.0;
    }
  }

  return 0;
}

int Averaging(MMatrix_t wtrho, float dV, int totalStep)
{
  int i, j, k;

  for(i=0; i<wtrho->rows; i++)
  {
    for(j=0; j<wtrho->cols; j++)
    {
      wtrho->m[i][j] = wtrho->m[i][j]*(MASSO+2*MASSH)/totalStep/dV/0.602;
    }
  }

  return 0;
}

int Output(MMatrix_t wtrho, float dA, float zbin, float zg, char *fn, char *fnGridR, char *fnGridZ)
{
  FILE *fp, *fpr, *fpz;
  int i,j,k;
  
  fp = fopen(fn, "w");
  fpr = fopen(fnGridR, "w");
  fpz = fopen(fnGridZ, "w");

  for(i=0; i<wtrho->rows; i++)
  {
    for(j=0; j<wtrho->cols; j++)
    {
      fprintf(fp, "%15.6f", wtrho->m[i][j]);
    }
    fprintf(fp, "\n");
  }

  for(i=0; i<wtrho->cols; i++)
  {
    fprintf(fpr, "%15.6f\n", (sqrt(i*dA/PI)+sqrt((i+1)*dA/PI))/2);
  }

  for(i=0; i<wtrho->rows; i++)
  {
    fprintf(fpz, "%15.6f\n", zbin*(i+i+1)/2+zg);
  }

  fclose(fp);
  fclose(fpr);
  fclose(fpz);

  return 0;
}
