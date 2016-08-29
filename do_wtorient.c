#include "do_wtorient.h"

//#define __MG_DEBUG__

void do_wtorient(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real gridsize, real dtheta)
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

  int nrbin=(int)(L[2]/gridsize)+1, nabin=(int)(180.0/dtheta);
  MMatrix_t wtorient = CreateMatrix(nrbin, nabin);

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
    analyze_WTORIENT(posOS, posOW, posH1, posH2, wtorient, nOS, nWt, gridsize, dtheta, L);
  }
  Averaging(wtorient);
            
  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(wtorient, fnRDF, gridsize, dtheta);

  FreeMatrix(wtorient);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posH1);
  FreeMatrix(posH2);
}

int analyze_WTORIENT(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posH1, MMatrix_t posH2,
		     MMatrix_t wtorient,
		     int nOS, int nWt, real gridsize, real dtheta, float L[3])
{
  int i,j,k;
  float zSheet, dohc[3], zdir[3]={0, 0, 1.0};
  int irbin, iabin;
  double CosAngle, theta;

  printf("nOS=%d, nWt=%d, gridsize=%f, dtheta=%f\n", nOS, nWt, gridsize, dtheta);

  if(posOS->m[0][2]-posOS->m[1][2]>1e-5)
  {
    printf("Plate z coordinate not the same\n");
    exit(0);
  }
  zSheet = posOS->m[0][2];

  for(i=0; i<nWt; i++)
  {
    irbin = (int)((posOW->m[i][2]-zSheet)/gridsize);

    for(j=0; j<3; j++)
    {
      dohc[j]=periodicity((posH1->m[i][j]+posH2->m[i][j])/2-posOW->m[i][j], L[j]);
    }
    CosAngle = cos_angle_new(dohc, zdir);
    theta = acos(CosAngle)/PI*180;
    iabin = (int)(theta/dtheta);

    if(irbin>(int)(L[2]/gridsize) || iabin>(int)(180/dtheta))
    {
      printf("irbin=%d, %d, iabin=%d, %d\n", irbin, (int)(L[2]/gridsize), iabin, (int)(180/dtheta));
      exit(0);
    }

    printf("irbin=%d, CosAngle=%f, theta=%f, iabin=%d\n", irbin, CosAngle, theta, iabin);

    wtorient->m[irbin][iabin] += 1.0;
  }

  return 0;
}

int Averaging(MMatrix_t wtorient)
{
  float sum;
  int i, j;

  for(i=0; i<wtorient->rows; i++)
  {
    sum = 0.0;
    for(j=0; j<wtorient->cols; j++)
    {
      sum = sum + wtorient->m[i][j];
    }
    for(j=0; j<wtorient->cols; j++)
    {
      if(sum)
      {
	wtorient->m[i][j] /= sum;
      }
    }

    printf("i=%d, sum=%f\n", i, sum);
  }

  return 0;
}

int Output(MMatrix_t wtorient, char *fn, real gridsize, real dtheta)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  fprintf(fp,"%15.6f",0.0);
  for(j=0;j<wtorient->cols;j++)
  {
    fprintf(fp,"%15.6f",(j+0.5)*dtheta);
  }
  fprintf(fp,"\n");

  for(i=0;i<wtorient->rows;i++)
  {
    fprintf(fp,"%15.6f",(i+0.5)*gridsize);
    for(j=0;j<wtorient->cols;j++)
    {
      fprintf(fp,"%15.6f",wtorient->m[i][j]);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  return 0;
}
