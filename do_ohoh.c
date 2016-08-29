#include "do_ohoh.h"

//#define __MG_DEBUG__

void do_ohoh(char *fnTPS, char *fnTRX, char *fnRDF, 
             int ng, int *isize, char **grpname, atom_id **index,
	     int snp, real dT, real gridsize, real dtheta, real zdis)
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

  int nabin=(int)(180.0/dtheta);
  MMatrix_t ohoh = CreateMatrix(nabin, nabin);

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
    analyze_OHOH(posOS, posOW, posH1, posH2, ohoh, nOS, nWt, gridsize, dtheta, L, zdis);
  }
//  Averaging(ohoh);
            
  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(ohoh, fnRDF);

  FreeMatrix(ohoh);
  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posH1);
  FreeMatrix(posH2);
}

int analyze_OHOH(MMatrix_t posOS, MMatrix_t posOW, MMatrix_t posH1, MMatrix_t posH2,
		 MMatrix_t ohoh,
		 int nOS, int nWt, real gridsize, real dtheta, float L[3], real zdis)
{
  int i,j,k;
  float zSheet, dohc[3], zdir[3]={0, 0, 1.0};
  int iabin1, iabin2;
  double CosAngle1, CosAngle2, theta1, theta2;

  //printf("nOS=%d, nWt=%d, gridsize=%f, zdis=%f, dtheta=%f\n", nOS, nWt, gridsize, zdis, dtheta);

  if(posOS->m[0][2]-posOS->m[1][2]>1e-5)
  {
    printf("Plate z coordinate not the same\n");
    exit(0);
  }
  zSheet = posOS->m[0][2];

  for(i=0; i<nWt; i++)
  {
    if(posOW->m[i][2]-zSheet<zdis-gridsize || posOW->m[i][2]-zSheet>zdis+gridsize)
    {
      continue;
    }

    for(j=0; j<3; j++)
    {
      dohc[j]=periodicity(posH1->m[i][j]-posOW->m[i][j], L[j]);
    }
    CosAngle1 = cos_angle_new(dohc, zdir);
    theta1 = acos(CosAngle1)/PI*180;

    for(j=0; j<3; j++)
    {
      dohc[j]=periodicity(posH2->m[i][j]-posOW->m[i][j], L[j]);
    }
    CosAngle2 = cos_angle_new(dohc, zdir);
    theta2 = acos(CosAngle2)/PI*180;

    if(theta1>185 || theta1<-5 || theta2>185 || theta2<-5)
    {
      printf("CosAngle1=%f, theta1=%f, CosAngle2=%f, theta2=%f\n", CosAngle1, theta1, CosAngle2, theta2);
      exit(0);
    }
    if(theta1>180) theta1=179.9999;
    if(theta2>180) theta2=179.9999;
    if(theta1<0) theta1=0.0;
    if(theta2<0) theta2=0.0;
      
    iabin1 = (int)(theta1/dtheta);
    iabin2 = (int)(theta2/dtheta);

//  printf("CosAngle1=%f, theta1=%f, iabin1=%d, CosAngle2=%f, theta2=%f, iabin2=%d\n", CosAngle1, theta1, iabin1, CosAngle2, theta2, iabin2);

    ohoh->m[iabin1][iabin2] += 1.0;
    ohoh->m[iabin2][iabin1] += 1.0;
  }

  return 0;
}

//int Averaging(MMatrix_t ohoh)
//{
//  float sum;
//  int i, j;
//
//  for(i=0; i<ohoh->rows; i++)
//  {
//    sum = 0.0;
//    for(j=0; j<ohoh->cols; j++)
//    {
//      sum = sum + ohoh->m[i][j];
//    }
//    for(j=0; j<ohoh->cols; j++)
//    {
//      if(sum)
//      {
//	ohoh->m[i][j] /= sum;
//      }
//    }
//
//    printf("i=%d, sum=%f\n", i, sum);
//  }
//
//  return 0;
//}

int Output(MMatrix_t ohoh, char *fn)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  for(i=0;i<ohoh->rows;i++)
  {
    for(j=0;j<ohoh->cols;j++)
    {
      fprintf(fp,"%10.1f",ohoh->m[i][j]);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  return 0;
}
