#include "do_cntnion.h"

//#define __MG_DEBUG__

void do_cntnion(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real zlcnt, real zucnt)
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

  int nWt;
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoOW], index[emoOW], grpname[emoOW]);
  }else
  {
    nWt = isize[emoOW];
  }

/**
 ** Matrix
 **/
  MMatrix_t posOW;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

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

  real ntoccupied, ntroccupied, nion;
  MMatrix_t cntnion = CreateMatrix(nWt, snp);
  MMatrixInt_t zion_start = CreateMatrixInt(nWt, 1);
  MMatrixInt_t zion_now = CreateMatrixInt(nWt, 1);
  MMatrixInt_t zion_prev = CreateMatrixInt(nWt, 1);
  MMatrixInt_t tion_instart = CreateMatrixInt(nWt, 1);

/**
 ** ANALYSIS
 **/
  ntoccupied = 0.0;
  nion = 0.0;
  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize)
      {
	StorePosition(fr, top, posOW, index[emoOW], isize[emoOW], nWt, i*nWt, FALSE);
      }
    }
    //printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_CNTNION(posOW, cntnion, &ntoccupied, &ntroccupied, &nion, zion_start, zion_now, zion_prev, tion_instart, nWt, section, L, zlcnt, zucnt, dT);
  }
  //Averaging(cntnion, gridsize, L, totalStep);
            
  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(cntnion, fnRDF);
  printf("Ion occupancy %10.3f, Translocated ion occupancy %10.7f, # of translocation events %10.0f\n", 
      ntoccupied/totalStep, ntroccupied/totalStep, nion);

  FreeMatrix(cntnion);
  FreeMatrix(posOW);
}

int analyze_CNTNION(MMatrix_t posOW,
		   MMatrix_t cntnion, real *ntoccupied, real *ntroccupied, real *nion,
		   MMatrixInt_t zion_start, MMatrixInt_t zion_now, MMatrixInt_t zion_prev, MMatrixInt_t tion_instart,
		   int nWt, int section, float L[3], real zlcnt, real zucnt, real dT)
{
  int i,j,k;

  for(i=0; i<nWt; i++)
  {
    cntnion->m[i][section] = posOW->m[i][2];
    while(cntnion->m[i][section]>=L[2])
    {
      cntnion->m[i][section] -= L[2];
    }
    while(cntnion->m[i][section]<0.0)
    {
      cntnion->m[i][section] += L[2];
    }
    if(cntnion->m[i][section]<zlcnt)
    {
      cntnion->m[i][section] = zlcnt-0.00001;
      zion_now->m[i][0] = -1;
    }else if(cntnion->m[i][section]>zucnt)
    {
      cntnion->m[i][section] = zucnt+0.00001;
      zion_now->m[i][0] = 1;
    }else
    {
      zion_now->m[i][0] = 0;
    }
  }

  for(i=0; i<nWt; i++)
  {
    if(cntnion->m[i][section]>zlcnt && cntnion->m[i][section]<zucnt)
    {
      *ntoccupied = *ntoccupied + 1.0;
      break;
    }
  }

  if(!section)
  {
    for(i=0; i<nWt; i++)
    {
      zion_start->m[i][0] = zion_now->m[i][0];
      zion_prev->m[i][0] = zion_now->m[i][0];
    }
  }

  if(section)
  {
    for(i=0; i<nWt; i++)
    {
      if(zion_now->m[i][0]!=0)
      {
	if(zion_prev->m[i][0]==0 && zion_now->m[i][0]!=zion_start->m[i][0])
	{
	  *nion = *nion + 1.0;
	  *ntroccupied = *ntroccupied + (section-tion_instart->m[i][0]);
	  printf("Translocation Event, %10d, T=%10.4f, deltaT_inside_ion=%10.4f, %10s, Total %10.0f\n", 
	      i+1, section*dT, (section-tion_instart->m[i][0])*dT, (zion_start->m[i][0]<0? "UP":"DOWN"), *nion);
	}
	zion_start->m[i][0] = zion_now->m[i][0];
	tion_instart->m[i][0] = section;
      }
      zion_prev->m[i][0] = zion_now->m[i][0];
    }
  }

  return 0;
}

//int Averaging(MMatrix_t cntnion, real gridsize, float L[3], int totalStep)
//{
//  float dV = L[0]*L[1]*gridsize;
//  int i, j;
//
//  for(i=0; i<cntnion->rows; i++)
//  {
//    cntnion->m[i][0] /= (dV*totalStep);
//  }
//
//  return 0;
//}

int Output(MMatrix_t cntnion, char *fn)
{
  FILE *fp;
  int i,j,k;
  
  fp = fopen(fn, "w");

  for(i=0;i<cntnion->rows;i++)
  {
    for(j=0;j<cntnion->cols;j++)
    {
//      if(j%10)
//      {
//	continue;
//      }
      fprintf(fp,"%15.6f",cntnion->m[i][j]);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  return 0;
}
