#include "do_dis.h"

//#define __MG_DEBUG__

void do_dis(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real gridsize)
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
//  LJSigmaByAtomName(char *atomname);

  int nOS, nWt;
  nOS = isize[emoPEG];
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoSOL], index[emoSOL], grpname[emoSOL]);
  }else
  {
    nWt = isize[emoSOL];
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

  MMatrixInt_t dis = CreateMatrixInt(totalStep,nWt);

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
	StorePosition(fr, top, posOW, index[emoSOL], isize[emoSOL], nWt, i*nWt, FALSE);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_DIS(posOS, posOW, dis, nOS, nWt, gridsize, section, L);
  }
//  Averaging(densWt, hBond, hBond2, OH_Oangle, nrbin, nAng);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(dis, nOS, fnRDF, totalStep, nWt, gridsize);

  FreeMatrixInt(dis);
}

int analyze_DIS(MMatrix_t posOS, MMatrix_t posOW, 
                MMatrixInt_t dis,
                int nOS, int nWt, real gridsize, int step, float L[3])
{
  int i,j,k;
  float r_dist_O;
  int iOS;

  printf("nOS=%d, nWt=%d, gridsize=%f\n", nOS, nWt, gridsize);

  for(i=0; i<nWt; i++)
  {
    GetDisOSOW2(posOW, i, posOS, nOS, &r_dist_O, &iOS, 0, L);
    dis->m[step][i]=(int)(r_dist_O/gridsize);
  }

//  for(i=0; i<nWt; i++)
//  {
//    printf("step=%10d, dis=%10d\n", step, dis->m[step][i]);
//  }

  return 0;
}

int Averaging(MMatrixInt_t dis, int totalStep, int nWt)
{
  return 0;
}

int Output(MMatrixInt_t dis, int nOS, char *fn, int totalStep, int nWt, real gridsize)
{
  FILE *fp;
  int i,j,k;
  
  if(nOS>0)
  {
    fp = fopen(fn, "w");
    for(i=0;i<totalStep;i++)
    {
      for(j=0;j<nWt;j++)
      {
	fprintf(fp,"%15.6f",dis->m[i][j]*gridsize);
      }
      fprintf(fp,"\n");
    }

    fclose(fp);
  }

  return 0;
}
