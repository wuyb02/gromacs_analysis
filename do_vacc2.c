#include "do_vacc2.h"

//#define __MG_DEBUG__

void do_vacc2(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real gridsize, real solutesize)
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

  int nOS;
  nOS = isize[emoPEG];
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
  posOW = CreateMatrix(1,3);

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
  int Ln[3]={(int)(fr.box[0][0]/gridsize), (int)(fr.box[1][1]/gridsize), (int)(fr.box[2][2]/gridsize)};

  MMatrixInt_t vacc = CreateMatrixInt(Ln[0]*Ln[1]*Ln[2],1);
  MMatrix_t fvacc = CreateMatrix(totalStep,1);

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
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_VACC(posOS, posOW, vacc, fvacc, top, isize[emoPEG], index[emoPEG], grpname[emoPEG],
                 nOS, solutesize, gridsize, section, L, Ln);
  }
//  Averaging(densWt, hBond, hBond2, OH_Oangle, nrbin, nAng);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(fvacc, nOS, fnRDF, totalStep);

  FreeMatrix(posOS);
  FreeMatrixInt(vacc);
}

int analyze_VACC(MMatrix_t posOS, MMatrix_t posOW, 
                 MMatrixInt_t vacc, MMatrix_t fvacc,
		 t_topology top, int size, atom_id *index, char *grpname,
                 int nOS, real solutesize, real gridsize, int step, float L[3], int Ln[3])
{
  int i,j,k;
  float r_dist_O, sigma;
  int iOS;
                                                                                                             
  for(i=0;i<Ln[0]*Ln[1]*Ln[2];i++)
  {
    posOW->m[0][0]=(i/(Ln[1]*Ln[2])+0.5)*gridsize;
    posOW->m[0][1]=((i/Ln[2])%Ln[1]+0.5)*gridsize;
    posOW->m[0][2]=(i%Ln[2]+0.5)*gridsize;
 
    GetDisOSOW2(posOW, 0, posOS, nOS, &r_dist_O, &iOS, 0, L);

    if(iOS>=size)
    {
      printf("do_vacc.c, Analyze_VACC, iOS>=size, %d, %d\n", iOS, size);
      exit(0);
    }
    sigma = LJSigmaByAtomName(*(top.atoms.atomtype[index[iOS]]));
//    printf("%10s%10.3f\n", *(top.atoms.atomtype[index[iOS]]), sigma);

    if(r_dist_O>0.5*(sigma+solutesize))
    {
      vacc->m[i][0]=1;
      fvacc->m[step][0]+=1.0;
    }else
    {
      vacc->m[i][0]=0;
    }

    printf("Finish%7dof %7d,%6.3f,%6.3f,%6.3f, From index[%5d]%5d%5s,%6.3f,%6.3f,%6.3f,%6.3f,result is %3d\n", 
	   i, Ln[0]*Ln[1]*Ln[2], 
	   posOW->m[0][0], posOW->m[0][1], posOW->m[0][2],
	   iOS, index[iOS], *(top.atoms.atomtype[index[iOS]]), posOS->m[iOS][0], posOS->m[iOS][1], posOS->m[iOS][2],
	   r_dist_O, vacc->m[i][0]);
  }

  fvacc->m[step][0] /= Ln[0]*Ln[1]*Ln[2];

  return 0;
}

int Averaging(MMatrix_t densWt, MMatrix_t hBond, MMatrix_t h2Bond, MMatrix_t OH_Oangle,
            int nbin, int nAng)
{
  int i,j;

  for(i=0;i<nbin;i++)
  {
    if(hBond->m[i][0]!=0)
    {
      for(j=0;j<nAng;j++)
      {
	OH_Oangle->m[i][j] /= hBond->m[i][0];
      }
    }
    if(densWt->m[i][0]!=0)
    {
      hBond->m[i][0] /= densWt->m[i][0];
      h2Bond->m[i][0] /= densWt->m[i][0];
    }
  }

  return 0;
}

int Output(MMatrix_t fvacc, int nOS, char *fn, int totalStep)
{
  FILE *fp;
  int i,j,k;
  
  if(nOS>0)
  {
    fp = fopen(fn, "w");
    for(i=0;i<totalStep;i++)
    {
      fprintf(fp,"%15.6f\n",
	      fvacc->m[i][0]);
    }

    fclose(fp);
  }

  return 0;
}
