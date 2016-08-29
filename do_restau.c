#include "do_restau.h"

void do_restau(char *fnTPS, char *fnTRX, char *fnRDF, 
               int ng, int *isize, char **grpname, atom_id **index,
	       int snp, real dT, real chunk, int interval)
{
  int i, j, k;

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

  int nWt, nOS;
  if(bMW==TRUE)
  {
    nWt = MolInGRP(top, isize[emoOW], index[emoOW], grpname[emoOW]);
  }else
  {
    nWt = isize[emoOW];
  }
  nOS = isize[emoPEG];
  
/**
 ** Matrix
 **/
  MMatrix_t posOS, posOW, posH1, posH2;

  MMatrix_t rbin;
  int nrbin = ReadRBin(&rbin, "RESRBIN.ppa");

  int totalStep = snp;
  real timeInterval = dT;
  int sectionSize  = (int)(chunk/timeInterval);
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int origins = effSectSize/2;        
  int nmax = origins;              
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);
  printf("# of origins : %d\n",origins);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);
  MMatrix_t resLife_bin = CreateMatrix(nrbin, nmax);

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
	StorePosition(fr, top, posOW, index[emoOW], isize[emoOW], nWt, i*nWt, bMW);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    accumulate_RESTAU(posOS, posOW, resLife_bin,
		      nrbin, rbin, nWt, nOS, nmax, origins, interval, L);
//    printf("HBTAU Computing, Section %d, done\n", section);
  }
  Averaging(resLife_bin, nWt, nrbin, nmax);
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(resLife_bin, rbin, dT, nrbin, nmax, nWt, fnRDF);
}

int Averaging(MMatrix_t resLife_bin, int nWt, int nbin, int nmax)
{
  int i,j;

  for(i=0;i<nbin;i++)
  {
    if(resLife_bin->m[i][0]>1e-4)
    {
      printf("bin%d, %10.5f\n", i, resLife_bin->m[i][0]);
      for(j=nmax-1;j>=0;j--)
      {
	resLife_bin->m[i][j] /= resLife_bin->m[i][0];
      }
    }
  }

  return 0;
}

int Output(MMatrix_t resLife, MMatrix_t rbin, real dt, int nbin, int nmax, int nWt, char *fn)
{
  FILE *fp;
  int i,j;

  if(nWt>0)
  {
    fp = fopen(fn, "w");
    for(i=0;i<nbin;i++)
    {
      for(j=0;j<nmax;j++)
      {
	fprintf(fp,"%15.6e %15.6e %15.6e\n",
	        rbin->m[i][0],
		j*dt,
	        resLife->m[i][j]);
      }
    }
    fclose(fp);
  }

  return 0;
}

int accumulate_RESTAU(MMatrix_t posOS, MMatrix_t posO, 
		      MMatrix_t resLife_bin,
		      int nbin, MMatrix_t rbin, int nWt, int nOS,
		      int nmax, int origins, int interval, float L[3])
{
  int i,j,water,me,me1,start;
  float r_dist_O, r_dist_O1;

  for(start=0; start<nmax; start=start+interval)
  {
    for(water=0; water<nWt; water++)
    {
      r_dist_O = GetDisOSOW(posO, start*nWt+water, posOS, nOS, 0, L);
      me = WhichBin(rbin, nbin, r_dist_O);
      resLife_bin->m[me][0]+=1.0;

      for(i=1;i<nmax;i++)
      {
	r_dist_O1 = GetDisOSOW(posO, (start+i)*nWt+water, posOS, nOS, (start+i)*nOS, L);
	me1 = WhichBin(rbin, nbin, r_dist_O);

	printf("start=%d, water=%d, me=%d, r=%10.4f, me1=%d, r1=%10.4f",
	       start, water, me, r_dist_O, me1, r_dist_O1);

	if(me1==me)
	{
	  resLife_bin->m[me][i]+=1.0;
	}//else
//	{
//	  break;
//	}
      }
    }

    printf("Accumulate HBond for frame  %d \n",start);
  }

  return 0;
}
