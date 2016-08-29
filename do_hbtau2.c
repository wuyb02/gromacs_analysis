//#define __MG_DEBUG__

#include "do_hbtau2.h"

void do_hbtau2(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, real dT, real chunk, int interval, real rcut, real acut, bool bMW)
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
  bool bTop;

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
  int nrbin = ReadRBin(&rbin, "HBTAURBIN2.ppa");

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
  posH1 = CreateMatrix(effSectSize*nWt,3);
  posH2 = CreateMatrix(effSectSize*nWt,3);
  MMatrix_t donorLife_bin = CreateMatrix(nrbin, nmax);
  MMatrix_t acceptorLife_bin = CreateMatrix(nrbin, nmax);

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
	StorePosition(fr, top, posH1, index[emoHW1], isize[emoHW1], nWt, i*nWt, bMW);
	StorePosition(fr, top, posH2, index[emoHW2], isize[emoHW2], nWt, i*nWt, bMW);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    accumulate_Hbonding(posOS, posOW, posH1, posH2, donorLife_bin, acceptorLife_bin,
			rcut, acut, nrbin, rbin, nWt, nOS, nmax, origins, interval, L);
//    printf("HBTAU Computing, Section %d, done\n", section);
  }
  Averaging(acceptorLife_bin, donorLife_bin, nWt, nrbin, nmax);
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(acceptorLife_bin, donorLife_bin, rbin, dT, nrbin, nmax, nWt, fnRDF);

  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posH1);
  FreeMatrix(posH2);
  FreeMatrix(donorLife_bin);
  FreeMatrix(acceptorLife_bin);
  FreeMatrix(rbin);
}

int Averaging(MMatrix_t acptLife_bin, MMatrix_t donorLife_bin, int nWt, int nbin, int nmax)
{
  int i,j;

  for(i=0;i<nbin;i++)
  {
    if(donorLife_bin->m[i][0]>1e-4)
    {
      for(j=nmax-1;j>=0;j--)
      {
	donorLife_bin->m[i][j] /= donorLife_bin->m[i][0];
      }
    }
    if(acptLife_bin->m[i][0]>1e-4)
    {
      for(j=nmax-1;j>=0;j--)
      {
	acptLife_bin->m[i][j] /= acptLife_bin->m[i][0];
      }
    }
  }

  return 0;
}

int Output(MMatrix_t acptrLife, MMatrix_t donorLife, MMatrix_t rbin, real dt, int nbin, int nmax, int nWt, char *fn)
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
	fprintf(fp,"%15.6e %15.6e %15.6e %15.6e\n",
	        rbin->m[i][0],
		j*dt,
	        donorLife->m[i][j], 
		acptrLife->m[i][j]);
      }
    }
    fclose(fp);
  }

  return 0;
}

int accumulate_Hbonding(MMatrix_t posOS, MMatrix_t posO, MMatrix_t posH1, MMatrix_t posH2,
                        MMatrix_t donorLife_bin, MMatrix_t acceptorLife_bin,
			float rcut, float acut,
			int nbin,MMatrix_t rbin, int nWt, int nOS,
			int nmax, int origins, int interval, float L[3])
{
  int i,j,water,me,k,start;
  float donor[3],hydrogen1[3],hydrogen2[3],acceptor[3];
  int donorList[12],acceptorList[12];    
  float r_dist_O;
  me = 0;

  for(start=0; start<nmax; start=start+interval)
  {
    for(water=0; water<nOS; water++)
    {
      acceptorList[0] = 0;

      // this water molecule as acceptor
      fillAcceptorInformation(acceptor,posOS,start*nOS+water);
      for(j=0;j<nWt;j++)
      {
	fillDonorInformation(donor,hydrogen1,hydrogen2,posO,posH1,posH2,start*nWt+j);
	// check first O-H1 donor pair
	if(isHbond(donor,hydrogen1,acceptor,rcut,acut,L)!=2003
	   || isHbond(donor,hydrogen2,acceptor,rcut,acut,L)!=2003)
	{
	  acceptorLife_bin->m[me][0] += 1.0;
	  acceptorList[++acceptorList[0]] = j;
	}
      }
//      printf("Accumulate List for frame %d, water %d\n",start,water);

      // 2. check for the life time of each HB of the water molecule for each frame
      for(i=1;i<nmax;i++)
      {
	// second check when present water as acceptor;
	fillAcceptorInformation(acceptor,posOS,(i+start)*nOS+water);
	for(k=1;k<acceptorList[0]+1;k++)
	{
	  fillDonorInformation(donor,hydrogen1,hydrogen2,posO,posH1,posH2,(i+start)*nWt+acceptorList[k]);
	  if(isHbond(donor,hydrogen1,acceptor,rcut,acut,L)!=2003 
	     || isHbond(donor,hydrogen2,acceptor,rcut,acut,L)!=2003)
	  {
	    acceptorLife_bin->m[me][i]+=1;
	  }
	}
      }
//      printf("Accumulate HBond for frame %d, water %d\n",start,water);
    }
      
    printf("Accumulate HBond for frame  %d \n",start);
  }

  return 0;
}
