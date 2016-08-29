#include "do_hbond.h"

//#define __MG_DEBUG__

void do_hbond(char *fnTPS, char *fnTRX, char *fnRDF, 
              int ng, int *isize, char **grpname, atom_id **index,
	      int snp, bool bMW, int nAng, real dT, real tUpdate, real rcut, real acut, real rlist)
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
  int nrbin = ReadRBin(&rbin, "HBRBIN.ppa");

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;
  int freqUpdate = (int)(tUpdate/dT);
  if(freqUpdate==0)
  {
    freqUpdate=1;
  }

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

  posOS = CreateMatrix(effSectSize*nOS,3);
  posOW = CreateMatrix(effSectSize*nWt,3);
  posH1 = CreateMatrix(effSectSize*nWt,3);
  posH2 = CreateMatrix(effSectSize*nWt,3);
  MMatrix_t densWt = CreateMatrix(nrbin,1);
  MMatrix_t hBond = CreateMatrix(nrbin,1);
  MMatrix_t hBond2 = CreateMatrix(nrbin,1);
  MMatrix_t OH_Oangle = CreateMatrix(nrbin,nAng);
  MMatrixInt_t list = CreateMatrixInt(nWt, maxList);
  real rlist2 = rlist*rlist;

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
	StorePosition(fr, top, posH1, index[emoH1], isize[emoH1], nWt, i*nWt, bMW);
	StorePosition(fr, top, posH2, index[emoH2], isize[emoH2], nWt, i*nWt, bMW);
	StorePosition(fr, top, posOS, index[emoPEG], isize[emoPEG], nOS, i*nOS, FALSE);
//	printf("StorePosition, Section %d, frame %d, done\n", section, i);
      }
      if((section*sectionSize+i)%freqUpdate==0)
	neighbourList(nWt,list,rlist2,posOW,L);
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    analyze_Hbonding(posOS, posOW, posH1, posH2,
                     densWt, hBond, OH_Oangle, hBond2, list,
                     nWt, nOS, nrbin, rbin, L, rcut, acut, nAng);
//    printf("MSD Computing, Section %d, done\n", section);
  }
  Averaging(densWt, hBond, hBond2, OH_Oangle, nrbin, nAng);
            
//  printf("Averaging done\n");

  printf("Analysis done\n");
//  printf("nrbin = %d\n", nrbin);
  Output(rbin, densWt, hBond, hBond2, OH_Oangle, nrbin, nWt, fnRDF, nAng);

  FreeMatrix(posOS);
  FreeMatrix(posOW);
  FreeMatrix(posH1);
  FreeMatrix(posH2);
  FreeMatrix(rbin);
  FreeMatrix(densWt);
  FreeMatrix(hBond);
  FreeMatrix(hBond2);
  FreeMatrix(OH_Oangle);
  FreeMatrixInt(list);
}

int analyze_Hbonding(MMatrix_t posOS, MMatrix_t posO,MMatrix_t posH1,MMatrix_t posH2,
                     MMatrix_t densWt, MMatrix_t hBond, MMatrix_t OH_Oangle, MMatrix_t hBond2, MMatrixInt_t list,
                     int nWt, int nOS, int nbin, MMatrix_t rbin, float L[3], float rcut,float acut, int nAng)
{
  int i,j,me,meAng,meAccept,m;
  float dAngHB = acut/nAng;
  float donor[3],hydrogen1[3],hydrogen2[3],acceptor[3];
  float yes,angle;
  float r_dist_O;
                                                                                                             
  for(i=0;i<nWt;i++)
  {
    fillDonorInformation(donor,hydrogen1,hydrogen2,posO,posH1,posH2,i);
 
    r_dist_O = GetDisOSOW(posO, i, posOS, nOS, 0, L);

    // bins water
    // by position from OS
    me = WhichBin(rbin, nbin, r_dist_O);
    densWt->m[me][0] += 1.0;
    for(m=0; m<list->m[i][0]; m++)
    {
      j = list->m[i][m+1];
      fillAcceptorInformation(acceptor, posO, j);

      // first O-H1 donor pair
      yes = isHbond(donor,hydrogen1,acceptor,rcut,acut,L);
      if(yes!=2003)
      {	      
        // bins donor 
        // by donor position along radial axis
        hBond->m[me][0]+=1.0;

        // bins OH--O angle or Hbond
        // by donor position along radial aixs
        // by OH--O angle
        meAng = yes/dAngHB;
        OH_Oangle->m[me][meAng]+=1.0;
      
        // bins acceptor
        // by acceptor position radial axis
	r_dist_O = GetDisOSOW(posO, j, posOS, nOS, 0, L);
	meAccept = WhichBin(rbin, nbin, r_dist_O);
        hBond2->m[meAccept][0]+=1.0;		
      }
  
      // second O-H2 donor pair
      yes = isHbond(donor,hydrogen2,acceptor,rcut,acut,L);
      if(yes!=2003)
      {
        hBond->m[me][0]+=1.0;
        meAng = (int)(yes/dAngHB);
        OH_Oangle->m[me][meAng]+=1.0;

	r_dist_O = GetDisOSOW(posO, j, posOS, nOS, 0, L);
	meAccept = WhichBin(rbin, nbin, r_dist_O);
        hBond2->m[meAccept][0]+=1.0;		
      }
    }
  }

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

int Output(MMatrix_t rbin, MMatrix_t densW, MMatrix_t hBond, MMatrix_t hBond2, MMatrix_t OH_Oangle,
           int nbin, int nWt, char *fn, int nAng)
{
  FILE *fp;
  int i,j,k;
  
  if(nWt>0)
  {
    fp = fopen(fn, "w");
    for(i=0;i<nbin;i++)
    {
      for(j=0; j<nAng; j++)
      {
	fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f\n",
	        rbin->m[i][0],
		densW->m[i][0],
		hBond->m[i][0],
		hBond2->m[i][0],
		OH_Oangle->m[i][j]);
      }
    }

    fclose(fp);
  }

  return 0;
}
