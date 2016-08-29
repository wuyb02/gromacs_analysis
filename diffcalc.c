/* 
   Note that we exchange the x- and z- position!

   This code 
   1. read in the trajectory of water, Ion1 and Ion2;
   2. compute the mean square displacement;
   3. output the m.s.d to three different files

   Later Matlab will be used to analyze the diffusivity 
   based on information obtained here.


   Requirements:
   1. Trajectory data must be in format of Na->Cl->SOL, 
      no wall atoms allowed.


   Rui Qiao
   
   May, 10, 2002
   Jan, 01, 2003
*/

#include "string.h"
#include "stdlib.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#define   bsz  1001
#define   bszz 1002

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "when writing your own analysis tools. The advantage of ",
    "using gromacs for this is that you have access to all ",
    "information in the topology, and your program will be ",
    "able to handle all types of coordinates and trajectory ",
    "files supported by gromacs. Go ahead and try it! ",
    "This test version just writes the coordinates of an ",
    "arbitrary atom to standard out for each frame. You can ",
    "select which atom you want to examine with the -n argument."
  };
  
  FILE   *anaResults,*dipoleResults;

  static int   nWt          = 0.0;
  static int   nNa          = 0.0;
  static int   nCl          = 0.0;
  static float Lx           = 4.0;
  static float Ly           = 4.0;
  static float Lz           = 7.7;
  static float Lz1          = 1.0;
  static float Lz2          = 4.4;
  static float R            = 1.0;
  static float Cx           = 3.25;
  static float Cy           = 3.35;
  static float timeInterval = 0.04;
  static float chunk        = 50;
  static int   nbin         = 1;
  static int   totalStep    = 0;

  static int bin1[bsz],bin2[bsz],bin3[bsz];
  static float rbin[bszz],rSpacing;
  static float vbin1[bsz][3],vbin2[bsz][3],vbin3[bsz][3],z,vx,vy,vz;
  static int i,j,nn;
  int sectionSize,effSectSize,totalSection,origins,nmax,section;
  float **posWt,**posNa,**posCl,***msdWt,***msdNa,***msdCl,**msdWtTemp,**msdNaTemp,**msdClTemp;

  int reset1(float ***msdWt,float ***msdNa,float ***msdCl,int n1,int n2,int n3);
  int reset2(int bin1[bsz],int bin2[bsz],int bin3[bsz],int nbin);
  float **CreateMatrix(int rows,int cols);
  float ***CreateMatrix2(int level,int rows,int cols);
  int storePosition(t_trxframe fr,float **posWt,float **posNa,float **posCl,int nWt,int nNa,int nCl,int start,float lb, float bsizeX, float bsizeY, float bsizeZ);
  int averaging(int bin1[bsz],int bin2[bsz],int bin3[bsz],int nmax,float ***msdWt,float ***msdNa,float ***msdCl,int nbin);
  int msdComputing(float **Wt,float **Na,float **Cl,int bin1[bsz],int bin2[bsz],int bin3[bsz],
		    int nWt,int nNa,int nCl,float ***msdWt,float ***msdNa,float ***msdCl,
                    float **msdWtTemp,float **msdNaTemp,float **msdClTemp,
		    int origins,float rSpacing,int nmax,int nbin,float Lx,float Ly,float Lz,float Cx,float Cy,float R,float Lz1,float Lz2);
  int output(float ***msdWt,float ***msdNa,float ***msdCl,float ybin[bszz],float dt,int origins,int nbin,int nNa,int nCl,int nWt);

  t_pargs pa[] = {
    { "-wt",FALSE, etINT, {&nWt},
      "Number of water molecular"
    },
    { "-Na",FALSE, etINT, {&nNa},
      "Number of Na+ ions"
    },
    { "-Cl",FALSE, etINT, {&nCl},
      "Number of Cl- ions"
    },
    { "-lx",FALSE, etREAL, {&Lx},
      "x-Length of simulation box (nm)"
    },
    { "-ly",FALSE, etREAL, {&Ly},
      "y-Length of Simulation box (nm)"
    },
    { "-lz",FALSE, etREAL, {&Lz},
      "z-Length of simulation box (nm)"
    },
    { "-Lz1",FALSE, etREAL, {&Lz1},
      "Lower limit of diffusion domain (nm)"
    },
    { "-Lz2",FALSE, etREAL, {&Lz2},
      "Upper limit of diffusion domain (nm)"
    },
    { "-R",FALSE, etREAL, {&R},
      "Radius of pore (nm)"
    },
    { "-Cx",FALSE, etREAL, {&Cx},
      "X-coordinate of pore center (nm)"
    },
    { "-Cy",FALSE, etREAL, {&Cy},
      "Y-coordinate of pore center (nm)"
    },
    { "-dT",FALSE, etREAL, {&timeInterval},
      "Time interval between trjactory frames (ps)"
    },
    { "-chunk",FALSE, etREAL, {&chunk},
      "Total length of each chunk of trajectory considered (ps)"
    },
    { "-bin",FALSE, etINT, {&nbin},
      "Number of bin across channel"
    },
    { "-snp",FALSE, etINT, {&totalStep},
      "Number of trajectory set"
    }
  };
  
  t_topology top;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X;   /* read only position */
  
  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
 
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);
  sfree(xtop);

  sectionSize  = chunk/timeInterval;
  effSectSize  = sectionSize;       
  totalSection = totalStep/sectionSize;
  origins      = effSectSize/2;        
  nmax         = origins;              

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);
  printf("# of origins : %d\n",origins);

  posWt     = CreateMatrix(effSectSize*nWt+1,3);
  posNa     = CreateMatrix(effSectSize*nNa+1,3);
  posCl     = CreateMatrix(effSectSize*nCl+1,3);
  msdWtTemp = CreateMatrix(3,origins);
  msdNaTemp = CreateMatrix(3,origins);
  msdClTemp = CreateMatrix(3,origins);
  msdWt     = CreateMatrix2(nbin,3,origins);
  msdNa     = CreateMatrix2(nbin,3,origins);
  msdCl     = CreateMatrix2(nbin,3,origins);

  rSpacing = R/nbin;
  for(i=0;i<nbin+1;i++) rbin[i] = i*rSpacing;
  reset1(msdWt,msdNa,msdCl,nbin,3,origins);
  reset2(bin1,bin2,bin3,nbin);

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  storePosition(fr,posWt,posNa,posCl,nWt,nNa,nCl,0,0, Lx, Ly, Lz);
  for(section=0;section<totalSection;section++) 
  {
    for(i=0;i<sectionSize;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize) storePosition(fr,posWt,posNa,posCl,nWt,nNa,nCl,i,0, Lx, Ly, Lz);
    }
    printf("Section %d of total section %d loaded\n",section,totalSection);
    msdComputing(posWt,posNa,posCl,bin1,bin2,bin3,nWt,nNa,nCl,msdWt,msdNa,msdCl,msdWtTemp,msdNaTemp,msdClTemp,
      origins,rSpacing,nmax,nbin,Lx,Ly,Lz,Cx,Cy,R,Lz1,Lz2);
  }
  averaging(bin1,bin2,bin3,nmax,msdWt,msdNa,msdCl,nbin);

  printf("Analysis done\n");
  printf("nbin = %d\n", nbin);
  output(msdWt,msdNa,msdCl,rbin,timeInterval,origins,nbin,nNa,nCl,nWt);
  thanx(stderr);
  
  return 0;
}


int output(float ***msdWt,float ***msdNa,float ***msdCl,float rbin[bszz],float dt,int origins,int nbin,int nNa,int nCl,int nWt)
{
  FILE *fp;
  int i,j,k;

  if(nWt>0)
  {
    fp = fopen("msdWtPore.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.6e %15.6e %15.6e %15.6e\n", j*dt, msdWt[i][0][j], msdWt[i][1][j], msdWt[i][2][j]);
///
/*        fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",(rbin[i]+rbin[i+1])/2,j*dt,msdWt[i][0][j]+msdWt[i][1][j],msdWt[i][2][j],
          msdWt[i][0][j]+msdWt[i][1][j]+msdWt[i][2][j]);
*/
    fclose(fp);
  }
  
  if(nNa>0)
  {
    fp = fopen("msdNaPore.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",(rbin[i]+rbin[i+1])/2,j*dt,msdNa[i][0][j]+msdNa[i][1][j],msdNa[i][2][j],
          msdNa[i][0][j]+msdNa[i][1][j]+msdNa[i][2][j]);
    fclose(fp);
  }

  if(nCl>0)
  {
    fp = fopen("msdClPore.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",(rbin[i]+rbin[i+1])/2,j*dt,msdCl[i][0][j]+msdCl[i][1][j],msdCl[i][2][j],
          msdCl[i][0][j]+msdCl[i][1][j]+msdCl[i][2][j]);
    fclose(fp);
  }

  return 0;
}


int storePosition(t_trxframe fr,float **posWt,float **posNa,float **posCl,int nWt,int nNa,int nCl,int start,float lb, float bsizeX, float bsizeY, float bsizeZ)
{
  int i, j, k;
  int ll, mm, nn;
  int begin,nOther;
  float zO,zH1,zH2,xO,xH1,xH2,yO,yH1,yH2,xxx,yyy,zzz;
  float bsize[3] = {bsizeX, bsizeY, bsizeZ};
  for(i=0; i<3; i++)
  {
    if(bsize[i]<3.0)
    {
      printf("bsize[%d] too small\n", i);
    }
  }
  
  begin = start*nWt;
  nOther = 0;
//  nOther = 318;

  for(i=0;i<nWt;i++)
  {
    //Remove PBC
    for(j=1; j<64; j++)
    {
      ll = nOther+i*64+j;
      for(k=0; k<3; k++)
      {
	printf("%f, ", fr.x[ll][k]);
        if(fabs(fr.x[ll][k]-fr.x[ll-1][k])>bsize[k]/4*3)
        {
          if(fr.x[ll][k]>fr.x[ll-1][k])
          {
            fr.x[ll][k] -= bsize[k];
          }else
          {
            fr.x[ll][k] += bsize[k];
          }
	  printf("*");
        }
	printf("%f      ", fr.x[ll][k]);
      }
      printf("\n");
    }
    printf("*******endofFrame********\n\n\n\n\n");

///
/*
      posWt[begin+i][0] = fr.x[nOther+i*3][XX]; 
      posWt[begin+i][1] = fr.x[nOther+i*3][YY];
      posWt[begin+i][2] = fr.x[nOther+i*3][ZZ]-lb;
*/
	posWt[begin+i][0] = 1	/	443.5676	*
			(
				fr.x[nOther+0+i*64][0]	*	 12.0110	+
				fr.x[nOther+1+i*64][0]	*	  1.0080	+
				fr.x[nOther+2+i*64][0]	*	  1.0080	+
				fr.x[nOther+3+i*64][0]	*	  1.0080	+
				fr.x[nOther+4+i*64][0]	*	 12.0110	+
				fr.x[nOther+5+i*64][0]	*	  1.0080	+
				fr.x[nOther+6+i*64][0]	*	  1.0080	+
				fr.x[nOther+7+i*64][0]	*	 15.9994	+
				fr.x[nOther+8+i*64][0]	*	 12.0110	+
				fr.x[nOther+9+i*64][0]	*	 15.9994	+
				fr.x[nOther+10+i*64][0]	*	 12.0110	+
				fr.x[nOther+11+i*64][0]	*	 12.0110	+
				fr.x[nOther+12+i*64][0]	*	  1.0080	+
				fr.x[nOther+13+i*64][0]	*	 12.0110	+
				fr.x[nOther+14+i*64][0]	*	  1.0080	+
				fr.x[nOther+15+i*64][0]	*	 12.0110	+
				fr.x[nOther+16+i*64][0]	*	  1.0080	+
				fr.x[nOther+17+i*64][0]	*	 12.0110	+
				fr.x[nOther+18+i*64][0]	*	  1.0080	+
				fr.x[nOther+19+i*64][0]	*	 12.0110	+
				fr.x[nOther+20+i*64][0]	*	 12.0110	+
				fr.x[nOther+21+i*64][0]	*	 12.0110	+
				fr.x[nOther+22+i*64][0]	*	 12.0110	+
				fr.x[nOther+23+i*64][0]	*	  1.0080	+
				fr.x[nOther+24+i*64][0]	*	 12.0110	+
				fr.x[nOther+25+i*64][0]	*	 12.0110	+
				fr.x[nOther+26+i*64][0]	*	  1.0080	+
				fr.x[nOther+27+i*64][0]	*	  1.0080	+
				fr.x[nOther+28+i*64][0]	*	  1.0080	+
				fr.x[nOther+29+i*64][0]	*	 12.0110	+
				fr.x[nOther+30+i*64][0]	*	 14.0067	+
				fr.x[nOther+31+i*64][0]	*	  1.0080	+
				fr.x[nOther+32+i*64][0]	*	 12.0110	+
				fr.x[nOther+33+i*64][0]	*	  1.0080	+
				fr.x[nOther+34+i*64][0]	*	  1.0080	+
				fr.x[nOther+35+i*64][0]	*	 12.0110	+
				fr.x[nOther+36+i*64][0]	*	  1.0080	+
				fr.x[nOther+37+i*64][0]	*	  1.0080	+
				fr.x[nOther+38+i*64][0]	*	  1.0080	+
				fr.x[nOther+39+i*64][0]	*	 12.0110	+
				fr.x[nOther+40+i*64][0]	*	  1.0080	+
				fr.x[nOther+41+i*64][0]	*	 12.0110	+
				fr.x[nOther+42+i*64][0]	*	 15.9994	+
				fr.x[nOther+43+i*64][0]	*	 12.0110	+
				fr.x[nOther+44+i*64][0]	*	 12.0110	+
				fr.x[nOther+45+i*64][0]	*	  1.0080	+
				fr.x[nOther+46+i*64][0]	*	 12.0110	+
				fr.x[nOther+47+i*64][0]	*	 12.0110	+
				fr.x[nOther+48+i*64][0]	*	  1.0080	+
				fr.x[nOther+49+i*64][0]	*	 12.0110	+
				fr.x[nOther+50+i*64][0]	*	 12.0110	+
				fr.x[nOther+51+i*64][0]	*	  1.0080	+
				fr.x[nOther+52+i*64][0]	*	  1.0080	+
				fr.x[nOther+53+i*64][0]	*	  1.0080	+
				fr.x[nOther+54+i*64][0]	*	 12.0110	+
				fr.x[nOther+55+i*64][0]	*	 14.0067	+
				fr.x[nOther+56+i*64][0]	*	  1.0080	+
				fr.x[nOther+57+i*64][0]	*	 12.0110	+
				fr.x[nOther+58+i*64][0]	*	  1.0080	+
				fr.x[nOther+59+i*64][0]	*	  1.0080	+
				fr.x[nOther+60+i*64][0]	*	 12.0110	+
				fr.x[nOther+61+i*64][0]	*	  1.0080	+
				fr.x[nOther+62+i*64][0]	*	  1.0080	+
				fr.x[nOther+63+i*64][0]	*	  1.0080
			);
	posWt[begin+i][1] = 1	/	443.5676	*
			(
				fr.x[nOther+0+i*64][1]	*	 12.0110	+
				fr.x[nOther+1+i*64][1]	*	  1.0080	+
				fr.x[nOther+2+i*64][1]	*	  1.0080	+
				fr.x[nOther+3+i*64][1]	*	  1.0080	+
				fr.x[nOther+4+i*64][1]	*	 12.0110	+
				fr.x[nOther+5+i*64][1]	*	  1.0080	+
				fr.x[nOther+6+i*64][1]	*	  1.0080	+
				fr.x[nOther+7+i*64][1]	*	 15.9994	+
				fr.x[nOther+8+i*64][1]	*	 12.0110	+
				fr.x[nOther+9+i*64][1]	*	 15.9994	+
				fr.x[nOther+10+i*64][1]	*	 12.0110	+
				fr.x[nOther+11+i*64][1]	*	 12.0110	+
				fr.x[nOther+12+i*64][1]	*	  1.0080	+
				fr.x[nOther+13+i*64][1]	*	 12.0110	+
				fr.x[nOther+14+i*64][1]	*	  1.0080	+
				fr.x[nOther+15+i*64][1]	*	 12.0110	+
				fr.x[nOther+16+i*64][1]	*	  1.0080	+
				fr.x[nOther+17+i*64][1]	*	 12.0110	+
				fr.x[nOther+18+i*64][1]	*	  1.0080	+
				fr.x[nOther+19+i*64][1]	*	 12.0110	+
				fr.x[nOther+20+i*64][1]	*	 12.0110	+
				fr.x[nOther+21+i*64][1]	*	 12.0110	+
				fr.x[nOther+22+i*64][1]	*	 12.0110	+
				fr.x[nOther+23+i*64][1]	*	  1.0080	+
				fr.x[nOther+24+i*64][1]	*	 12.0110	+
				fr.x[nOther+25+i*64][1]	*	 12.0110	+
				fr.x[nOther+26+i*64][1]	*	  1.0080	+
				fr.x[nOther+27+i*64][1]	*	  1.0080	+
				fr.x[nOther+28+i*64][1]	*	  1.0080	+
				fr.x[nOther+29+i*64][1]	*	 12.0110	+
				fr.x[nOther+30+i*64][1]	*	 14.0067	+
				fr.x[nOther+31+i*64][1]	*	  1.0080	+
				fr.x[nOther+32+i*64][1]	*	 12.0110	+
				fr.x[nOther+33+i*64][1]	*	  1.0080	+
				fr.x[nOther+34+i*64][1]	*	  1.0080	+
				fr.x[nOther+35+i*64][1]	*	 12.0110	+
				fr.x[nOther+36+i*64][1]	*	  1.0080	+
				fr.x[nOther+37+i*64][1]	*	  1.0080	+
				fr.x[nOther+38+i*64][1]	*	  1.0080	+
				fr.x[nOther+39+i*64][1]	*	 12.0110	+
				fr.x[nOther+40+i*64][1]	*	  1.0080	+
				fr.x[nOther+41+i*64][1]	*	 12.0110	+
				fr.x[nOther+42+i*64][1]	*	 15.9994	+
				fr.x[nOther+43+i*64][1]	*	 12.0110	+
				fr.x[nOther+44+i*64][1]	*	 12.0110	+
				fr.x[nOther+45+i*64][1]	*	  1.0080	+
				fr.x[nOther+46+i*64][1]	*	 12.0110	+
				fr.x[nOther+47+i*64][1]	*	 12.0110	+
				fr.x[nOther+48+i*64][1]	*	  1.0080	+
				fr.x[nOther+49+i*64][1]	*	 12.0110	+
				fr.x[nOther+50+i*64][1]	*	 12.0110	+
				fr.x[nOther+51+i*64][1]	*	  1.0080	+
				fr.x[nOther+52+i*64][1]	*	  1.0080	+
				fr.x[nOther+53+i*64][1]	*	  1.0080	+
				fr.x[nOther+54+i*64][1]	*	 12.0110	+
				fr.x[nOther+55+i*64][1]	*	 14.0067	+
				fr.x[nOther+56+i*64][1]	*	  1.0080	+
				fr.x[nOther+57+i*64][1]	*	 12.0110	+
				fr.x[nOther+58+i*64][1]	*	  1.0080	+
				fr.x[nOther+59+i*64][1]	*	  1.0080	+
				fr.x[nOther+60+i*64][1]	*	 12.0110	+
				fr.x[nOther+61+i*64][1]	*	  1.0080	+
				fr.x[nOther+62+i*64][1]	*	  1.0080	+
				fr.x[nOther+63+i*64][1]	*	  1.0080
			);
	posWt[begin+i][2] = 1	/	443.5676	*
			(
				fr.x[nOther+0+i*64][2]	*	 12.0110	+
				fr.x[nOther+1+i*64][2]	*	  1.0080	+
				fr.x[nOther+2+i*64][2]	*	  1.0080	+
				fr.x[nOther+3+i*64][2]	*	  1.0080	+
				fr.x[nOther+4+i*64][2]	*	 12.0110	+
				fr.x[nOther+5+i*64][2]	*	  1.0080	+
				fr.x[nOther+6+i*64][2]	*	  1.0080	+
				fr.x[nOther+7+i*64][2]	*	 15.9994	+
				fr.x[nOther+8+i*64][2]	*	 12.0110	+
				fr.x[nOther+9+i*64][2]	*	 15.9994	+
				fr.x[nOther+10+i*64][2]	*	 12.0110	+
				fr.x[nOther+11+i*64][2]	*	 12.0110	+
				fr.x[nOther+12+i*64][2]	*	  1.0080	+
				fr.x[nOther+13+i*64][2]	*	 12.0110	+
				fr.x[nOther+14+i*64][2]	*	  1.0080	+
				fr.x[nOther+15+i*64][2]	*	 12.0110	+
				fr.x[nOther+16+i*64][2]	*	  1.0080	+
				fr.x[nOther+17+i*64][2]	*	 12.0110	+
				fr.x[nOther+18+i*64][2]	*	  1.0080	+
				fr.x[nOther+19+i*64][2]	*	 12.0110	+
				fr.x[nOther+20+i*64][2]	*	 12.0110	+
				fr.x[nOther+21+i*64][2]	*	 12.0110	+
				fr.x[nOther+22+i*64][2]	*	 12.0110	+
				fr.x[nOther+23+i*64][2]	*	  1.0080	+
				fr.x[nOther+24+i*64][2]	*	 12.0110	+
				fr.x[nOther+25+i*64][2]	*	 12.0110	+
				fr.x[nOther+26+i*64][2]	*	  1.0080	+
				fr.x[nOther+27+i*64][2]	*	  1.0080	+
				fr.x[nOther+28+i*64][2]	*	  1.0080	+
				fr.x[nOther+29+i*64][2]	*	 12.0110	+
				fr.x[nOther+30+i*64][2]	*	 14.0067	+
				fr.x[nOther+31+i*64][2]	*	  1.0080	+
				fr.x[nOther+32+i*64][2]	*	 12.0110	+
				fr.x[nOther+33+i*64][2]	*	  1.0080	+
				fr.x[nOther+34+i*64][2]	*	  1.0080	+
				fr.x[nOther+35+i*64][2]	*	 12.0110	+
				fr.x[nOther+36+i*64][2]	*	  1.0080	+
				fr.x[nOther+37+i*64][2]	*	  1.0080	+
				fr.x[nOther+38+i*64][2]	*	  1.0080	+
				fr.x[nOther+39+i*64][2]	*	 12.0110	+
				fr.x[nOther+40+i*64][2]	*	  1.0080	+
				fr.x[nOther+41+i*64][2]	*	 12.0110	+
				fr.x[nOther+42+i*64][2]	*	 15.9994	+
				fr.x[nOther+43+i*64][2]	*	 12.0110	+
				fr.x[nOther+44+i*64][2]	*	 12.0110	+
				fr.x[nOther+45+i*64][2]	*	  1.0080	+
				fr.x[nOther+46+i*64][2]	*	 12.0110	+
				fr.x[nOther+47+i*64][2]	*	 12.0110	+
				fr.x[nOther+48+i*64][2]	*	  1.0080	+
				fr.x[nOther+49+i*64][2]	*	 12.0110	+
				fr.x[nOther+50+i*64][2]	*	 12.0110	+
				fr.x[nOther+51+i*64][2]	*	  1.0080	+
				fr.x[nOther+52+i*64][2]	*	  1.0080	+
				fr.x[nOther+53+i*64][2]	*	  1.0080	+
				fr.x[nOther+54+i*64][2]	*	 12.0110	+
				fr.x[nOther+55+i*64][2]	*	 14.0067	+
				fr.x[nOther+56+i*64][2]	*	  1.0080	+
				fr.x[nOther+57+i*64][2]	*	 12.0110	+
				fr.x[nOther+58+i*64][2]	*	  1.0080	+
				fr.x[nOther+59+i*64][2]	*	  1.0080	+
				fr.x[nOther+60+i*64][2]	*	 12.0110	+
				fr.x[nOther+61+i*64][2]	*	  1.0080	+
				fr.x[nOther+62+i*64][2]	*	  1.0080	+
				fr.x[nOther+63+i*64][2]	*	  1.0080
			);
  }
  
  begin = start*nNa;
  nOther = nWt*3;
  for(i=0;i<nNa;i++)
    {
      posNa[begin+i][0] = fr.x[nOther+i][XX];
      posNa[begin+i][1] = fr.x[nOther+i][YY];
      posNa[begin+i][2] = fr.x[nOther+i][ZZ]-lb;
    }
  
  begin = start*nCl;
  nOther = nWt*3+nNa;
  for(i=0;i<nCl;i++)
    {
      posCl[begin+i][0] = fr.x[nOther+i][XX];
      posCl[begin+i][1] = fr.x[nOther+i][YY];
      posCl[begin+i][2] = fr.x[nOther+i][ZZ]-lb;
    }
  return 0;
}


float periodicity(float seem,float box)
{
  float actual;
  if(fabs(seem)>box/2)
    actual = box-fabs(seem);
  else
    actual = seem;
  return actual;
}


int averaging(int bin1[bsz],int bin2[bsz],int bin3[bsz],int nmax,float ***msdWt,float ***msdNa,float ***msdCl,int nbin)
{
  int i,j,k;

  for(i=0;i<nbin;i++)
    for(j=0;j<3;j++)
      for(k=0;k<nmax;k++)
	if(bin1[i]!=0)
	  msdWt[i][j][k] /= (bin1[i]*nmax);

  for(i=0;i<nbin;i++)
    for(j=0;j<3;j++)
      for(k=0;k<nmax;k++)
	if(bin2[i]!=0)
	  msdNa[i][j][k] /= (bin2[i]*nmax);
  
  for(i=0;i<nbin;i++)
    for(j=0;j<3;j++)
      for(k=0;k<nmax;k++)
	if(bin3[i]!=0)
	  msdCl[i][j][k] /= (bin3[i]*nmax);
  
  return 0;
}


int msdComputing(float **Wt,float **Na,float **Cl,int bin1[bsz],int bin2[bsz],int bin3[bsz],
                 int nWt,int nNa,int nCl,float ***msdWt,float ***msdNa,float ***msdCl,
                 float **msdWtTemp,float **msdNaTemp,float **msdClTemp,
                 int origins,float rSpacing,int nmax,int nbin,float Lx,float Ly,float Lz,float Cx,float Cy,float R,float Lz1,float Lz2)
{
  int i,j,k,me,l,m,n,indicator;
  int jstart,kend;

  if(nWt != 0)
  {
    for(i=0;i<nWt;i++)
    {
      if(sqrt(pow(Wt[i][0]-Cx,2)+(pow(Wt[i][1]-Cy,2)))<R)
      {
        indicator = 1;
        me = sqrt(pow(Wt[i][0]-Cx,2)+(pow(Wt[i][1]-Cy,2)))/rSpacing;
        if(me > nbin) me = nbin;

        for(l=0;l<3;l++)
        for(m=0;m<nmax;m++)
          msdWtTemp[l][m] = 0;

        for(j=0;j<origins;j++)
        {
          jstart = j*nWt+i;
          for(k=0;k<nmax;k++)
 	  {
	    kend = jstart+k*nWt;
            //if(Wt[kend][2]<Lz1 || Wt[kend][2]>Lz2) indicator = 0;
	    msdWtTemp[0][k] += pow(periodicity(Wt[kend][0]-Wt[jstart][0],Lx),2);
	    msdWtTemp[1][k] += pow(periodicity(Wt[kend][1]-Wt[jstart][1],Ly),2);
            msdWtTemp[2][k] += pow(periodicity(Wt[kend][2]-Wt[jstart][2],Lz),2);
	  }
        }

        if (indicator == 1)
        {
	  bin1[me]++;
          for(k=0;k<nmax;k++)
 	  {
             msdWt[me][0][k] += msdWtTemp[0][k];
	     msdWt[me][1][k] += msdWtTemp[1][k];
	     msdWt[me][2][k] += msdWtTemp[2][k];
          }
        }
      }
    }
  }
    
  if(nNa != 0)
  {
    for(i=0;i<nNa;i++)
    { 
      if(sqrt(pow(Na[i][0]-Cx,2)+(pow(Na[i][1]-Cy,2)))<R)
      {
        indicator = 1;
        me = sqrt(pow(Na[i][0]-Cx,2)+(pow(Na[i][1]-Cy,2)))/rSpacing;
        if(me > nbin) me = nbin;

        for(l=0;l<3;l++)
        for(m=0;m<nmax;m++)
          msdNaTemp[l][m] = 0;

        for(j=0;j<origins;j++)
        {
          jstart = j*nNa+i;
          for(k=0;k<nmax;k++)
 	  {
	    kend = jstart+k*nNa;
            if(Na[kend][2]<Lz1 || Na[kend][2]>Lz2) indicator = 0;
	    msdNaTemp[0][k] += pow(periodicity(Na[kend][0]-Na[jstart][0],Lx),2);
	    msdNaTemp[1][k] += pow(periodicity(Na[kend][1]-Na[jstart][1],Ly),2);
            msdNaTemp[2][k] += pow(periodicity(Na[kend][2]-Na[jstart][2],Lz),2);
	  }
        }

        if (indicator == 1)
        {
	  bin2[me]++;
          for(k=0;k<nmax;k++)
 	  {
            msdNa[me][0][k] += msdNaTemp[0][k];
	    msdNa[me][1][k] += msdNaTemp[1][k];
	    msdNa[me][2][k] += msdNaTemp[2][k];
          }
        }
      }
    }
  }
 
  if(nCl != 0)
  {
    for(i=0;i<nCl;i++)
    {
      if(sqrt(pow(Cl[i][0]-Cx,2)+(pow(Cl[i][1]-Cy,2)))<R)
      {
        indicator = 1; 
        me = sqrt(pow(Cl[i][0]-Cx,2)+(pow(Cl[i][1]-Cy,2)))/rSpacing;
        if (me > nbin) me = nbin;

        for(l=0;l<3;l++)
        for(m=0;m<nmax;m++)
          msdClTemp[l][m] = 0;

        for(j=0;j<origins;j++)
        {
          jstart = j*nCl+i;
          for(k=0;k<nmax;k++)
  	  {
 	    kend = jstart+k*nCl;
            if(Cl[kend][2]<Lz1 || Cl[kend][2]>Lz2) indicator = 0;
	    msdClTemp[0][k] += pow(periodicity(Cl[kend][0]-Cl[jstart][0],Lx),2);
	    msdClTemp[1][k] += pow(periodicity(Cl[kend][1]-Cl[jstart][1],Ly),2);
            msdClTemp[2][k] += pow(periodicity(Cl[kend][2]-Cl[jstart][2],Lz),2);
          }
        }

        if (indicator == 1)
        {
	  bin3[me]++;
          for(k=0;k<nmax;k++)
 	  {
            msdCl[me][0][k] += msdClTemp[0][k];
	    msdCl[me][1][k] += msdClTemp[1][k];
	    msdCl[me][2][k] += msdClTemp[2][k];
          }
        }
      }
    }
  }

  return 0;
}

 
int reset1(float ***msdWt,float ***msdNa,float ***msdCl,int n1,int n2,int n3)
{
  int i,j,k;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	msdWt[i][j][k] = msdNa[i][j][k] = msdCl[i][j][k] = 0;
  
  return 0;
}


int reset2(int *bin1,int *bin2,int *bin3,int nbin)
{
  int i;
  for(i=0;i<nbin;i++)
    bin1[i] = bin2[i] = bin3[i] = 0;
  
  return 0;
}


float ** CreateMatrix(int rows,int cols)
{
   int	i;
   float  **m;

   m = calloc((unsigned int) rows,sizeof(float *));
   for (i=0; i < rows; i++) {
      m[i] = calloc((unsigned int) cols,sizeof(float ));
   }

   return m;
}


float *** CreateMatrix2(int level,int rows,int cols)
{
   int	i;
   float  ***m;

   m = calloc((unsigned int) level,sizeof(float **));
   for (i=0; i < level; i++) 
     m[i] = CreateMatrix(rows,cols);
   
   return m;
}
