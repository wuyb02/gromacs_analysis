/* 
   This code calculates
     1. time-autocorrelations of end-to-end vectors
     2. time-autocorrelations of principal-axis vectors
        - Mondello and Grest, JCP, 103, 7156 (1995)
        - Takeuchi, JCP, 109, 5614 (1998)
     3. time-autocorrelations of segment vectors
        - Takeuchi and Roe, JCP, 94, 7446 (1991)
        - Mondello and Grest, JCP, 103, 7156 (1995)
     4. time-autocorrelations of dihedral angles
        - Takeuchi and Roe, JCP, 94, 7446 (1991)
        - Brown et al, JCP, 100, 1684 (1994)
        - Dysthe et al, JCP, 110, 4047 (1999)

   - Programmed by Jae H. Park
   - Spetember 8, 2006  

*/

#include "string.h"
#include "stdlib.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "complex.h"
#define   bsz  1001

float *CreateVector(int rows);
float **CreateMatrix(int rows,int cols);
float periodicity(float seem,float box);
int storePositionOM(t_trxframe fr,float **posOM,float **vecEE,float **vecPA,int nOM,int nAT,int start);
int storePositionSG(t_trxframe fr,float **posSG,float **vecSGa,float **vecSGb,float **vecSGc,int nOM,int nAT,int nSG,int start);
int storePositionDA(t_trxframe fr,float **posDA,float *phiDA,int nOM,int nAT,int nDA,int start);
int autoCorrelationOM(float **OM,float **vecEE,float **vecPA,int binOM[bsz],int nOM,float **corEE01,float **corEE02,float **corPA01,float **corPA02,int origins,int nmax,float z1,float z2,float zsize);
int autoCorrelationSG(float **SG,float **vecSGa,float **vecSGb,float **vecSGc,int binSG[bsz],int nSG,float **corSG01a,float **corSG02a,float **corSG01b,float **corSG02b,float **corSG01c,float **corSG02c,int origins,int nmax,float z1,float z2, float zsize);
int autoCorrelationDA(float **DA,float *phiDA,int binDA[bsz],int nDA,float **corDA01,float **corDA02,float **corDA03,int origins,int nmax,float z1,float z2, float zsize);


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
  
  static int   nOM       = 0;
  static int   nAT       = 0;
  static float z1        = 0.0;
  static float z2        = 0.0;
  static float deltaT    = 0.0;
  static float chunk     = 0;
  static int   nbin      = 0;
  static int   totalStep = 0;

  FILE *fp;
  int i,j,k,binOM[bsz],binSG[bsz],binDA[bsz];
  int sectionSize,effSectSize,totalSection,origins,nmax,section;
  float zsize,*zbin,**posOM,**vecEE,**vecPA,**corEE01,**corEE02,**corPA01,**corPA02;

  int nSG;
  float **posSG,**vecSGa,**vecSGb,**vecSGc,**corSG01a,**corSG01b,**corSG01c,**corSG02a,**corSG02b,**corSG02c;

  int nDA;
  float **posDA,*phiDA,**corDA01,**corDA02,**corDA03,**corDA;

  t_pargs pa[] = {
    { "-nOM",FALSE, etINT, {&nOM},
      "Number of organic molecules"
    },
    { "-nAT",FALSE, etINT, {&nAT},
      "Number of atoms in single organic molecules"
    },
    { "-z1",FALSE, etREAL, {&z1},
      "lower limit for compuation (nm)"
    },
    { "-z2",FALSE, etREAL, {&z2},
      "upper limit for computation (nm)"
    },
    { "-dT",FALSE, etREAL, {&deltaT},
      "Time interval between trjactory frames (ps)"
    },
    { "-chunk",FALSE, etREAL, {&chunk},
      "Length of chunk in trajectory considered (ps)"
    },
    { "-nbin",FALSE, etINT, {&nbin},
      "Number of bins in sampling domain (nm)"
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

  sectionSize  = chunk/deltaT;
  effSectSize  = sectionSize;       
  totalSection = totalStep/sectionSize;
  origins      = effSectSize/2;        
  nmax         = origins;              

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);
  printf("# of origins : %d\n",origins);

  nSG = nOM*(nAT-2);
  nDA = nOM*(nAT-3);

  zbin      = CreateVector(nbin+1);

  posOM     = CreateMatrix(effSectSize*nOM+1,3);
  vecEE     = CreateMatrix(effSectSize*nOM+1,3);
  vecPA     = CreateMatrix(effSectSize*nOM+1,3);
  corEE01   = CreateMatrix(nbin,origins);
  corEE02   = CreateMatrix(nbin,origins);
  corPA01   = CreateMatrix(nbin,origins);
  corPA02   = CreateMatrix(nbin,origins);

  posSG     = CreateMatrix(effSectSize*nSG+1,3);
  vecSGa    = CreateMatrix(effSectSize*nSG+1,3);
  vecSGb    = CreateMatrix(effSectSize*nSG+1,3);
  vecSGc    = CreateMatrix(effSectSize*nSG+1,3);
  corSG01a  = CreateMatrix(nbin,origins);
  corSG01b  = CreateMatrix(nbin,origins);
  corSG01c  = CreateMatrix(nbin,origins);
  corSG02a  = CreateMatrix(nbin,origins);
  corSG02b  = CreateMatrix(nbin,origins);
  corSG02c  = CreateMatrix(nbin,origins);

  posDA     = CreateMatrix(effSectSize*nDA+1,3);
  phiDA     = CreateVector(effSectSize*nDA+1);
  corDA01   = CreateMatrix(nbin,origins);
  corDA02   = CreateMatrix(nbin,origins);
  corDA03   = CreateMatrix(nbin,origins);
  corDA     = CreateMatrix(nbin,origins);

  zsize = (z2-z1)/nbin;
  for(i=0;i<nbin+1;i++) zbin[i] = zsize*i;
  for(i=0;i<nbin;i++) 
  {
    binOM[i] = 0;
    binSG[i] = 0;
    for(j=0;j<origins;j++){
      corEE01[i][j]  = 0.0;
      corEE02[i][j]  = 0.0;
      corPA01[i][j]  = 0.0;
      corPA02[i][j]  = 0.0;
      corSG01a[i][j] = 0.0;
      corSG01b[i][j] = 0.0;
      corSG01c[i][j] = 0.0;
      corSG02a[i][j] = 0.0;
      corSG02b[i][j] = 0.0;
      corSG02c[i][j] = 0.0;
      corDA01[i][j]  = 0.0;
      corDA02[i][j]  = 0.0;
      corDA03[i][j]  = 0.0;
      corDA[i][j]  = 0.0;
    }
  }

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
  storePositionOM(fr,posOM,vecEE,vecPA,nOM,nAT,0); 
  storePositionSG(fr,posSG,vecSGa,vecSGb,vecSGc,nOM,nAT,nSG,0);
  storePositionDA(fr,posDA,phiDA,nOM,nAT,nDA,0);
  for(section=0;section<totalSection;section++) 
  {
    for(i=0;i<sectionSize;i++)
    {
      read_next_frame(status,&fr);
      if(i<effSectSize)
      { 
        storePositionOM(fr,posOM,vecEE,vecPA,nOM,nAT,i);
        storePositionSG(fr,posSG,vecSGa,vecSGb,vecSGc,nOM,nAT,nSG,i);
        storePositionDA(fr,posDA,phiDA,nOM,nAT,nDA,i);
      }
    }  
    printf("Section %d of total section %d loaded\n",section,totalSection);
    autoCorrelationOM(posOM,vecEE,vecPA,binOM,nOM,corEE01,corEE02,corPA01,corPA02,origins,nmax,z1,z2,zsize);
    autoCorrelationSG(posSG,vecSGa,vecSGb,vecSGc,binSG,nSG,corSG01a,corSG02a,corSG01b,corSG02b,corSG01c,corSG02c,origins,nmax,z1,z2,zsize);
    autoCorrelationDA(posDA,phiDA,binDA,nDA,corDA01,corDA02,corDA03,origins,nmax,z1,z2,zsize);

  }

  printf("%d \n",binOM[0]*nmax);
  for(i=0;i<nbin;i++) for(j=0;j<nmax;j++) if(binOM[i] != 0.0)
  {
    corEE01[i][j] /= binOM[i]*nmax;  
    corEE02[i][j] /= binOM[i]*nmax;  
    corPA01[i][j] /= binOM[i]*nmax;  
    corPA02[i][j] /= binOM[i]*nmax;  
  }

  printf("%d \n",binSG[0]*nmax);
  for(i=0;i<nbin;i++) for(j=0;j<nmax;j++) if(binSG[i] != 0.0)
  {
    corSG01a[i][j] /= binSG[i]*nmax;  
    corSG02a[i][j] /= binSG[i]*nmax;  
    corSG01b[i][j] /= binSG[i]*nmax;  
    corSG02b[i][j] /= binSG[i]*nmax;  
    corSG01c[i][j] /= binSG[i]*nmax;  
    corSG02c[i][j] /= binSG[i]*nmax;  
  }

  printf("%d \n",binDA[0]*nmax);
  for(i=0;i<nbin;i++) for(j=0;j<nmax;j++) if(binDA[i] != 0.0)
  {
    corDA01[i][j] /= binDA[i]*nmax;
    corDA02[i][j] /= binDA[i]*nmax;
    corDA03[i][j] /= binDA[i]*nmax;
    corDA[i][j] = (corDA01[i][j]-pow(corDA02[i][j],2))/(corDA03[i][j]-pow(corDA02[i][j],2));
  }

  printf("Analysis done\n");
  if(nOM>0)
  {
    fp = fopen("angCorrOM.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n",(zbin[i]+zbin[i+1])/2,j*deltaT,corEE01[i][j],corEE02[i][j],corPA01[i][j],corPA02[i][j]);
    fclose(fp);
  }
  
  if(nSG>0)
  {
    fp = fopen("angCorrSG.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n",(zbin[i]+zbin[i+1])/2,j*deltaT,corSG01a[i][j],corSG02a[i][j],corSG01b[i][j],corSG02b[i][j],corSG01c[i][j],corSG02c[i][j]);
    fclose(fp);
  }

  if(nDA>0)
  {
    fp = fopen("angCorrDA.dat","w");
    for(i=0;i<nbin;i++)
      for(j=0;j<origins;j++)
        fprintf(fp,"%15.5e%15.5e%15.5e\n",(zbin[i]+zbin[i+1])/2,j*deltaT,corDA[i][j]);
    fclose(fp);
  }
  

  return 0;
}


int storePositionOM(t_trxframe fr,float **posOM,float **vecEE,float **vecPA,int nOM,int nAT,int start)
{
  int i,j;
  int begin,nOther;
  float temp;
  float **posAT,*xmass;
  float massCH2 = 14.027;
  float massCH3 = 15.035;  
  float totMASS = massCH3*2+massCH2*(nAT-2);
  float vecMAG;
  float xx,yy,zz,inertia[3][3],b,c,d,q,r;
  complex s,t;
  complex x1,x2,x3;
  float eigenvalue[3],eigenmin;
  float bb[3][3];

  posAT = CreateMatrix(nAT,3);
  xmass = CreateVector(nAT);

  xmass[    0] = massCH3;
  for(i=1;i<nAT-1;i++)  xmass[i] = massCH2;
  xmass[nAT-1] = massCH3;
  
  begin = start*nOM;
  nOther = 0;
  for(i=0;i<nOM;i++)
  {
    for(j=0;j<nAT;j++)
    {
      posAT[j][0] = fr.x[nOther+i*nAT+j][XX];
      posAT[j][1] = fr.x[nOther+i*nAT+j][YY];
      posAT[j][2] = fr.x[nOther+i*nAT+j][ZZ];
    }

    posOM[begin+i][0] = 0.0; 
    posOM[begin+i][1] = 0.0; 
    posOM[begin+i][2] = 0.0; 
    for(j=0;j<nAT;j++)
    {
      posOM[begin+i][0] += xmass[j]*posAT[j][0];
      posOM[begin+i][1] += xmass[j]*posAT[j][1];
      posOM[begin+i][2] += xmass[j]*posAT[j][2];
    }
    posOM[begin+i][0] /= totMASS; 
    posOM[begin+i][1] /= totMASS;
    posOM[begin+i][2] /= totMASS;

    vecEE[begin+i][0] = posAT[nAT-1][0]-posAT[0][0]; 
    vecEE[begin+i][1] = posAT[nAT-1][1]-posAT[0][1]; 
    vecEE[begin+i][2] = posAT[nAT-1][2]-posAT[0][2]; 
    vecMAG=sqrt(pow(vecEE[begin+i][0],2)+pow(vecEE[begin+i][1],2)+pow(vecEE[begin+i][2],2));

    vecEE[begin+i][0] /= vecMAG;
    vecEE[begin+i][1] /= vecMAG;
    vecEE[begin+i][2] /= vecMAG;

    inertia[0][0] = 0.0;
    inertia[0][1] = 0.0;
    inertia[0][2] = 0.0;
    inertia[1][0] = 0.0;
    inertia[1][1] = 0.0;
    inertia[1][2] = 0.0;
    inertia[2][0] = 0.0;
    inertia[2][1] = 0.0;
    inertia[2][2] = 0.0;
    for(j=0;j<nAT;j++)
    {
      xx = posAT[j][0]-posOM[begin+i][0];
      yy = posAT[j][1]-posOM[begin+i][1];  
      zz = posAT[j][2]-posOM[begin+i][2];  
      inertia[0][0] += xmass[j]*(pow(yy,2)+pow(zz,2));
      inertia[0][1] -= xmass[j]*xx*yy;
      inertia[0][2] -= xmass[j]*xx*zz;
      inertia[1][1] += xmass[j]*(pow(xx,2)+pow(zz,2));
      inertia[1][2] -= xmass[j]*yy*zz;
      inertia[2][2] += xmass[j]*(pow(xx,2)+pow(yy,2));
    }
    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];

    b = -(inertia[0][0]+inertia[1][1]+inertia[2][2]);
    c = inertia[0][0]*inertia[1][1]+inertia[1][1]*inertia[2][2]+inertia[2][2]*inertia[0][0]
       -inertia[0][1]*inertia[1][0]-inertia[1][2]*inertia[2][1]-inertia[0][2]*inertia[2][0];
    d = inertia[1][0]*inertia[0][1]*inertia[2][2]+inertia[1][2]*inertia[2][1]*inertia[0][0]+inertia[2][0]*inertia[0][2]*inertia[1][1]
       -inertia[0][0]*inertia[1][1]*inertia[2][2]-inertia[0][1]*inertia[1][2]*inertia[2][0]-inertia[0][2]*inertia[1][0]*inertia[2][1];
    q = (3.0*c-pow(b,2))/9.0;
    r = (9.0*b*c-27*d-2*pow(b,3))/54.0;
    s = cpow( r+csqrt(pow(q,3)+pow(r,2)),1.0/3.0 );
    t = cpow( r-csqrt(pow(q,3)+pow(r,2)),1.0/3.0 );
    
    x1 = -1.0/3.0*b+(s+t);
    x2 = -1.0/3.0*b-0.5*(s+t)+0.5*I*sqrt(3.0)*(s-t);
    x3 = -1.0/3.0*b-0.5*(s+t)-0.5*I*sqrt(3.0)*(s-t);

    eigenvalue[0] = creal(x1);
    eigenvalue[1] = creal(x2);
    eigenvalue[2] = creal(x3);
    
    eigenmin = eigenvalue[0];
    for(j=0;j<3;j++) if(eigenmin > eigenvalue[j]) eigenmin = eigenvalue[j];
    
    bb[0][0] = inertia[0][0] - eigenmin;
    bb[0][1] = inertia[0][1];
    bb[0][2] = inertia[0][2];
    bb[1][0] = inertia[1][0];
    bb[1][1] = inertia[1][1] - eigenmin;
    bb[1][2] = inertia[1][2];
    bb[2][0] = inertia[2][0];
    bb[2][1] = inertia[2][1];
    bb[2][2] = inertia[2][2] - eigenmin;

    if(fabs(vecEE[begin+i][0]) > 1.e-10)
    {
      vecPA[begin+i][0] = vecEE[begin+i][0];
      vecPA[begin+i][1] = -vecPA[begin+i][0]*(bb[0][0]*bb[1][2]-bb[0][2]*bb[1][0])/(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1]);
      vecPA[begin+i][2] = -vecPA[begin+i][0]*(bb[0][1]*bb[1][0]-bb[0][0]*bb[1][1])/(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1]);
    }
    else if(fabs(vecEE[begin+i][1]) > 1.e-10)
    {
      vecPA[begin+i][1] = vecEE[begin+i][1];
      vecPA[begin+i][0] = -vecPA[begin+i][1]*(bb[0][1]*bb[1][2]-bb[0][2]*bb[1][1])/(bb[0][0]*bb[1][2]-bb[0][2]*bb[1][0]);
      vecPA[begin+i][2] = -vecPA[begin+i][1]*(bb[0][0]*bb[1][2]-bb[0][2]*bb[1][0])/(bb[0][0]*bb[1][2]-bb[0][2]*bb[1][0]);
    }
    else if(fabs(vecEE[begin+i][2]) > 1.e-10)
    {
      vecPA[begin+i][2] = vecEE[begin+i][2];
      vecPA[begin+i][0] = -vecPA[begin+i][2]*(bb[0][2]*bb[1][1]-bb[0][1]*bb[1][2])/(bb[0][0]*bb[1][1]-bb[0][1]*bb[1][0]);
      vecPA[begin+i][1] = -vecPA[begin+i][2]*(bb[0][0]*bb[1][2]-bb[0][2]*bb[1][0])/(bb[0][0]*bb[1][1]-bb[0][1]*bb[1][0]);
    }
    vecMAG=sqrt(pow(vecPA[begin+i][0],2)+pow(vecPA[begin+i][1],2)+pow(vecPA[begin+i][2],2));

    vecPA[begin+i][0] /= vecMAG;
    vecPA[begin+i][1] /= vecMAG;
    vecPA[begin+i][2] /= vecMAG;

  }

  return 0;
}


int storePositionSG(t_trxframe fr,float **posSG,float **vecSGa,float **vecSGb,float **vecSGc,int nOM,int nAT,int nSG,int start)
{
  int i,j,k;
  int begin,nOther;
  float temp;
  float **posAT;
  float lvec01[3],lvec02[3],vecMAG;

  posAT = CreateMatrix(nAT,3);

  begin = start*nSG;
  nOther = 0;
  k = 0;
  
  for(i=0;i<nOM;i++)
  {
    for(j=0;j<nAT;j++)
    {
      posAT[j][0] = fr.x[nOther+i*nAT+j][XX];
      posAT[j][1] = fr.x[nOther+i*nAT+j][YY];
      posAT[j][2] = fr.x[nOther+i*nAT+j][ZZ];
    }

    for(j=1;j<nAT-1;j++)
    {
      posSG[begin+k][0] = posAT[j][0];
      posSG[begin+k][1] = posAT[j][1];
      posSG[begin+k][2] = posAT[j][2];
      lvec01[0] = posAT[j][0]-posAT[j-1][0];
      lvec01[1] = posAT[j][1]-posAT[j-1][1];
      lvec01[2] = posAT[j][2]-posAT[j-1][2];
      lvec02[0] = posAT[j+1][0]-posAT[j][0];
      lvec02[1] = posAT[j+1][1]-posAT[j][1];
      lvec02[2] = posAT[j+1][2]-posAT[j][2];
       
      vecSGb[begin+k][0] = lvec01[1]*lvec02[2]-lvec01[2]*lvec02[1];
      vecSGb[begin+k][1] = lvec01[2]*lvec02[0]-lvec01[0]*lvec02[2];
      vecSGb[begin+k][2] = lvec01[0]*lvec02[1]-lvec01[1]*lvec02[0];
      vecMAG = sqrt(pow(vecSGb[begin+k][0],2)+pow(vecSGb[begin+k][1],2)+pow(vecSGb[begin+k][2],2));
      vecSGb[begin+k][0] /= vecMAG; 
      vecSGb[begin+k][1] /= vecMAG;
      vecSGb[begin+k][2] /= vecMAG;

      vecSGc[begin+k][0] = 0.5*(lvec01[0]+lvec02[0]);
      vecSGc[begin+k][1] = 0.5*(lvec01[1]+lvec02[1]);
      vecSGc[begin+k][2] = 0.5*(lvec01[2]+lvec02[2]);
      vecMAG = sqrt(pow(vecSGc[begin+k][0],2)+pow(vecSGc[begin+k][1],2)+pow(vecSGc[begin+k][2],2));
      vecSGc[begin+k][0] /= vecMAG; 
      vecSGc[begin+k][1] /= vecMAG;
      vecSGc[begin+k][2] /= vecMAG;
 
      vecSGa[begin+k][0] = vecSGb[begin+k][1]*vecSGc[begin+k][2]-vecSGb[begin+k][2]*vecSGc[begin+k][1];
      vecSGa[begin+k][1] = vecSGb[begin+k][2]*vecSGc[begin+k][0]-vecSGb[begin+k][0]*vecSGc[begin+k][2];
      vecSGa[begin+k][2] = vecSGb[begin+k][0]*vecSGc[begin+k][1]-vecSGb[begin+k][1]*vecSGc[begin+k][0];
      vecMAG = sqrt(pow(vecSGa[begin+k][0],2)+pow(vecSGa[begin+k][1],2)+pow(vecSGa[begin+k][2],2));
      vecSGa[begin+k][0] /= vecMAG; 
      vecSGa[begin+k][1] /= vecMAG;
      vecSGa[begin+k][2] /= vecMAG;

      k++;
     
    }
  }

  return 0;
}


int storePositionDA(t_trxframe fr,float **posDA,float *phiDA,int nOM,int nAT,int nDA,int start)
{
  int i,j,k;
  int begin,nOther;
  float temp;
  float **posAT;

  float a[3],b[3],vec01[3],vec02[3],vecMAG01,vecMAG02;

  posAT = CreateMatrix(nAT,3);

  begin = start*nDA;
  nOther = 0;
  k = 0;
  for(i=0;i<nOM;i++)
  {
    for(j=0;j<nAT;j++)
    {
      posAT[j][0] = fr.x[nOther+i*nAT+j][XX];
      posAT[j][1] = fr.x[nOther+i*nAT+j][YY];
      posAT[j][2] = fr.x[nOther+i*nAT+j][ZZ];
    }
    for(j=0;j<nAT-3;j++)
    {
      a[0] = posAT[j  ][0]-posAT[j+1][0];
      a[1] = posAT[j  ][1]-posAT[j+1][1];
      a[2] = posAT[j  ][2]-posAT[j+1][2];
      b[0] = posAT[j+2][0]-posAT[j+1][0];
      b[1] = posAT[j+2][1]-posAT[j+1][1];
      b[2] = posAT[j+2][2]-posAT[j+1][2];
      vec01[0] = a[1]*b[2]-a[2]*b[1];
      vec01[1] = a[2]*b[0]-a[0]*b[2];
      vec01[2] = a[0]*b[1]-a[1]*b[0];
      vecMAG01 = sqrt(pow(vec01[0],2)+pow(vec01[1],2)+pow(vec01[2],2));
      vec01[0] /= vecMAG01;
      vec01[1] /= vecMAG01;
      vec01[2] /= vecMAG01;
 
      a[0] = posAT[j+1][0]-posAT[j+2][0];
      a[1] = posAT[j+1][1]-posAT[j+2][1];
      a[2] = posAT[j+1][2]-posAT[j+2][2];
      b[0] = posAT[j+3][0]-posAT[j+2][0];
      b[1] = posAT[j+3][1]-posAT[j+2][1];
      b[2] = posAT[j+3][2]-posAT[j+2][2];
      vec02[0] = a[1]*b[2]-a[2]*b[1];
      vec02[1] = a[2]*b[0]-a[0]*b[2];
      vec02[2] = a[0]*b[1]-a[1]*b[0];
      vecMAG02 = sqrt(pow(vec02[0],2)+pow(vec02[1],2)+pow(vec02[2],2));
      vec02[0] /= vecMAG02;
      vec02[1] /= vecMAG02;
      vec02[2] /= vecMAG02;

      posDA[begin+k][0] = 0.25*(posAT[j][0]+posAT[j+1][0]+posAT[j+2][0]+posAT[j+3][0]); 
      posDA[begin+k][1] = 0.25*(posAT[j][1]+posAT[j+1][1]+posAT[j+2][1]+posAT[j+3][1]); 
      posDA[begin+k][2] = 0.25*(posAT[j][2]+posAT[j+1][2]+posAT[j+2][2]+posAT[j+3][2]); 
      
      phiDA[begin+k] = vec01[0]*vec02[0]+vec01[1]*vec02[1]+vec01[2]*vec02[2];

      k++;
     
    }
  }

  return 0;
}
int autoCorrelationOM(float **OM,float **vecEE,float **vecPA,int binOM[bsz],int nOM,float **corEE01,float **corEE02,float **corPA01,float **corPA02,int origins,int nmax,float z1,float z2, float zsize)
{
  int i,j,k,mz;
  int jstart,kend;
  float temp;

  if(nOM != 0){
    for(i=0; i<nOM; i++)
    {
      if(OM[i][2]>z1 && OM[i][2]<z2) 
      {
        mz = (OM[i][2]-z1)/zsize;
        binOM[mz]++;
      
        for(j=0; j<origins; j++)
        {
          jstart = j*nOM+i;
          for(k=0; k<nmax; k++)
          {
            kend = jstart+k*nOM;
            temp =  vecEE[kend][0]*vecEE[jstart][0]+vecEE[kend][1]*vecEE[jstart][1]+vecEE[kend][2]*vecEE[jstart][2]; 
            corEE01[mz][k] += temp; 
            corEE02[mz][k] += 0.5*(3.0*temp*temp-1.0); 

            temp =  vecPA[kend][0]*vecPA[jstart][0]+vecPA[kend][1]*vecPA[jstart][1]+vecPA[kend][2]*vecPA[jstart][2]; 
            corPA01[mz][k] += temp; 
            corPA02[mz][k] += 0.5*(3.0*temp*temp-1.0); 
          }
        }
      }
    }
  }
  return 0;
}


int autoCorrelationSG(float **SG,float **vecSGa,float **vecSGb,float **vecSGc,int binSG[bsz],int nSG,float **corSG01a,float **corSG02a,float **corSG01b,float **corSG02b,float **corSG01c,float **corSG02c,int origins,int nmax,float z1,float z2, float zsize)
{
  int i,j,k,mz;
  int jstart,kend;
  float temp;

  if(nSG != 0){
    for(i=0; i<nSG; i++)
    {
      if(SG[i][2]>z1 && SG[i][2]<z2) 
      {
        mz = (SG[i][2]-z1)/zsize;
        binSG[mz]++;
      
        for(j=0; j<origins; j++)
        {
          jstart = j*nSG+i;
          for(k=0; k<nmax; k++)
          {
            kend = jstart+k*nSG;
            temp =  vecSGa[kend][0]*vecSGa[jstart][0]+vecSGa[kend][1]*vecSGa[jstart][1]+vecSGa[kend][2]*vecSGa[jstart][2]; 
            corSG01a[mz][k] += temp; 
            corSG02a[mz][k] += 0.5*(3.0*temp*temp-1.0); 

            kend = jstart+k*nSG;
            temp =  vecSGb[kend][0]*vecSGb[jstart][0]+vecSGb[kend][1]*vecSGb[jstart][1]+vecSGb[kend][2]*vecSGb[jstart][2]; 
            corSG01b[mz][k] += temp; 
            corSG02b[mz][k] += 0.5*(3.0*temp*temp-1.0); 

            kend = jstart+k*nSG;
            temp =  vecSGc[kend][0]*vecSGc[jstart][0]+vecSGc[kend][1]*vecSGc[jstart][1]+vecSGc[kend][2]*vecSGc[jstart][2]; 
            corSG01c[mz][k] += temp; 
            corSG02c[mz][k] += 0.5*(3.0*temp*temp-1.0); 
          }
        }
      }
    }
  }
  return 0;
}


int autoCorrelationDA(float **DA,float *phiDA,int binDA[bsz],int nDA,float **corDA01,float **corDA02,float **corDA03,int origins,int nmax,float z1,float z2, float zsize)
{
  int i,j,k,mz;
  int jstart,kend;
  float temp;

  if(nDA != 0){
    for(i=0; i<nDA; i++)
    {
      if(DA[i][2]>z1 && DA[i][2]<z2) 
      {
        mz = (DA[i][2]-z1)/zsize;
        binDA[mz]++;
      
        for(j=0; j<origins; j++)
        {
          jstart = j*nDA+i;
          for(k=0; k<nmax; k++)
          {
            kend = jstart+k*nDA;
            corDA01[mz][k] += phiDA[kend]*phiDA[jstart]; 
            corDA02[mz][k] += phiDA[jstart]; 
            corDA03[mz][k] += pow(phiDA[jstart],2); 
          }
        }
      }
    }
  }
  return 0;
}


float* CreateVector(int rows) { return calloc((unsigned int) rows,sizeof(float)); }


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
