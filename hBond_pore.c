// this code calculates 
// 1. neighbour[]
//    bins possible HBond Pair not considering angle constraints
//    by donor position along radial axis
//    per water molecule
// 2. densWt[]
//    bins water
//    by oxygen position along radial axis
//    per bins volume
// 3. hBond[]
//    bins donor 
//    by donor position along radial axis
//    per water molecule
// 4. hBond2[]
//    bins acceptor
//    by acceptor position radial axis
//    per water molecule
// 5. OH_Oangle[][]
//    bins OH--O angle or Hbond
//    by donor position along radial aixs
//    by OH--O angle
//    no normalization
// 6. O_Oangle[][]
//    bins O--O orientation
//    by donor position along radial axis
//    by O--O angle w. r. to the channel wall
//    no normalization
//
// Original, Rui Qiao, Jan, 15, 2004
// Modified, Yanbin WU, June, 3, 2007
// 
/*
./template -s topol.tpr -f traj.xtc -nice 0 -b 10 -e 690 -dt 0.01 -nWa 0 -nWt 359 -nNa 0 -nCl 0 -nbin 10 -Lx 3.2 -Ly 3.2 -Lz 7.982 -R 1.0 -Rx 1.6 -Ry 1.6 -rcut 0.35 -acut 30 -rlist 1.0 -nsp 68001 -tUpdate 0.01

*/
//
 
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#define  maxList 200
#define  PI 3.1415926

int nAng = 60;
double normNEW(float a[3]);
float pbc(float dx,float Lx);
float* CreateVector(int rows);
double* CreateVector2(int rows);
double i_prod(float a[3],float b[3]);
float cos_angle_new(float a[3],float b[3]);
float **CreateMatrix(int rows,int cols);
int **CreateMatrixInt(int rows,int cols);
double ** CreateMatrix2(int rows,int cols);
float isHbond(float d[3],float h[3],float a[3],float rcut,float acut,float L[3]);
int storePosition(t_trxframe fr,float **positionO,float **positionH1,float **positionH2,
		  float **positionNa,float **positionCl,int nWt,int nNa,int nCl,int start);
float dist2(float **pos,int i,int j,float L[3]);
int neighbourList(int nWt,int **List,float rcut2,float **posO,float L[3]);
int init(double *densWt,double *densNa,double *densCl,double *hBond,double *hBond2,
	 double **OH_Oangle,double **O_Oangle,int nbin);
int analyze_neighbour(float rcut,float **posO,int nWt,int **list,double *neighbour,float L[3],
		      int nbin,double *densWt,float R, float ccx, float ccy);
int analyze_Hbonding(float rcut,float acut,float **posO,float **posH1,float **posH2,
		     int nWt,int nbin,float L[3],double *densWt,double *hBond,
		     double **OH_Oangle,double **O_Oangle,double *hBond2,int **list,float R, float ccx, float ccy);
int scaling(double *densWt,double *neighbour,double *hBond,double *h2Bond,double **OH_Oangle,
	    double **O_Oangle,int nbin,float L[3],int nsp,float * volBin);
int report(FILE *hbF,FILE *hbAngDisF,FILE *hbAngWalF,double *densW,double *neighbour,
	   double *hBond,double *hBond2,double **OH_Oangle,double **O_Oangle,
	   int nbin,float acut,float *rbin);


int main(int argc,char *argv[])
{
  static char *desc[] = {"Compute the hydrogen bond density across the cahannel"};
  
  static int nWt = 0;
  static int nWa = 0;
  static int nNa = 0;
  static int nCl = 0;
  static int nbin  = 121;
  static int nsp = 0;
  
  static float Lx = 4.655;
  static float Ly = 4.428;
  static float Lz = 3.487;
  static float R = 1.0;
  static float ccx= 2.15;
  static float ccy= 2.25;
  
  static float rcut = 0.35; // cutoff distance between hydrogen and acceptor
  static float acut = 30;   // cutoff angle between donor-hydro and hydro-accept vectors
  static float rlist = 0.8; // cutoff for the neighbour list
  static int   freqUpdate = 200; // frequency of update
  static float   tUpdate = 0.400;    // time interval between updates of neighbour list(ps) 
  static float   dt = 0.01;          // time interval between saved frames (ps)

  t_pargs pa[] = {
    { "-nWa", FALSE, etINT, {&nWa},
      "wall atom number"
    },
    { "-nWt", FALSE, etINT, {&nWt},
      "water molecule number"
    },
    { "-nNa", FALSE, etINT, {&nNa},
      "Na+ ion number"
    },
    { "-nCl", FALSE, etINT, {&nCl},
      "Cl- ion number"
    },
    { "-nbin", FALSE, etINT, {&nbin},
      "number of bins across the channel ((z-direction}"
    },
    { "-R", FALSE, etREAL, {&R},
      "pore radius (nm)"
    },
    { "-Rx", FALSE, etREAL, {&ccx},
      "pore center x (nm)"
    },
    { "-Ry", FALSE, etREAL, {&ccy},
      "pore center y (nm)"
    },
    { "-Lx", FALSE, etREAL, {&Lx},
      "Box in x-direction (nm)"
    },
    { "-Ly", FALSE, etREAL, {&Ly},
      "Box in y-direction (nm)"
    },
    { "-Lz", FALSE, etREAL, {&Lz},
      "channel width (nm)"
    },
    { "-rcut", FALSE, etREAL, {&rcut},
      "cutoff distance between hydro-acceptor (nm)"
    },
    { "-acut", FALSE, etREAL, {&acut},
      "cutoff angle between donor/hydro - hydro-acceptor vector (degree)"
    },
    { "-nsp", FALSE, etINT, {&nsp},
      "number of frames available"
    },
    { "-rlist", FALSE, etREAL, {&rlist},
      "cutoff distance between hydro-acceptor (nm)"
    },
    { "-tUpdate", FALSE, etREAL, {&tUpdate},
      "time interval between updates of neighbour list(ps)"
    },
    { "-dt", FALSE, etREAL, {&dt},
      "interval between saved frames (ps)"
    },
  };
  
  
  // variables specifically for this code
  float  L[3];
  float  *rbin,*volBin;
  float  **posO,**posH1,**posH2,**posNa,**posCl;
  double  *densWt,*densNa,*densCl;
  double  *hBond,**OH_Oangle,**O_Oangle;    // home Oxygen as donor
  double  *hBond2; // home Oxygen as acceptor
  double  *neighbour;
  float  rlist2;
  int    **list;
  FILE   *hbFile = fopen("hBondDensity.dat","w");
  FILE   *hbAngDisFile = fopen("hBond-AngleDist.dat","w");
  FILE   *hbAngWalFile = fopen("O-O-AngleWall.dat","w");

  // others
  int i,j;

  // Gromacs stuff
  t_topology top;   char       title[STRLEN];   t_trxframe fr;
  rvec       *xtop; matrix     box;   int        status;
  int        flags = TRX_READ_X; 
  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  }; 
  #define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  posO   = CreateMatrix(nWt+1,3);
  posH1  = CreateMatrix(nWt+1,3);
  posH2  = CreateMatrix(nWt+1,3);
  posNa  = CreateMatrix(nNa+1,3);
  posCl  = CreateMatrix(nCl+1,3);
  rbin   = CreateVector(nbin+1000);
  volBin  = CreateVector(nbin+1000);

  densWt     = CreateVector2(nbin+1000);
  densNa     = CreateVector2(nbin+1000);
  densCl     = CreateVector2(nbin+1000);
  hBond      = CreateVector2(nbin+1000);
  neighbour  = CreateVector2(nbin+1000);
  OH_Oangle  = CreateMatrix2(nbin+1000,nAng); // OH_Oangle: angle of the H-bond, divide into 60 sections
  O_Oangle   = CreateMatrix2(nbin+1000,nAng); // H--O angle with respect to the nearest channel wall 
  hBond2     = CreateVector2(nbin+1000);
  list       = CreateMatrixInt(nWt,maxList); // neighbourlist
  rlist2     = rlist*rlist;
  freqUpdate = tUpdate/dt;

  for(i=0;i<nbin+1;i++)
    rbin[i] = i*R/nbin;
  for(i=0;i<nbin;i++)
    volBin[i] = PI*(pow(rbin[i+1],2)-pow(rbin[i],2))*Lz;

  for(i=0;i<nbin;i++)
    rbin[i] = (i+0.5)*R/nbin;//we can use the same variable to position rbin in the center of the bin

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
   
  L[0] = Lx; L[1] = Ly; L[2] = Lz;
  init(densWt,densNa,densCl,hBond,hBond2,OH_Oangle,O_Oangle,nbin);
  for(i=0;i<nbin;i++)
    neighbour[i] = 0;

  printf("Start to processing data\n");
  for(i=0;i<nsp-1;i++)
  {
    read_next_frame(status,&fr);
    storePosition(fr,posO,posH1,posH2,posNa,posCl,nWt,nNa,nCl,nWa);

    if(i%freqUpdate==0)
      neighbourList(nWt,list,rlist2,posO,L);
      
    analyze_Hbonding(rcut,acut,posO,posH1,posH2,nWt,nbin,L,densWt,hBond,
                     OH_Oangle,O_Oangle,hBond2,list,R,ccx,ccy);
    //bins possible HBond Pair not considering angle constraints
    //by donor position along radial axis
    analyze_neighbour(rcut,posO,nWt,list,neighbour,L,nbin,densWt,R,ccx,ccy);
  }
  
  scaling(densWt,neighbour,hBond,hBond2,OH_Oangle,O_Oangle,nbin,L,nsp,volBin);
  report(hbFile,hbAngDisFile,hbAngWalFile,densWt,neighbour,hBond,hBond2,OH_Oangle,
         O_Oangle,nbin,acut,rbin);

  fclose(hbFile);
  fclose(hbAngDisFile);
  fclose(hbAngWalFile);

  thanx(stderr);
  
  return 0;
}

int analyze_neighbour(float rcut,float **posO,int nWt,int **list,double *neighbour,float L[3],
                      int nbin,double *densWt,float R, float ccx, float ccy)
{
  int i,j,m,me;
  float rSpacing, r_dist_O;
  float rcut2;

  rSpacing = R/nbin;
  rcut2 = rcut*rcut;
  
  for(i=0;i<nWt;i++)
  {
    /*  distance goes here */
    r_dist_O = sqrt((posO[i][0]-ccx)*(posO[i][0]-ccx)+(posO[i][1]-ccy)*(posO[i][1]-ccy));
    /* which bin the Oxygen is in */

    me = r_dist_O/rSpacing;
    //      densWt[me] +=1.0;
    for(m=0; m<list[i][0]; m++)
    {
      j = list[i][m+1];
      if(dist2(posO,i,j,L)<rcut2)
        neighbour[me]+= 1.0;
    }
  }

  return 0;
}
float dist2(float **pos,int i,int j,float L[3])
{
  float r2 = 0,dr;
  int k;
  for(k=0;k<3;k++)
    {
      dr = pbc(pos[i][k]-pos[j][k],L[k]);
      r2+= dr*dr;
    }
  return r2;
}

int neighbourList(int nWt,int **List,float rcut2,float **posO,float L[3])
{
  int i,j;
  for(i=0;i<nWt;i++)
  {
    List[i][0] = 0;
    for(j=0;j<nWt;j++)
    {
      if(i!=j && dist2(posO,i,j,L)<rcut2)
      {
        List[i][0]++;
        List[i][List[i][0]] = j;
      }
      if(List[i][0]>maxList)
      {
        printf("Number of atoms in neighbour exceeds %d\n",maxList);
        exit(0);
      }
    }
  }
  return 0;
}

int fillAcceptorInformation(float acceptor[3],float **posO,int i)
{
  acceptor[0] = posO[i][0];
  acceptor[1] = posO[i][1];
  acceptor[2] = posO[i][2];
 
  return 0;
}
                                                                                                         
int fillDonorInformation(float donor[3],float hydrogen1[3],float hydrogen2[3],float **posO,
                         float **posH1,float **posH2,int i)
{
  donor[0] = posO[i][0];
  donor[1] = posO[i][1];
  donor[2] = posO[i][2];
 
  hydrogen1[0] = posH1[i][0];
  hydrogen1[1] = posH1[i][1];
  hydrogen1[2] = posH1[i][2];
 
  hydrogen2[0] = posH2[i][0];
  hydrogen2[1] = posH2[i][1];
  hydrogen2[2] = posH2[i][2];
 
  return 0;
}

int analyze_Hbonding(float rcut,float acut,float **posO,float **posH1,float **posH2,
                     int nWt,int nbin,float L[3],double *densWt,double *hBond,
                     double **OH_Oangle,double **O_Oangle,double *hBond2,int **list,float R, float ccx, float ccy)
{
  int i,j,me,meAng,meAccept,m;
  float dAngHB = acut/nAng;
  float dAngWa = 180.0/nAng;
  float donor[3],hydrogen1[3],hydrogen2[3],acceptor[3],vecOH[3];
  float radialvec[3]={0,0,0};
  float yes,angle;
  float rSpacing,r_dist_O;
                                                                                                             
  rSpacing = R/nbin;
                                                                                                             
  for(i=0;i<nWt;i++)
  {
    fillDonorInformation(donor,hydrogen1,hydrogen2,posO,posH1,posH2,i);
 
    r_dist_O = sqrt((posO[i][0]-ccx)*(posO[i][0]-ccx)+(posO[i][1]-ccy)*(posO[i][1]-ccy));
    radialvec[0]=posO[i][0]-ccx;
    radialvec[1]=posO[i][1]-ccy;
    radialvec[2]=0;

    // bins water
    // by position along radial axis
    me = r_dist_O/rSpacing;
    densWt[me] += 1.0;
    for(m=0; m<list[i][0]; m++)
    {
      j = list[i][m+1];
      fillAcceptorInformation(acceptor, posO, j);

      // first O-H1 donor pair
      yes = isHbond(donor,hydrogen1,acceptor,rcut,acut,L);
      if(yes!=2003)
      {	      
        // bins donor 
        // by donor position along radial axis
        hBond[me]+=1.0;

        // bins OH--O angle or Hbond
        // by donor position along radial aixs
        // by OH--O angle
        meAng = yes/dAngHB;
        OH_Oangle[me][meAng]+=1.0;
      
        // bins O--O orientation
        // by donor position along radial axis
        // by O--O angle w. r. to the channel wall
        vecOH[0] = pbc(posO[j][0]-posO[i][0],L[0]);
        vecOH[1] = pbc(posO[j][1]-posO[i][1],L[1]);
        vecOH[2] = pbc(posO[j][2]-posO[i][2],L[2]);	      

       	angle = acos(cos_angle_new(vecOH,radialvec))/3.1415926*180;
      
        meAng = angle/(180.0/nAng);	      
        O_Oangle[me][meAng]+=1.0;

        // bins acceptor
        // by acceptor position radial axis
        r_dist_O = sqrt((posO[j][0]-ccx)*(posO[j][0]-ccx)+(posO[j][1]-ccy)*(posO[j][1]-ccy));
        meAccept = r_dist_O/rSpacing;
        hBond2[meAccept]+=1.0;		
      }
  
      // second O-H2 donor pair
      yes = isHbond(donor,hydrogen2,acceptor,rcut,acut,L);
      if(yes!=2003)
      {
        hBond[me]+=1.0;
        meAng = yes/dAngHB;
        OH_Oangle[me][meAng]+=1.0;

        // compute the O--H angle with respect to the channel wall
        vecOH[0] = pbc(posO[j][0]-posO[i][0],L[0]);;
        vecOH[1] = pbc(posO[j][1]-posO[i][1],L[1]);;
        vecOH[2] = pbc(posO[j][2]-posO[i][2],L[2]);;
      
       	angle = acos(cos_angle_new(vecOH,radialvec))/3.1415926*180;

        meAng = angle/(180.0/nAng);
        O_Oangle[me][meAng]+=1.0;

        r_dist_O = sqrt((posO[j][0]-ccx)*(posO[j][0]-ccx)+(posO[j][1]-ccy)*(posO[j][1]-ccy));
        meAccept = r_dist_O/rSpacing;
                                                                                                            
        hBond2[meAccept]+=1.0;
      }
    }
  }

  return 0;
}

// the output of OH_Oangle has the number of hits in each angle bin, when visualizing
// the results, add all the contributions together to find the total hits!
// O_Oangle follows the same rule
int scaling(double *densWt,double *neighbour,double *hBond,double *h2Bond,double **OH_Oangle,
            double **O_Oangle,int nbin,float L[3],int nsp,float* volBin )
{
  int i,j;

  for(i=0;i<nbin;i++)
  {
    if(densWt[i]!=0)
    {
      hBond[i] /= densWt[i];
      h2Bond[i] /= densWt[i];
      neighbour[i] /= densWt[i];
    }
    densWt[i] /= volBin[i]*(nsp-1);
  }

  for(i=0;i<nbin;i++)
  {
    for(j=0;j<nAng;j++)
    {
      OH_Oangle[i][j] /= nsp-1;
      O_Oangle[i][j]  /= nsp-1;
    }
  }

  return 0;
}

int report(FILE *hbF,FILE *hbAngDisF,FILE *hbAngWalF,double *densW,double *neighbour,
	   double *hBond,double *hBond2,double **OH_Oangle,double **O_Oangle,
	   int nbin,float acut,float *rbin)
{
  int i,j;
  
  for(i=0;i<nbin;i++)
    fprintf(hbF,"%6.3f %6.3f %6.3f %6.3f %6.3f\n",rbin[i],densW[i],hBond[i],hBond2[i],neighbour[i]);

  for(i=0;i<nbin;i++)
  {
    for(j=0;j<nAng;j++)
      fprintf(hbAngDisF,"%9.6e ",OH_Oangle[i][j]);
    fprintf(hbAngDisF,"\n");
  }

  for(i=0;i<nbin;i++)
  {
    for(j=0;j<nAng;j++)
      fprintf(hbAngWalF,"%9.6e ",O_Oangle[i][j]);
    fprintf(hbAngWalF,"\n");
  }

  return 0;
}

double normNEW(float a[3])
{
  double aa=a[0];
  double aaa=a[1];
  double aaaa=a[2];

  double res = aa*aa+aaa*aaa+aaaa*aaaa;

  return sqrt(res);
}
double i_prod(float a[3],float b[3])
{
  double aa=a[0];
  double aaa=a[1];
  double aaaa=a[2];

  double bb=b[0];
  double bbb=b[1];
  double bbbb=b[2];

  double res = aa*bb+aaa*bbb+aaaa*bbbb;

  return res;
}
// cosine of the angle between vector a and b
float cos_angle_new(float a[3],float b[3])
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  double r_ab_2,r_a,r_b;
  double cosAng;
  
  r_ab_2 = i_prod(a,b);
  r_a  = normNEW(a);
  r_b  = normNEW(b);
  
  cosAng = r_ab_2/(r_a*r_b);
  
  if(cosAng<-1 || cosAng > 1.0)
    {
      printf("Houston, we got a problem: |cos(alpha)| > 1 !\n");
      exit(0);
    }  
  else
    return cosAng;
}

float pbc(float dx,float Lx)
{
  if(dx>Lx/2)
    dx -= Lx;
  else if(dx<-Lx/2)
    dx += Lx;
  
  return dx;
}

// if Hbond is found, return the angle of Hbond, otherwise return 2003
float isHbond_old(float d[3],float h[3],float a[3],
		  float rcut,float acut,float L[3]) // use Gromacs implementation
{    
  int i,j;
  float r_ha[3],r_dh[3];
  float r_ha2=0,r_dh2=0;
  double CosAngle;

  acut = cos(acut/180.0*3.1415926); // convert from angle to cos(angle)
  
  // compute vector: hydrogen->acceptor and hydrogen->donor
  for(i=0;i<3;i++)
    {
      r_ha[i] = pbc(a[i]-h[i],L[i]); // r_ha = a - h
      r_dh[i] = pbc(h[i]-d[i],L[i]); // r_dh = h - d
      
      if(r_ha[i]>rcut || r_ha[i]< -rcut)
	return 2003;
    }
  
  CosAngle = cos_angle_new(r_ha,r_dh);
  
  if(CosAngle < acut) // small angle gives big cos(angle)
    return 2003;
  else
    return acos(CosAngle)/3.1415926*180;  
}

float isHbond(float d[3],float h[3],float a[3],
              float rcut,float acut,float L[3])
{    
  int i,j;
  float r_da[3],r_dh[3];
  float r_da2=0,r_dh2=0;
  double CosAngle;

  acut = cos(acut/180.0*3.1415926); // convert from angle to cos(angle)
  
  // compute vector: donor->acceptor and hydrogen->donor
  for(i=0;i<3;i++)
  {
    r_da[i] = pbc(a[i]-d[i],L[i]); // r_da = a - d
    r_dh[i] = pbc(h[i]-d[i],L[i]); // r_dh = h - d
    
    if(r_da[i]>rcut || r_da[i]< -rcut)
      return 2003;
  }
  
  if(i_prod(r_da,r_da)>rcut*rcut)
    return 2003;
  else
  {
    CosAngle = cos_angle_new(r_da,r_dh);
    
    if(CosAngle < acut) // small angle gives big cos(angle)
      return 2003;
    else
      return acos(CosAngle)/3.1415926*180;  
  }
}

int **CreateMatrixInt(int rows,int cols)
{
  int i;
  int **m;
  m = calloc((unsigned int) rows,sizeof(int *));
  for (i=0; i < rows; i++) {
    m[i] = calloc((unsigned int) cols,sizeof(int ));
  }
  
  return m;
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
double ** CreateMatrix2(int rows,int cols)
{
   int	i;
   double  **m;

   m = calloc((unsigned int) rows,sizeof(double *));
   for (i=0; i < rows; i++) {
      m[i] = calloc((unsigned int) cols,sizeof(double ));
   }

   return m;
}
double* CreateVector2(int rows) { return calloc((unsigned int) rows,sizeof(double)); }
float* CreateVector(int rows) { return calloc((unsigned int) rows,sizeof(float)); }
int storePosition(t_trxframe fr,float **positionO,float **positionH1,float **positionH2,
		  float **positionNa,float **positionCl,int nWt,int nNa,int nCl,int start
		  )
{
  int i;
  int begin,nOther;
  float zO,zH1,zH2,xO,xH1,xH2,yO,yH1,yH2,xxx,yyy,zzz;


  nOther= start;
  for(i=0;i<nWt;i++)
    {
      zO = fr.x[nOther+i*3][ZZ] ;
      zH1= fr.x[nOther+i*3+1][ZZ];
      zH2= fr.x[nOther+i*3+2][ZZ];

      xO = fr.x[nOther+i*3][XX] ;
      xH1= fr.x[nOther+i*3+1][XX];
      xH2= fr.x[nOther+i*3+2][XX];      

      yO = fr.x[nOther+i*3][YY] ;
      yH1= fr.x[nOther+i*3+1][YY];
      yH2= fr.x[nOther+i*3+2][YY];

      positionO[i][0] = xO;
      positionO[i][1] = yO;
      positionO[i][2] = zO;

      positionH1[i][0] = xH1;
      positionH1[i][1] = yH1;
      positionH1[i][2] = zH1;

      positionH2[i][0] = xH2;
      positionH2[i][1] = yH2;
      positionH2[i][2] = zH2;

    }  
nOther=start+nWt*3;
  for(i=0;i<nNa;i++)
    {
      //      printf("i = %d of %d %f %f %f\n",i,nNa,fr.x[start+i][XX],fr.x[start+i][YY],fr.x[start+i][ZZ]-lb);
      positionNa[i][0] = fr.x[nOther+i][XX];
      positionNa[i][1] = fr.x[nOther+i][YY];
      positionNa[i][2] = fr.x[nOther+i][ZZ];
    }
  //  printf("Na+ ok ...\n");
  for(i=0;i<nCl;i++)
    {
      positionCl[i][0] = fr.x[nOther+nNa+i][XX];
      positionCl[i][1] = fr.x[nOther+nNa+i][YY];
      positionCl[i][2] = fr.x[nOther+nNa+i][ZZ];
    }

  return 0;
}

int init(double *densWt,double *densNa,double *densCl,double *hBond,double *hBond2,
	 double **OH_Oangle,double **O_Oangle,int nbin)
{
  int i,j;
  for(i=0;i<nbin;i++)
    {
      densNa[i] = densCl[i] = densWt[i] = hBond[i] = hBond2[i] = 0;
      for(j=0;j<nAng;j++)
	OH_Oangle[i][j] = O_Oangle[i][j] = 0;
    }

  return 0;
}
