/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_ca.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  static char *desc[] = {
    "Original, GROMACS 3.3\n"
    "Modified, By Rui Qiao, 05/10/2002, 01/01/2003\n"
    "Modified, By Yanbin, 11/14/2010\n"
    "Contact Angle, Water Droplet Center of Mass\n"
  };

  static int  ng=1, snp=1;
  static real dT = 0.05;
  static real zg = 0.44;
  static real Lz = 5.0, R = 5.0;
  static real zbin = 0.01, dA = 0.95;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-ng",       FALSE, etINT,  {&ng},       "# of index grp to read" },
    { "-dT",       FALSE, etREAL, {&dT},       "frame step"},
    { "-zg",       FALSE, etREAL, {&zg},       "graphene sheet z coordinates"},
    { "-Lz",       FALSE, etREAL, {&Lz},       "z upper limit"},
    { "-R",        FALSE, etREAL, {&R},        "R upper limit"},
    { "-zbin",     FALSE, etREAL, {&zbin},     "z bin size"},
    { "-dA",       FALSE, etREAL, {&dA},       "r bin size"},
  };
  #define NPA asize(pa)

  char       *fnTPS,*fnNDX, *fnTRX, *fnRDF, *fnGridR, *fnGridZ;
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,     ffREAD },
    { efTPS, "-s",  NULL,     ffREAD },
    { efNDX, "-n",  "index",    ffREAD },
    { efXVG, "-o",  "wtrho",    ffWRITE },
    { efXVG, "-ox",  "gridR",    ffWRITE },
    { efXVG, "-oy",  "gridZ",    ffWRITE },
  };
  #define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnTRX = ftp2fn_null(efTRX,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  fnRDF = ftp2fn_null(efXVG,NFILE,fnm);
  fnGridR = opt2fn_null("-ox",NFILE,fnm);
  fnGridZ = opt2fn_null("-oy",NFILE,fnm);


  if(emoNR!=ng)
  {
    printf("mg_hbond, emoNR!=ng, %d, %d\n", emoNR, ng);
    exit(0);
  }

/**
 ** INDEX
 **/
  int *isize;
  char **grpname;
  atom_id **index;
  snew(grpname,ng);
  snew(isize,ng);
  snew(index,ng);
  fprintf(stderr,"\nSelect %d group%s\n", ng,ng==1?"":"s");

  read_index(fnNDX, &ng, isize, index, grpname);

  do_ca(fnTPS, fnTRX, fnRDF, fnGridR, fnGridZ, ng, isize, grpname, index, snp, dT, zg, Lz, R, zbin, dA);

  thanx(stderr);
  
  return 0;
}
