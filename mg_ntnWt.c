/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_ntnWt.h"

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
  static real z1 = 0.0, z2 = 10.0, xnt = 1.7121, ynt = 1.7201, rnt = 0.35;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-ng",       FALSE, etINT,  {&ng},       "# of index grp to read" },
    { "-dT",       FALSE, etREAL, {&dT},       "frame step"},
    { "-z1",       FALSE, etREAL, {&z1},       "z1"},
    { "-z2",       FALSE, etREAL, {&z2},       "z2"},
    { "-xnt",      FALSE, etREAL, {&xnt},      "xnt"},
    { "-ynt",      FALSE, etREAL, {&ynt},      "ynt"},
    { "-rnt",      FALSE, etREAL, {&rnt},      "rnt"},
  };
  #define NPA asize(pa)

  char       *fnTPS,*fnNDX, *fnTRX, *fnRDF;
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,     ffREAD },
    { efTPS, "-s",  NULL,     ffREAD },
    { efNDX, "-n",  "ndx",    ffREAD },
    { efXVG, "-o",  "rdf",    ffWRITE },
  };
  #define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnTRX = ftp2fn_null(efTRX,NFILE,fnm);
  fnRDF = ftp2fn_null(efXVG,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);

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

/**
 ** Orientation autocorrelation
 **/
  do_ntnWt(fnTPS, fnTRX, fnRDF, ng, isize, grpname, index, snp, dT, z1, z2, xnt, ynt, rnt);

  thanx(stderr);
  
  return 0;
}
