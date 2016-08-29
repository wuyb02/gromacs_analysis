/*
   Original, GROMACS 3.2
   Modified, By Yanbin, 11/19/2007
*/

#include "mg_bbcn.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  static char *desc[] = {
    "Original, GROMACS 3.2\n"
    "Modified, By Rui Qiao, 05/10/2002, 01/01/2003\n"
    "Modified, By Sony\n"
    "Modified, By Yanbin, 02/25/2008\n"
    "Template for analysis\n"
  };

  static int  ng=1, snp=1;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-ng",       FALSE, etINT,  {&ng},       "# of index grp to read" },
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
    printf("mg_bbcn, emoNR!=ng, %d, %d\n", emoNR, ng);
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
  do_bbcn(fnTPS, fnTRX, fnRDF,
	  ng, isize, grpname, index,
	  snp);

  thanx(stderr);
  
  return 0;
}