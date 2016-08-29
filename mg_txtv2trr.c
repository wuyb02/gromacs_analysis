/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_txtv2trr.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  static char *desc[] = {
    "Original, GROMACS 3.3\n"
    "Modified, By Yanbin, 09/07/2012\n"
    "TRR Combine\n"
  };

  static int  snp=1;
  static real dT = 0.05;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-dT",       FALSE, etREAL, {&dT},       "frame step"},
  };
  #define NPA asize(pa)

  char *fnTRX, *fnTRXx, *fnTRXv, *fnTRXo;
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD },
    { efXVG, "-fv",  NULL,     ffOPTRD },
    { efTRX, "-fo",  "trjout",     ffWRITE }
  };
  #define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTRX  = opt2fn_null("-f",NFILE,fnm);
  fnTRXv = opt2fn_null("-fv",NFILE,fnm);
  fnTRXo = opt2fn_null("-fo",NFILE,fnm);

  do_txtv2trr(fnTRX, fnTRXv, fnTRXo, snp, dT);

  thanx(stderr);
  
  return 0;
}
