/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_trjcombine.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  const char *desc[] = {
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

  char *fnTRX1, *fnTRX2, *fnTRXo;
  t_filenm   fnm[] = {
    { efTRX, "-f1",  NULL,     ffREAD },
    { efTRX, "-f2",  NULL,     ffREAD },
    { efTRX, "-fo",  NULL,     ffWRITE }
  };
  #define NFILE asize(fnm)
  output_env_t oenv;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);

  fnTRX1 = opt2fn_null("-f1",NFILE,fnm);
  fnTRX2 = opt2fn_null("-f2",NFILE,fnm);
  fnTRXo = opt2fn_null("-fo",NFILE,fnm);

  do_trjcombine(fnTRX1, fnTRX2, fnTRXo, snp, dT);

  thanx(stderr);
  
  return 0;
}
