/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_addvr.h"

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
  static real Xcom = 2.5;
  static real Ycom = 2.5;
  static real vr = 0.2;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-dT",       FALSE, etREAL, {&dT},       "frame step"},
    { "-Xcom",     FALSE, etREAL, {&Xcom},     "Xcom"},
    { "-Ycom",     FALSE, etREAL, {&Ycom},     "Ycom"},
    { "-vr",       FALSE, etREAL, {&vr},       "vr"},
  };
  #define NPA asize(pa)

  char *fnTRX, *fnTRXo;
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,     ffREAD },
    { efTRX, "-fo",  "trjout",     ffWRITE }
  };
  #define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTRX = opt2fn_null("-f",NFILE,fnm);
  fnTRXo = opt2fn_null("-fo",NFILE,fnm);

  do_addvr(fnTRX, fnTRXo, snp, dT, Xcom, Ycom, vr);

  thanx(stderr);
  
  return 0;
}
