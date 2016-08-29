/*
   Original, GROMACS 3.3
   Modified, By Yanbin, 11/14/2010
*/

#include "mg_txt2trr.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  static char *desc[] = {
    "lmp trj file to gmx trr file\n"
    "\n"
    "mv dump.all.trrx dump.all.trrx.xvg\n"
    "mv dump.all.trrv dump.all.trrv.xvg\n"
    "./mg_txt2trr -f md011-tmplat.trr -fx dump.all.trrx.xvg -fv dump.all.trrv.xvg -fo md011o.trr -snp 200000 -icnt 0 -dT 0.001\n"
    "\n"
    "lmp input (real unit)\n"
    "variable vx1000 atom vx*1000\n"
    "variable vy1000 atom vy*1000\n"
    "variable vz1000 atom vz*1000\n"
    "dump    1 all custom 500 dump.lammpstrj id mol type x y z\n"
    "dump    2 fluidmols custom ${dumpvel_itv} dump.all.trrx id mol type x y z\n"
    "dump    3 fluidmols custom ${dumpvel_itv} dump.all.trrv id mol type v_vx1000 v_vy1000 v_vz1000\n"
  };

  static int  snp=1;
  static int icnt=0;
  static real dT = 0.05;
  static bool bVels = TRUE;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-icnt",     FALSE, etINT,  {&icnt},     "environment atoms before the molecule" },
    { "-dT",       FALSE, etREAL, {&dT},       "frame step"},
    { "-bVels",    FALSE, etBOOL, {&bVels},    "velocity reading"},
  };
  #define NPA asize(pa)

  char *fnTRX, *fnTRXx, *fnTRXv, *fnTRXo;
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD },
    { efXVG, "-fx",  NULL,     ffREAD },
    { efXVG, "-fv",  NULL,     ffOPTRD },
    { efTRX, "-fo",  "trjout",     ffWRITE }
  };
  #define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTRX  = opt2fn_null("-f",NFILE,fnm);
  fnTRXx = opt2fn_null("-fx",NFILE,fnm);
  fnTRXv = opt2fn_null("-fv",NFILE,fnm);
  fnTRXo = opt2fn_null("-fo",NFILE,fnm);

  do_txt2trr(fnTRX, fnTRXx, fnTRXv, fnTRXo, snp, icnt, dT, bVels);

  thanx(stderr);
  
  return 0;
}
