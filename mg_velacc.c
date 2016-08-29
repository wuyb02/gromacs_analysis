#include "mg_velacc.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
  static char *desc[] = {
    "g_velacc with bins\n"
    "mg_velacc -f md01.trr -s md01-331.tpr -n index.ndx -o vac_Hmg.xvg -dt 0.001 -dT 0.001 -snp 10000 -ng 1 -nskipt0 50 -chunk 5.0001 -z1 0 -z2 10 -nobMW\n"
  };

  static bool bMW=TRUE;
  static real binwidth = 0.001;
  static real dT = 0.04;
  static real chunk = 50.0;
  static real z1 = 0.0, z2 = 10.0;
  static int  ng=1, snp=1, nskipt0=1;
  t_pargs pa[] = {
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
    { "-ng",       FALSE, etINT,  {&ng},       "# of index grp to read" },
    { "-dT",       FALSE, etREAL, {&dT},       "time between frames, same as dt, ps" },
    { "-chunk",    FALSE, etREAL, {&chunk},    "chunk size, ps" },
    { "-z1",       FALSE, etREAL, {&z1},       "bin z lower limit" },
    { "-z2",       FALSE, etREAL, {&z2},       "bin z upper limit" },
    { "-nskipt0",  FALSE, etINT,  {&nskipt0},  "# of t0 to skip" },
    { "-bMW",      FALSE, etBOOL, {&bMW},      "Mass Weighted MSD" },
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
    printf("mg_velacc, emoNR!=ng, %d, %d\n", emoNR, ng);
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
  do_velacc(fnTPS, fnTRX, fnRDF,
         ng, isize, grpname, index,
         snp, dT, chunk, bMW, z1, z2, nskipt0);

  thanx(stderr);
  
  return 0;
}
