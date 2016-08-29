/*
   Original, GROMACS 3.2
   Modified, By Rui Qiao, 05/10/2002, 01/01/2003
   Modified, By Yanbin, 11/22/2007
*/

#include "mg_rdf.h"

int main(int argc,char *argv[])
{
  int i,j,k,ll,mm,nn;
/*
  static char *desc[] = {
    "Original, GROMACS 3.2"
    "Modified, By Yanbin, 11/19/2007"
    "Template for analysis"
  };

  static bool bCM=FALSE;
  static real binwidth=0.001;
  static real Lx, Ly, Lz;
  static int  natom, snp;
  t_pargs pa[] = {
    { "-bin",      FALSE, etREAL, {&binwidth}, "Binwidth (nm)" },
    { "-Lx",       FALSE, etREAL, {&Lx},       "Lx (nm)" },
    { "-Ly",       FALSE, etREAL, {&Ly},       "Ly (nm)" },
    { "-Lz",       FALSE, etREAL, {&Lz},       "Lz (nm)" },
    { "-natom",    FALSE, etINT,  {&natom},    "# of atoms in each molecule" },
    { "-snp",      FALSE, etINT,  {&snp},      "# of frames" },
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
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  fnTRX = ftp2fn_null(efTRX,NFILE,fnm);
  fnRDF = ftp2fn_null(efXVG,NFILE,fnm);
*/

/**
 ** INDEX
  int ng, *isize;
  char **grpname;
  atom_id **index;

  read_index(fnNDX, &ng, isize, index, grpname);
  do_rdf(fnTPS, fnTRX, fnRDF, 
         ng, isize, grpname, index,
         binwidth, Lx, Ly, Lz, natom, snp);

  thanx(stderr);
 **/
 
/**
 ** Matrix operation
 **/
  int MDIM=3, NP=4;
  MMatrix_t A=CreateMatrix(NP, MDIM);
  A->m[0][0]=1; A->m[0][1]=2; A->m[0][2]=2;
  A->m[1][0]=5; A->m[1][1]=4; A->m[1][2]=2;
  A->m[2][0]=6; A->m[2][1]=7; A->m[2][2]=2;
  A->m[3][0]=4; A->m[3][1]=9; A->m[3][2]=2;
  MMatrix_t R=CreateMatrix(NP,1);
  for(i=0;i<R->rows;i++)
  {
    R->m[i][0]=1;
  }
  MMatrix_t AT=MatrixTranspose(A);
  MMatrix_t AA=MatrixProduct(AT,A);
  MMatrix_t AR=MatrixProduct(AT,R);
  PrintMatrix(A);
  PrintMatrix(AA);
  PrintMatrix(R);
//  printf("Det(A)=%11.3f\n", MatrixDet(A));
  MMatrix_t aa=EqnSolver(AA,AR);
  PrintMatrix(aa);
  MatrixNorm(aa);
  PrintMatrix(aa);
  FreeMatrix(A);
  FreeMatrix(AA);
  FreeMatrix(AR);
  FreeMatrix(aa);
  
  return 0;
}
