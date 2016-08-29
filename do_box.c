#include "Mgromacs.h"
#include "Mmatrix.h"

int storePosition(t_trxframe);

void do_box(char *fnTPS, char *fnTRX, int snp)
{
  int        status;
  int        g,natoms,nbin,j0,j1,n,nframes;
  int        **count;
  int        isize_cm=0,nrdf=0,max_i;
  atom_id    *index_cm=NULL;
  unsigned long int *sum;
  t_block    *excl;
  t_trxframe fr;
  excl=NULL;
  int i, j, k;

  t_topology top;
  char title[256];
  rvec *xdum;
  matrix box;
  bool bMW=FALSE, bTop;

/*
  bTop=read_tps_conf(fnTPS,title,&top,&xdum,NULL,box,bMW);
  for(i=0; i<top.atoms.nr; i++)
  {
    printf("%10.5f %10.5f %10d %10d\n",
	top.atoms.atom[i].m,
	top.atoms.atom[i].q,
	top.atoms.atom[i].type,
	top.atoms.atom[i].resnr
	);
  }
  printf("\n%10d\n",
	top.atoms.nres
	);
  printf("\n\n");
*/
  
/*
  int flags = TRX_READ_X;
  read_first_frame(&status, fnTRX, &fr, flags);
  int step;
  for(step=0; step<snp; step++)
  {
    read_next_frame(status,&fr);
    storePosition(fr);
  }
*/
}

int storePosition(t_trxframe fr)
{
  int i, j, k;

  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      printf("%11.3f", fr.box[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  return 0;
}
