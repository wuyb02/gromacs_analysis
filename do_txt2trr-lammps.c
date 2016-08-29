#include "do_txt2trr-lammps.h"

void do_txt2trr(char *fnTRX, char *fnTRXx, char *fnTRXv, char *fnTRXo,
	           int snp, int icnt, real dT, bool bVels)
{
  int i, j, k;
  int ll, mm, nn;

  int totalStep = snp;
  int sectionSize  = 1;
  int effSectSize  = sectionSize;       
  int totalSection = totalStep/sectionSize;
  int section;

  printf("Total section: %d\n",totalSection);
  printf("Section size : %d\n",sectionSize);

/**
 ** Trajectory
 **/
  t_trxframe fr, frout;
  int status;
  int flags = TRX_READ_X | TRX_READ_V;
  rvec *xmem=NULL,*vmem=NULL,*fmem=NULL;
  // gmx_conect gc=NULL;
  int nout;
  static bool bForce=FALSE;

  read_first_frame(&status, fnTRX, &fr, flags);

  nout = fr.natoms;
  snew(xmem,nout);
  if(bVels) 
  {
    snew(vmem,nout);
  }
  if(bForce) 
  {
    snew(fmem,nout);
  }

  char filemode[5];
  strcpy(filemode,"w");
  int trxout = open_trx(fnTRXo,filemode);

  FILE *FnTRXx = fopen(fnTRXx, "r");
  FILE *FnTRXv = NULL;
  if(bVels) 
  {
     FnTRXv = fopen(fnTRXv, "r");
  }


  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      frout = fr;
      frout.natoms = nout;
      frout.time = dT*section;

      read_next_txtframe(FnTRXx,xmem,9,nout,icnt);
      frout.x = xmem;
      if(bVels) 
      {
        read_next_txtframe(FnTRXv,vmem,9,nout,icnt);
        frout.v = vmem;
      }
      if(bForce) 
      {
        frout.f = fmem;
      }

      frout.bV=TRUE;
      frout.bF=FALSE;

      //write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc); @../../src/gmxlib/trxio.c
      //int write_trxframe(int fnum,t_trxframe *fr)
      write_trxframe(trxout,&frout);
    }
  }
}

int read_next_txtframe(FILE *FnTRXx, rvec *xmem, int nskip, int natom, int icnt)
{
  char *s = NULL;
  size_t bufsize = 0;

  int i, j, k;
  int id, type, mol;
  float x1, x2, x3;

  for(i=0; i<nskip; i++)
  {
    fgetline(&s, &bufsize, 32768, FnTRXx, 1);
  }

  for(i=0; i<natom; i++)
  {
    fgetline(&s, &bufsize, 32768, FnTRXx, 0);
    sscanf(s, "%d%d%d%f%f%f", &id, &type, &mol, &x1, &x2, &x3);
    if(id-icnt-1<0 || id-icnt-1>=natom)
    {
      printf("id=%d, icnt=%d, natom=%d\n", id, icnt, natom);
      exit(0);
    }
    xmem[id-icnt-1][0]=x1/10;
    xmem[id-icnt-1][1]=x2/10;
    xmem[id-icnt-1][2]=x3/10;
  }

  return 0;
}
