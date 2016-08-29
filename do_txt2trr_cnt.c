#include "do_txt2trr_cnt.h"

void do_txt2trr_cnt(char *fnTRX, char *fnTRXx, char *fnTRXv, char *fnTRXo,
	           int snp, real dT, bool bVels)
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

      read_next_txtframe(FnTRXx,xmem,9,nout);
      frout.x = xmem;
      if(bVels) 
      {
        read_next_txtframe(FnTRXv,vmem,9,nout);
        frout.v = vmem;
      }
      if(bForce) 
      {
        frout.f = fmem;
      }

      frout.bV=TRUE;
      frout.bF=FALSE;

      //write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc); @../../src/gmxlib/trxio.c
      write_trxframe(trxout,&frout);
    }
  }
}

int read_next_txtframe(FILE *FnTRXx, rvec *xmem, int nskip, int natom)
{
  char *s = NULL;
  size_t bufsize = 0;

  int i, j, k;
  int id, type;
  float x1, x2, x3;

  for(i=0; i<nskip; i++)
  {
    fgetline(&s, &bufsize, 32768, FnTRXx, 1);
  }

  for(i=0; i<natom; i++)
  {
    fgetline(&s, &bufsize, 32768, FnTRXx, 0);
    sscanf(s, "%d%d%f%f%f", &id, &type, &x1, &x2, &x3);
    xmem[id-1][0]=x1/10;
    xmem[id-1][1]=x2/10;
    xmem[id-1][2]=x3/10;
  }

  return 0;
}
