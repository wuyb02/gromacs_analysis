#include "do_trjcombine.h"

void do_trjcombine(char *fnTRX1, char *fnTRX2, char *fnTRXo,
	           int snp, real dT)
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
  t_trxframe fr1, fr2, frout;
  int status1, status2, status_out;
  int flags = TRX_READ_X;
  rvec *xmem=NULL,*vmem=NULL,*fmem=NULL;
  // gmx_conect gc=NULL;
  int nout;
  static gmx_bool bVels=TRUE, bForce=FALSE;
  output_env_t oenv;

  read_first_frame(oenv,&status1, fnTRX1, &fr1, flags);
  read_first_frame(oenv,&status2, fnTRX2, &fr2, flags);

  nout = fr1.natoms + fr2.natoms;
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
  strcpy(filemode,"a");
  t_trxstatus *trxout = open_trx(fnTRXo,filemode);

  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      frout = fr1;
      frout.natoms = nout;

      frout.x = xmem;
      if(bVels) 
      {
        frout.v = vmem;
      }
      if(bForce) 
      {
        frout.f = fmem;
      }

      for(i=0; i<fr1.natoms; i++) 
      {
        copy_rvec(fr1.x[i],frout.x[i]);
        if (bVels && fr1.bV) 
        {
          copy_rvec(fr1.v[i],frout.v[i]);
        }
        if (bForce && fr1.bF) 
        {
          copy_rvec(fr1.f[i],frout.f[i]);
        }
      }
      for(i=0; i<fr2.natoms; i++) 
      {
        copy_rvec(fr2.x[i],frout.x[i+fr1.natoms]);
        if (bVels && fr2.bV) 
        {
          copy_rvec(fr2.v[i],frout.v[i+fr1.natoms]);
        }
        if (bForce && fr2.bF) 
        {
          copy_rvec(fr2.f[i],frout.f[i+fr1.natoms]);
        }
      }
      frout.bV=TRUE;
      frout.bF=FALSE;

      //write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc); @../../src/gmxlib/trxio.c
      write_trxframe(trxout,&frout,NULL);

      read_next_frame(oenv,status1,&fr1);
      read_next_frame(oenv,status2,&fr2);
    }
  }
}
