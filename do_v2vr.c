#include "do_v2vr.h"

void do_v2vr(char *fnTRX, char *fnTRXo,
	     int snp, real dT, real Xcom, real Ycom)
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
  int status, status_out;
  int flags = TRX_READ_X | TRX_READ_V;
  rvec *xmem=NULL,*vmem=NULL,*fmem=NULL;
  // gmx_conect gc=NULL;
  int nout;
  static bool bVels=TRUE, bForce=FALSE;

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

  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      frout = fr;
      frout.natoms = nout;

      frout.x = xmem;
      frout.v = vmem;

      v2vr(&fr, &frout, Xcom, Ycom);

      frout.bV=TRUE;
      frout.bF=FALSE;

      //write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc); @../../src/gmxlib/trxio.c
      write_trxframe(trxout,&frout);

      read_next_frame(status,&fr);
    }
  }
}

int v2vr(t_trxframe *fr_p, t_trxframe *frout_p, real Xcom, real Ycom)
{
  int i, j, k;
  real Xcom_frame, Ycom_frame;
  rvec ir;
  real ir_norm, vr;

  Xcom_frame=0.0;
  Ycom_frame=0.0;
  for(i=0; i<fr_p->natoms; i++)
  {
    Xcom_frame+=fr_p->x[i][0];
    Ycom_frame+=fr_p->x[i][1];
  }
  Xcom_frame/=fr_p->natoms;
  Ycom_frame/=fr_p->natoms;
  
  if(Xcom_frame-Xcom>0.3 || Xcom_frame-Xcom<=-0.3)
  {
    printf("Warning: Xcom_frame=%f, Xcom=%f\n", Xcom_frame, Xcom);
  }
  if(Ycom_frame-Ycom>0.3 || Ycom_frame-Ycom<=-0.3)
  {
    printf("Warning: Ycom_frame=%f, Ycom=%f\n", Ycom_frame, Ycom);
  }

  for(i=0; i<fr_p->natoms; i++) 
  {
    copy_rvec(fr_p->x[i],frout_p->x[i]);

    ir[0]=fr_p->x[i][0]-Xcom_frame;
    ir[1]=fr_p->x[i][1]-Ycom_frame;
    ir[2]=0.0;
    ir_norm=sqrt(ir[0]*ir[0]+ir[1]*ir[1]);
    ir[0]=ir[0]/ir_norm;
    ir[1]=ir[1]/ir_norm;

    vr = ir[0]*fr_p->v[i][0]+ir[1]*fr_p->v[i][1];
    frout_p->v[i][0]=vr;
    frout_p->v[i][1]=0.0;
    frout_p->v[i][2]=0.0;
  }

  return 0;
}
