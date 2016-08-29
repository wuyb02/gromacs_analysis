#include "do_txtv2trr.h"

void do_txtv2trr(char *fnTRX, char *fnTRXv, char *fnTRXo,
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
  snew(vmem,nout);

  char filemode[5];
  strcpy(filemode,"w");
  int trxout = open_trx(fnTRXo,filemode);

  FILE *FnTRXv = fopen(fnTRXv, "r");
  fread_patternskip(FnTRXv, "Pzz totvr");

  for(section=0;section<totalSection;section++)
  {
    for(i=0;i<sectionSize;i++)
    {
      frout = fr;
      frout.natoms = nout;
      frout.time = dT*section;

      xmem[0][0] = 2.5;
      xmem[0][1] = 2.5;
      xmem[0][2] = 1.0;
      frout.x = xmem;

      read_next_txtvframe(FnTRXv,vmem);
      frout.v = vmem;

      frout.bV=TRUE;
      frout.bF=FALSE;

      //write_trxframe(t_trxstatus *status,t_trxframe *fr,gmx_conect gc); @../../src/gmxlib/trxio.c
      write_trxframe(trxout,&frout);
    }
  }
}

int read_next_txtvframe(FILE *FnTRXv, rvec *vmem)
{
  char *s = NULL;
  size_t bufsize = 0;

  int step; 
  float temp, lx, ly, lz, pxx, pyy, pzz, totvr;

  fgetline(&s, &bufsize, 32768, FnTRXv, 0);
  sscanf(s, "%d%f%f%f%f%f%f%f%f", &step, &temp, &lx, &ly, &lz, &pxx, &pyy, &pzz, &totvr);
  //printf("step=%10d, totvr=%10.5f\n", step, totvr);
  //vmem[0][0]=totvr/10;
  vmem[0][0]=totvr*10;
  vmem[0][1]=0.0;
  vmem[0][2]=0.0;

  return 0;
}

int fread_patternskip(FILE *FnTRXv, char *pattern)
{
  char *s = NULL;
  size_t bufsize = 0;
  int flag, flag_strctn;
  int pattern_len=strlen(pattern), s_len;

  int i, j, k;

  flag = 0;
  while(flag!=1)
  {
    flag=fgetline(&s, &bufsize, 32768, FnTRXv, 1);
    flag_strctn = 0;
    s_len=strlen(s);
    for(i=0; i<s_len-pattern_len+1; i++)
    {
      for(j=0; pattern[j] && s[j+i]==pattern[j]; j++);
      if(j==pattern_len)
      {
        flag_strctn=1;
      }
    }
    //printf("s=%s, pattern=%s, flag=%d, flag_strctn=%d, s_len=%d\n", s, pattern, flag, flag_strctn, s_len);
    if(flag_strctn)
    {
      break;
    }
  }

  if(flag==1)
  {
    printf("ERROR: pattern \"%s\" not found the vtxt file\n", pattern);
    exit(0);
  }

  return 0;
}
