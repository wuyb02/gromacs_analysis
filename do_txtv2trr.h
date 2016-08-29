/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_TXTV2TRR_H__
#define __DO_TXTV2TRR_H__

#include "do_xxx.h"
#include "pdbio.h"
#include "statutil.h"

void do_txtv2trr(char *fnTRX, char *fnTRXv, char *fnTRXo,
	         int snp, real dT);
int read_next_txtvframe(FILE *fntrxv, rvec *vmem);
int fread_patternskip(FILE *FnTRXv, char *pattern);

#endif
