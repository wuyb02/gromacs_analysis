/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_TXT2TRR_H__
#define __DO_TXT2TRR_H__

#include "do_xxx.h"
#include "pdbio.h"
#include "statutil.h"

void do_txt2trr_cnt(char *fnTRX, char *fnTRXx, char *fnTRXv, char *fnTRXo,
	           int snp, real dT, bool bVels);
int read_next_txtframe(FILE *FnTRXx, rvec *xmem, int nskip, int natom);


#endif
