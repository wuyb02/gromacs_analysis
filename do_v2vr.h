/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_V2VR_H__
#define __DO_V2VR_H__

#include "do_xxx.h"
#include "pdbio.h"
#include "statutil.h"

void do_v2vr(char *fnTRX, char *fnTRXo,
	     int snp, real dT, real Xcom, real Ycom);

int v2vr(t_trxframe *fr_p, t_trxframe *frout_p, real Xcom, real Ycom);


#endif
