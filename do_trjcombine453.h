/**
 ** do_wtcomxy.h
 ** Original: GROMACS
 ** Modified by Rui Qiao, 05/10/2002, 01/01/2003
 ** Modified by Yanbin, 11/14/2010, increasing statistics by 2
 **
 **/

#ifndef __DO_TRJCOMBINE_H__
#define __DO_TRJCOMBINE_H__

#include "do_xxx453.h"
#include "pdbio.h"
#include "statutil.h"

void do_trjcombine(char *fnTRX1, char *fnTRX2, char *fnTRXo,
	           int snp, real dT);


#endif
