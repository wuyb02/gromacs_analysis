#ifndef __MGROMACS_H__
#define __MGROMACS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "pdbio.h"
#include "index.h"
#include "smalloc.h"
//#include "fftgrid.h"
#include "calcgrid.h"
#include "nrnb.h"
//#include "shift_util.h"
//#include "pme.h"
#include "gstat.h"
#include "matio.h"

#define nint(a) ((int)(999.5+a)-999)

#endif
