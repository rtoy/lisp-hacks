
/* @(#)fdlibm.h 1.5 04/04/22 */
/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 */

/* Sometimes it's necessary to define __LITTLE_ENDIAN explicitly
   but these catch some common cases. */

#if defined(i386) || defined(i486) || \
	defined(intel) || defined(x86) || defined(i86pc) || \
	defined(__alpha) || defined(__osf__)
#define __LITTLE_ENDIAN
#endif

#ifdef __LITTLE_ENDIAN
enum { HIWORD = 1, LOWORD = 0 };
#else
enum { HIWORD = 0, LOWORD = 1 };
#endif

/*
 * ANSI/POSIX
 */
extern double fabs(double);
extern double floor(double);

/* ieee style elementary functions */
extern int    __ieee754_rem_pio2(double,double*);

/*
 * Functions callable from C, intended to support IEEE arithmetic.
 */
extern double scalbn(double, int);

/* fdlibm kernel function */
extern int    __kernel_rem_pio2(double*,double*,int,int,int,const int*);
