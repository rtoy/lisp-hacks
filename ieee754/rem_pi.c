/* Simplified interface to e_rem_pio2 for Lisp */

#include "fdlibm.h"

int rem_pio2 (double x, double *y0, double *y1)
{
    double y[2];
    int n;

    n = ieee754_rem_pio2(x, y);
    *y0 = y[0];
    *y1 = y[1];
    return n;
}
