#include <stdio.h>
#include <math.h>
#include "math_base.h"
#include "math_linalg.h"
//if nn={8,8}, ndim=2 then data[128] for 64 complex #s
EXP void fourn(double*, unsigned long *nn, int ndim, int isign);
//speq(nx,2*ny)
EXP void rlft3(double*, math::Matrix<double>& speq, unsigned long, unsigned long, unsigned long, int);
