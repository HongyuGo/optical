#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <math.h>
#include "matrix.h"

double h_response(double, double, double);
double readback(double ,double , Matrix * ,double ,double ,double);
void gen_firtaps_v2(Matrix* random, Matrix *sampled, MATRIX_TYPE* gpr_coeff, int fir_len, const char * constraint,char method);
#endif
