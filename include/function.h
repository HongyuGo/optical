#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <math.h>
#include "matrix.h"

double h_response(double, double, double);
double readback(double ,double , Matrix * ,double ,double ,double);
void gen_firtaps_v2(Matrix* random, Matrix *sample, MATRIX_TYPE* gpr_coeff, int fir_len, char constraint,char method);
Matrix* auto_corr(Matrix *x, int bot, int top);
//void auto_corr(Matrix * A, Matrix * B,int shift , Matrix * output, int row, int col);
#endif
