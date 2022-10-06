#include "function.h"
#include <math.h>
#include "matrix.h"

double h_response(double t, double TL, double S) { return erf(t / S / TL); }

double readback(double t, double jitter, Matrix* d, double S, double T, double TL) {
    int len = d->column;
    double rs = 0;
    double tmp = 0;
    for (int i = 0; i < len; i++) {
        tmp = d->data[i] * h_response(t - (i + 1) * T + jitter, TL, S);
        rs = rs + tmp;
    }
    return rs;
}

void gen_firtaps_v2(Matrix *random, Matrix *ampled, MATRIX_TYPE* gpr_coeff, int fir_len, const char * constraint,char method)
{

}
