#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include "commom.h"

typedef struct {
    int row, column;
    int* data;
} Matrix;

Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data);
Matrix* Matrix_copy(Matrix* _mat_sourse);
int M_free(Matrix* _mat);
Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus, Matrix* _mat_minus);

#endif
