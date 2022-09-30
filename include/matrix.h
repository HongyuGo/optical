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
Matrix* M_Cut(Matrix* _mat, int row_head, int row_tail, int column_head, int column_tail);
Matrix* M_full(Matrix* _mat, int row_up, int row_down, int column_left, int column_right, MATRIX_TYPE full_data);
Matrix* M_Zeros(int row, int column);
int M_print(Matrix* _mat);
Matrix* M_numul(Matrix* _mat, void* _num, char type);
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num);
Matrix* Matrix_Transition(Matrix* _mat);
#endif
