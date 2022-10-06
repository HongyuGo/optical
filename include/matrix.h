#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include "commom.h"

typedef struct {
    int row, column;
    MATRIX_TYPE * data;
} Matrix;

/* Generate Matrix Struct */
Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data);
/* Copy Mtrix(gen new one)*/
Matrix* Matrix_copy(Matrix* _mat_sourse);
/* Free Memory*/
int M_free(Matrix* _mat);
/* Add & Sub*/
Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus, Matrix* _mat_minus);
/*Cut_out_part_of_Matrix*/
Matrix* M_Cut(Matrix* _mat, int row_head, int row_tail, int column_head, int column_tail);
/*FULL*/
Matrix* M_full(Matrix* _mat, int row_up, int row_down, int column_left, int column_right, MATRIX_TYPE full_data);
/*Generate Zeros _matrix*/
Matrix* M_Zeros(int row, int column);
/*Print Matrix*/
int M_print(Matrix* _mat, const char *name);
/*Matrix Multiply*/
Matrix* M_numul(Matrix* _mat, double _num);
/*Matrix sub*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num);
/*Generation of transition matrix*/
Matrix* M_Transition(Matrix* _mat);
/*Find min value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Min_value(MATRIX_TYPE *data, int size); 
/*Find max value in a MATRIX_TYPE*/
MATRIX_TYPE M_Max_value(MATRIX_TYPE *data, int size);

#endif
