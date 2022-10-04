#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include "commom.h"

typedef struct {
    int row, column;
    MATRIX_TYPE * data;
} Matrix;
/*导入_生成矩阵*/
Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data);
/*复制矩阵（生成新矩阵）*/
Matrix* Matrix_copy(Matrix* _mat_sourse);
/*释放矩阵，释放内存*/
int M_free(Matrix* _mat);
/*矩阵加减法*/
Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus, Matrix* _mat_minus);
/*切取部分矩阵*/
Matrix* M_Cut(Matrix* _mat, int row_head, int row_tail, int column_head, int column_tail);
/*填充矩阵*/
Matrix* M_full(Matrix* _mat, int row_up, int row_down, int column_left, int column_right, MATRIX_TYPE full_data);
/*生成全零矩阵*/
Matrix* M_Zeros(int row, int column);
/*打印矩阵*/
int M_print(Matrix* _mat, const char *name);
/*矩阵数乘*/
Matrix* M_numul(Matrix* _mat, double _num);
/*矩阵数减*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num);
/*跃迁矩阵生成*/
Matrix* M_Transition(Matrix* _mat);
/*找到行矩阵中最小的值*/
MATRIX_TYPE M_Min_value(MATRIX_TYPE *data, int size); 
/*找到行矩阵中最大的值*/
MATRIX_TYPE M_Max_value(MATRIX_TYPE *data, int size);

#endif
