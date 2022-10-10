#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include "commom.h"

typedef struct {
    int row, column;
    MATRIX_TYPE * data;
} Matrix;
/*Store the result of matrix_inverse*/
typedef struct _matrix_inverse_struct {
    Matrix *_matrix;
    struct _Elementary_Transformation *_Etrans_head;
} M_inv_struct;

/*Store the Operation of Elementary_Transformation*/
typedef struct _Elementary_Transformation {
    
    int minuend_line;
    int subtractor_line;
    TRANS_TYPE scale;
    struct _Elementary_Transformation *forward_E_trans;
    struct _Elementary_Transformation *next_E_trans;
} Etrans_struct;

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
/*Assign a value to a position of the matrix*/
void M_value_one(Matrix *_mat, int row , int col, MATRIX_TYPE value);
/*Get a value to a position of the matrix*/
MATRIX_TYPE M_get_one(Matrix *_mat, int row, int col);
/*Matrix inverse*/
Matrix *M_Inverse(Matrix *_mat);
/*Upper_triangular_transformation_for_Inverse*/
M_inv_struct *M_Uptri_4inv(Matrix *_mat_source);
/*Element teransfor Matrix*/
/*lin3_sstting 设置是行初等变换还是列初等变换*/
int M_E_trans(Matrix *_mat, Etrans_struct *_Etrans_, int line_setting); 
/*Swap Line*/
int M_Swap(Matrix *_mat, int _line_1, int _line_2, int line_setting);
/*_Lower_triangular_transformation_for_Inverse*/
M_inv_struct *M_Lowtri_4inv(Matrix *_mat_source); 
/*M_Inv for Dia_matrix*/
Matrix *M_Dia_Inv(Matrix *_mat_source); 
/*Inverse_Element_trans_to_Matrix*/
Matrix *Etrans_4_Inverse(Matrix *_mat_result, Etrans_struct *_Etrans_, int line_setting); 

#endif
