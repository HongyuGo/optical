#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "params.h"

typedef struct {
    int row, column;
    MATRIX_TYPE** data;
} Matrix;
/*Store the result of matrix_inverse*/
typedef struct _matrix_inverse_struct {
    Matrix* _matrix;
    struct _Elementary_Transformation* _Etrans_head;
} M_inv_struct;

/*Store the Operation of Elementary_Transformation*/
typedef struct _Elementary_Transformation {
    int minuend_line;
    int subtractor_line;
    TRANS_TYPE scale;
    struct _Elementary_Transformation* forward_E_trans;
    struct _Elementary_Transformation* next_E_trans;
} Etrans_struct;

/*Store the result of matrix_eigen*/
typedef struct _matrix_eigen_struct_single {
    Matrix* eigen_matrix;
    double eigen_value;
} M_eigen_struct;

/*trellis*/
typedef struct {
    int input;
    MATRIX_TYPE output;
    // int cur;
    int next;
    // int counter; /*we can default that every state has two branch*/
} Trellis_nst;
typedef struct{
    int input[2];
    MATRIX_TYPE output[2];
    //int next;
    int pre[2];
    int counter;
}Trellis_pst;

MATRIX_TYPE** GetMemory(int row, int col);
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
int M_print(Matrix* _mat, const char* name);
/*Matrix Multiply*/
Matrix* M_nummul(Matrix* _mat, double _num);
/*Matrix sub*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num);
/*Generation of transition matrix*/
Matrix* M_Transition(Matrix* _mat);
/*Find min value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Min_value(MATRIX_TYPE* data, int size);
/*Find max value in a MATRIX_TYPE*/
MATRIX_TYPE M_Max_value(MATRIX_TYPE* data, int size);
/*Assign a value to a position of the matrix*/
void M_value_one(Matrix* _mat, int row, int col, MATRIX_TYPE value);
/*Get a value to a position of the matrix*/
MATRIX_TYPE M_get_one(Matrix* _mat, int row, int col);
/*Matrix inverse*/
Matrix* M_Inverse(Matrix* _mat);
/*Matrix Multiply*/
Matrix* M_mul(Matrix* _mat_left, Matrix* _mat_right);
/*Transpose*/
Matrix* M_T(Matrix* _mat_source);
Matrix* M_limit(Matrix* _mat);
/*Matrix Convolution*/
Matrix* M_Conv(Matrix* _mat1, Matrix* _mat2);
/*Matrix Reverse*/
Matrix* M_Rever(Matrix* _mat, int row);
/*Matrix Max Eigenvalue(vec)*/
M_eigen_struct* M_eigen_max(Matrix* _mat);
/*Generate Ones _matrix*/
Matrix* M_Ones(int row, int column);
/*compare the diff of two Matrix*/
int M_Compare(Matrix *mat1,Matrix *mat2);
/*viterbi_mlse*/
Matrix* viterbi_mlse(int gpr_length, Matrix* fk1, Matrix* gpr_coeff);
#endif
