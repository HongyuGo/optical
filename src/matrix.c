#include "matrix.h"

/* Generate Matrix Struct 导入_生成矩阵
 */
Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data) {
    Matrix* _mat = (Matrix*)malloc(sizeof(Matrix));
    if (_mat == NULL) return 0;
    _mat->row = row;
    _mat->column = column;
    int size = _mat->row * _mat->column;
    _mat->data = (MATRIX_TYPE*)malloc((size) * sizeof(MATRIX_TYPE));
    int i;
    for (i = 0; i < size; i++) {
        _mat->data[i] = data[i];
    }
    return _mat;
}

/* Copy Mtrix(gen new one)
         复制矩阵（生成新矩阵）*/
Matrix* Matrix_copy(Matrix* _mat_sourse) {
    Matrix* _mat_copy = Matrix_gen(_mat_sourse->row, _mat_sourse->column, _mat_sourse->data);
    return _mat_copy;
}

/* Free Memory
   释放矩阵，释放内存
*/
int M_free(Matrix* _mat) {
    free(_mat->data);
    //(_DETAILED_ >= 3) ? printf(">>Matrix_%x has been freed.\n", _mat) : 0;
    free(_mat);
    return 0;
}

/* Add & Sub
矩阵加减法
*/
Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus, Matrix* _mat_minus) {
    Matrix* _mat_result = NULL;
    if ((_mat_subed->column == _mat_minus->column) && (_mat_subed->row == _mat_minus->row)) {
        _mat_result = Matrix_copy(_mat_subed);
        int size = (_mat_subed->row) * (_mat_subed->column), i;
        for (i = 0; i < size; i++) {
            _mat_result->data[i] = (_mat_result->data[i]) * scale_mat_subed - (_mat_minus->data[i]) * scale_mat_minus;
        }
    } else {
        printf(M_add_sub_003);
    }
    return _mat_result;
}

/*Cut_out_part_of_Matrix
切取部分矩阵*/
Matrix* M_Cut(Matrix* _mat, int row_head, int row_tail, int column_head, int column_tail) {
    Matrix* mat_result = NULL;
    if (row_tail < 0) {
        if (row_tail == _END_) {
            row_tail = _mat->row;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (row_head < 0) {
        if (row_head == _END_) {
            row_head = _mat->row;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (column_tail < 0) {
        if (column_tail == _END_) {
            column_tail = _mat->column;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (column_head < 0) {
        if (column_head == _END_) {
            column_head = _mat->column;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if ((row_tail > _mat->row) || (column_tail > _mat->column)) {
        printf(M_Cut_005);
        system("pause");
    } else {
        if ((row_head > row_tail) || (column_head > column_tail)) {
            printf(M_Cut_006);
            system("pause");
        } else {
            row_head = row_head - 1;
            column_head = column_head - 1;
            mat_result = (Matrix*)malloc(sizeof(Matrix));
            mat_result->row = row_tail - row_head;
            mat_result->column = column_tail - column_head;
            mat_result->data = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * (mat_result->row) * (mat_result->column));
            int i, j;
            for (i = 0; i < (row_tail - row_head); i++) {
                for (j = 0; j < (column_tail - column_head); j++) {
                    mat_result->data[i * (mat_result->column) + j] =
                        _mat->data[(i + row_head) * (_mat->column) + (j + column_head)];
                }
            }
        }
    }
    return mat_result;
}

/*Full
填充矩阵*/
Matrix* M_full(Matrix* _mat, int row_up, int row_down, int column_left, int column_right, MATRIX_TYPE full_data) {
    Matrix* mat_result = NULL;
    mat_result = (Matrix*)malloc(sizeof(Matrix));
    mat_result->row = (_mat->row + row_up + row_down);
    mat_result->column = (_mat->column + column_left + column_right);
    mat_result->data = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * (mat_result->row) * (mat_result->column));
    int i, j;
    for (i = 0; i < mat_result->row; i++) {
        for (j = 0; j < mat_result->column; j++) {
            if ((i >= row_up) && (i < (row_up + _mat->row))) { /*这里的双判断，可以优化*/
                if ((j >= column_left) && (j < (column_left + _mat->column))) {
                    mat_result->data[i * (mat_result->column) + j] =
                        _mat->data[(_mat->column) * (i - row_up) + (j - column_left)];
                } else {
                    mat_result->data[i * (mat_result->column) + j] = full_data;
                }
            } else {
                mat_result->data[i * (mat_result->column) + j] = full_data;
            }
        }
    }
    return mat_result;
}
/*Generate Zeros _matrix
生成全零矩阵*/
Matrix* M_Zeros(int row, int column) {
    Matrix* Zero_mat = (Matrix*)malloc(sizeof(Matrix));
    Zero_mat->column = column;
    Zero_mat->row = row;
    int size = row * column, i;
    MATRIX_TYPE* data = (MATRIX_TYPE*)malloc((size) * sizeof(MATRIX_TYPE));
    for (i = 0; i < size; i++) {
        data[i] = 0;
    }
    Zero_mat->data = data;
    return Zero_mat;
}
/*Print Matrix
打印矩阵*/
int M_print(Matrix* _mat) {
    int i, j;
    printf("row:%d col:%d\n", _mat->row, _mat->column);
    for (i = 0; i < _mat->row; i++) {
        for (j = 0; j < _mat->column; j++) {
            printf(PRECISION, (int)_mat->data[i * (_mat->column) + j]);
        }
        printf("\n");
    }
    return 0;
}
/*Matrix Multiply
矩阵数乘*/
Matrix* M_numul(Matrix* _mat, double _num) {
    MATRIX_TYPE* data = _mat->data;
    int Size_mat = (_mat->row) * (_mat->column), i;
    for (i = 0; i < Size_mat; i++) {
        data[i] = data[i] * _num;
    }
    return _mat;
}
/*矩阵数减*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num) {
    MATRIX_TYPE* data = _mat->data;
    int Size_mat = (_mat->row) * (_mat->column), i;
    for (i = 0; i < Size_mat; i++) {
        data[i] = data[i] - _num;
    }
    return _mat;
}
/*跃迁矩阵生成*/
Matrix* Matrix_Transition(Matrix* _mat) {
    if (_mat == NULL) {
        printf(Matrix_Transition_001);
        return NULL;
    }
    if (_mat->column < 2) {
        printf(Matrix_Transition_002);
        return NULL;
    }
    Matrix* _mat_result = NULL;
    int row = _mat->row;
    int col = _mat->column;
    int size_data = row * col;
    MATRIX_TYPE* _data = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * (row * col - 1));
    for (int i = 0; i < size_data - 1; i++) {
        _data[i] = _mat->data[i + 1] - _mat->data[i];
    }
    _mat_result = (Matrix*)malloc(sizeof(Matrix));
    _mat_result->row = row;
    _mat_result->column = col - 1;
    _mat_result->data = _data;
    return _mat_result;
}
