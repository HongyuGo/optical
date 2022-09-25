#include "matrix.h"

Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data) { /*Generate Matrix Struct
        导入_生成矩阵*/
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

Matrix* Matrix_copy(Matrix* _mat_sourse) { /*Copy Mtrix(gen new one)
         复制矩阵（生成新矩阵）*/
    Matrix* _mat_copy = Matrix_gen(_mat_sourse->row, _mat_sourse->column, _mat_sourse->data);
    return _mat_copy;
}

int M_free(Matrix* _mat) { /*Free Memory
         释放矩阵，释放内存*/
    free(_mat->data);
    //(_DETAILED_ >= 3) ? printf(">>Matrix_%x has been freed.\n", _mat) : 0;
    free(_mat);
    return 0;
}

Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus,
                  Matrix* _mat_minus) { /*Add & Sub
矩阵加减法*/
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
