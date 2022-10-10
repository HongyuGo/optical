#include "matrix.h"

/* Generate Matrix Struct */
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
/* Copy Mtrix(gen new one)*/
Matrix* Matrix_copy(Matrix* _mat_sourse) {
    Matrix* _mat_copy = Matrix_gen(_mat_sourse->row, _mat_sourse->column, _mat_sourse->data);
    return _mat_copy;
}
/* Free Memory*/
int M_free(Matrix* _mat) {
    free(_mat->data);
    //(_DETAILED_ >= 3) ? printf(">>Matrix_%x has been freed.\n", _mat) : 0;
    free(_mat);
    return 0;
}
/* Add & Sub*/
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
/*Cut_out_part_of_Matrix*/
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
/*Full*/
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
/*Generate Zeros _matrix*/
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
/*Print Matrix*/
int M_print(Matrix* _mat, const char *name) {
    int i, j;
    printf("%s row:%d col:%d\n",name, _mat->row, _mat->column);
    for (i = 0; i < _mat->row; i++) {
        for (j = 0; j < _mat->column; j++) {
            printf(PRECISION, _mat->data[i * (_mat->column) + j]);
        }
        printf("\n");
    }
    return 0;
}
/*Matrix Multiply*/
Matrix* M_numul(Matrix* _mat, double _num) {
    MATRIX_TYPE* data = _mat->data;
    int Size_mat = (_mat->row) * (_mat->column), i;
    for (i = 0; i < Size_mat; i++) {
        data[i] = data[i] * _num;
    }
    return _mat;
}
/*Matrix sub*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num) {
    MATRIX_TYPE* data = _mat->data;
    int Size_mat = (_mat->row) * (_mat->column), i;
    for (i = 0; i < Size_mat; i++) {
        data[i] = data[i] - _num;
    }
    return _mat;
}
/*Generation of transition matrix*/
Matrix* M_Transition(Matrix* _mat) {
    if (_mat == NULL) {
        printf(M_Transition_001);
        return NULL;
    }
    if (_mat->column < 2) {
        printf(M_Transition_002);
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
/*Find min value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Min_value(MATRIX_TYPE *data, int size) {
    MATRIX_TYPE Val_min = data[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        if (data[i] <= Val_min) {
            Val_min = data[i];
        }
    }
    return Val_min;
}
/*Find max value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Max_value(MATRIX_TYPE *data, int size) {
    MATRIX_TYPE Val_max = data[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        if (data[i] >= Val_max) {
            Val_max = data[i];
        }
    }
    return Val_max;
}
/*Assign a value to a position of the matrix*/
void M_value_one(Matrix *_mat, int row , int col,MATRIX_TYPE value){
    int _m_row = _mat->row;
    int _m_column = _mat->column;
    if(row <= 0 || row > _m_row || col <= 0 || col > _m_column){
        printf("%s\n",M_value_one_001);
        return;
    }
    _mat->data[(row - 1) * _m_column + (col - 1)] = value;
}
/*Get a value to a position of the matrix*/
MATRIX_TYPE M_get_one(Matrix *_mat, int row, int col){
    int _m_row = _mat->row;
    int _m_column = _mat->column;
    if(row <= 0 || row > _m_row || col <= 0 || col > _m_column){
        printf("%s\n",M_get_one_001);
        return 0.0;
    }
    return _mat->data[(row - 1) * _m_column + (col - 1)];
}
/*Matrix inverse*/
Matrix *M_Inverse(Matrix *_mat) {
    M_inv_struct *_Uptri_ = M_Uptri_4inv(_mat);
    // M_print(_Uptri_->_matrix,"_Uptri_->_matrix");
    M_inv_struct *_Lowtri_ = M_Lowtri_4inv(_Uptri_->_matrix);
    //M_print(_Lowtri_->_matrix,"_Lowtri_->_matrix");
    Matrix *_mat_dia_inv = M_Dia_Inv(_Lowtri_->_matrix);
    M_print(_mat_dia_inv,"_mat_diq_inv");
    Matrix *_mat_inv = Etrans_4_Inverse(_mat_dia_inv, _Lowtri_->_Etrans_head, _ROW_);
    _mat_inv = Etrans_4_Inverse(_mat_inv, _Uptri_->_Etrans_head, _COLUMN_);
    // 释放内存
    M_free(_Uptri_->_matrix);
    M_free(_Lowtri_->_matrix);
    free(_Uptri_);
    free(_Lowtri_);
    return _mat_inv;
}
/*Upper_triangular_transformation_for_Inverse*/
M_inv_struct *M_Uptri_4inv(Matrix *_mat_source) {
    Matrix *_mat = Matrix_copy(_mat_source);
    int i, j, k, flag;
    Etrans_struct *_Etrans_temp_last = NULL;
    Etrans_struct *_Etrans_temp_head = NULL;

    /*初等变换*/
    for (i = 0; i < _mat->column; i++) {
        for (j = i + 1; j < _mat->row; j++) {
            flag = 0;
            Etrans_struct *_Etrans_temp = (Etrans_struct *) malloc(sizeof(Etrans_struct));
            _Etrans_temp->minuend_line = j + 1;
            _Etrans_temp->subtractor_line = i + 1;
            if ((_mat->data[(_mat->column) * i + i]) != 0) {
                _Etrans_temp->scale = (_mat->data[(_mat->column) * j + i]) / (_mat->data[(_mat->column) * i + i]);
            } else {
                _Etrans_temp->scale = 0;
                for (k = i + 1; k < _mat->row; k++) {
                    flag = 1;//无可替代行
                    if ((_mat->data[(_mat->column) * k + i]) != 0) {
                        _Etrans_temp->minuend_line = -(i + 1);
                        _Etrans_temp->subtractor_line = -(k + 1);
                        flag = 2;//表示能够替换行
                        break;
                    }
                }
                if (flag == 1) {
                    break;
                }
            }
            _Etrans_temp->forward_E_trans = NULL;
            _Etrans_temp->next_E_trans = NULL;
            //if (j==1){
            if (_Etrans_temp_head == NULL) {
                _Etrans_temp_head = _Etrans_temp;
                _Etrans_temp->forward_E_trans = NULL;
            } else {
                _Etrans_temp->forward_E_trans = _Etrans_temp_last;

            }
            if ((i + 1) == _mat->column) {
                _Etrans_temp->next_E_trans = NULL;
            } else {
                if (_Etrans_temp_last != NULL) {
                    _Etrans_temp_last->next_E_trans = _Etrans_temp;
                }
            }
            M_E_trans(_mat, _Etrans_temp, _ROW_);
            //M_print(_mat); //显示具体矩阵
            _Etrans_temp_last = _Etrans_temp;

            if (flag == 2) {
                i = i - 1;
                break;
            }
        }
    }
    M_inv_struct *_Uptri = (M_inv_struct *) malloc(sizeof(M_inv_struct));
    _Uptri->_matrix = _mat;
    _Uptri->_Etrans_head = _Etrans_temp_last;
    return _Uptri;
}
/*_Lower_triangular_transformation_for_Inverse*/
M_inv_struct *M_Lowtri_4inv(Matrix *_mat_source) {
    Matrix *_mat = Matrix_copy(_mat_source);
    int i, j, k, flag;
    Etrans_struct *_Etrans_temp_last = NULL;
    Etrans_struct *_Etrans_temp_head = NULL;
    for (i = 0; i < _mat->row; i++) {
        for (j = i + 1; j < _mat->column; j++) {
            flag = 0;
            Etrans_struct *_Etrans_temp = (Etrans_struct *) malloc(sizeof(Etrans_struct));
            _Etrans_temp->minuend_line = j + 1;
            _Etrans_temp->subtractor_line = i + 1;


            if ((_mat->data[(_mat->column) * i + i]) != 0) {
                _Etrans_temp->scale = (_mat->data[(_mat->column) * i + j]) / (_mat->data[(_mat->column) * i + i]);;
            } else {
                _Etrans_temp->scale = 0;
                for (k = i + 1; k < _mat->row; k++) {
                    flag = 1;//无可替代行
                    if ((_mat->data[(_mat->column) * k + i]) != 0) {
                        _Etrans_temp->minuend_line = -(i + 1);
                        _Etrans_temp->subtractor_line = -(k + 1);
                        flag = 2;//表示能够替换行
                        break;
                    }
                }
                if (flag == 1) {
                    break;
                }
            }

            _Etrans_temp->forward_E_trans = NULL;
            _Etrans_temp->next_E_trans = NULL;
            if (_Etrans_temp_head == NULL) {
                _Etrans_temp_head = _Etrans_temp;
                _Etrans_temp->forward_E_trans = NULL;
            } else {
                _Etrans_temp->forward_E_trans = _Etrans_temp_last;
            }
            if ((i + 1) == _mat->column) {
                _Etrans_temp->next_E_trans = NULL;
            } else {
                if (_Etrans_temp_last != NULL) {
                    _Etrans_temp_last->next_E_trans = _Etrans_temp;
                }
            }
            M_E_trans(_mat, _Etrans_temp, _COLUMN_);
            //M_print(_mat); //显示具体矩阵
            _Etrans_temp_last = _Etrans_temp;
            if (flag == 2) {
                i = i - 1;
                break;
            }
        }
    }
    M_inv_struct *_Lowtri = (M_inv_struct *) malloc(sizeof(M_inv_struct));
    _Lowtri->_matrix = _mat;
    _Lowtri->_Etrans_head = _Etrans_temp_last;
    return _Lowtri;
}
/*M_Inv for Dia_matrix*/
Matrix *M_Dia_Inv(Matrix *_mat_source) {
    Matrix *_mat_inv = NULL;
    if (_mat_source->row != _mat_source->column) {
        // printf(M_Dia_Inv_002);
        // system("pause");
    } else {
        _mat_inv = Matrix_copy(_mat_source);
        MATRIX_TYPE *data = _mat_inv->data;
        int i, order = _mat_source->column;
        for (i = 0; i < order; i++) {
            if((data)[i * (order + 1)] == 0){ // 不可逆
                // printf(M_Dia_Inv_023);
                // system("pause");
                // (data)[i * (order + 1)] = 1 / (data[i * (order + 1)]);
                (data)[i * (order + 1)]  = 0.0;
            }else{
                (data)[i * (order + 1)] = 1 / (data[i * (order + 1)]);
            }
        }
    }
    return _mat_inv;
}
/*Inverse_Element_trans_to_Matrix*/
Matrix *Etrans_4_Inverse(Matrix *_mat_result, Etrans_struct *_Etrans_, int line_setting) {
    Etrans_struct *temp_Etrans = _Etrans_, *temp_Etrans_pre = _Etrans_;
    int temp_num = 0;
    // 此处方案感谢 @1u2e, github.com/Amoiensis/Matrix_hub/issues/4
    while (temp_Etrans != NULL) {
        temp_num = temp_Etrans->minuend_line;
        temp_Etrans->minuend_line = temp_Etrans->subtractor_line;
        temp_Etrans->subtractor_line = temp_num;
        M_E_trans(_mat_result, temp_Etrans, line_setting);
        // 此处修改方案感谢 @1u2e, github.com/Amoiensis/Matrix_hub/issues/4
        temp_Etrans = temp_Etrans->forward_E_trans;
        free(temp_Etrans_pre);
        temp_Etrans_pre = temp_Etrans;
    }
    return _mat_result;
}
/*Element teransfor Matrix*/
/*lin3_sstting 设置是行初等变换还是列初等变换*/
int M_E_trans(Matrix *_mat, Etrans_struct *_Etrans_, int line_setting) {
    int line_num, i;
    if (line_setting == _ROW_) {
        /*行初等变换*/
        line_num = _mat->column;
        if (_Etrans_->scale) {
            for (i = 0; i < line_num; i++) {
                _mat->data[(_Etrans_->minuend_line - 1) * (_mat->column) + i] -=
                        (_Etrans_->scale) * (_mat->data[(_Etrans_->subtractor_line - 1) * (_mat->column) + i]);

            }
        } else {
            if ((_Etrans_->minuend_line < 0) && (_Etrans_->subtractor_line < 0)) {/*交换*/
                M_Swap(_mat, -(_Etrans_->minuend_line), -(_Etrans_->subtractor_line), line_setting);
            }
        }
    } else {
        /*列初等变换*/
        line_num = _mat->row;
        if (_Etrans_->scale) {
            for (i = 0; i < line_num; i++) {
                _mat->data[(_Etrans_->minuend_line - 1) + (_mat->column) * i] -=
                        (_Etrans_->scale) * (_mat->data[(_Etrans_->subtractor_line - 1) + (_mat->column) * i]);
            }
        } else {
            if ((_Etrans_->minuend_line < 0) && (_Etrans_->subtractor_line < 0)) {/*交换*/
                M_Swap(_mat, -(_Etrans_->minuend_line), -(_Etrans_->subtractor_line), line_setting);
            }
        }
    }
    return 0;
}
/*Swap Line*/
int M_Swap(Matrix *_mat, int _line_1, int _line_2, int line_setting) {
    _line_1 = _line_1 - 1;
    _line_2 = _line_2 - 1;
    int i;
    MATRIX_TYPE temp;
    if (line_setting == _ROW_) {
        if ((_line_1 < _mat->row) && (_line_2 < _mat->row)) {
            for (i = 0; i < (_mat->column); i++) {
                temp = _mat->data[_line_1 * (_mat->column) + i];
                _mat->data[_line_1 * (_mat->column) + i] = _mat->data[_line_2 * (_mat->column) + i];
                _mat->data[_line_2 * (_mat->column) + i] = temp;
            }
        } else {
            printf(M_swap_004);
            system("pause");
        }
    } else {
        if ((_line_1 < _mat->column) && (_line_2 < _mat->column)) {
            for (i = 0; i < (_mat->row); i++) {
                temp = _mat->data[_line_1 + (_mat->column) * i];
                _mat->data[_line_1 + (_mat->column) * i] = _mat->data[_line_2 + (_mat->column) * i];
                _mat->data[_line_2 + (_mat->column) * i] = temp;
            }
        } else {
            printf(M_swap_004);
            system("pause");
        }
    }
    return 0;
}

